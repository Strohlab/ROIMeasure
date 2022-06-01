%% ROI time measurement and dF/F tool
% Version 4.1
% Stroh-Lab
% Zeke Barger, Hirofumi Watari, Saleh Altahini

%% Preset
% Fill the following parameters to skip the question dialogs
% or leave empty to show the prompts

pathName='F:\2 Round\';    % keep track of selected directory

ftype=2; %1=TIFF stack  2=image sequence
img1=1; %1=choose existing image     2=create average image
imgWidthMicron=458; %image width in microns
hz=30.8; %imaging frequency
savePlots = 1; % change to 0 to not save the plots and figures

%% Data Loading


rois={};        % empty table to load ROIs into
chosenRoi=0;
loadrois=[];
choice=questdlg('How do you want to start?','Region of interest',...
    'Draw new ROIs','Load saved ROIs','Cancel','Draw new ROIs');
switch choice
    case 'Draw new ROIs'
        loadrois=0;
    case 'Load saved ROIs'
        loadrois=1;
    case 'Cancel'
        disp('User canceled');
        return
end

% get saved rois
if loadrois 
    disp('Select your saved ROI .mat file');
    isDone=0;
    while ~isDone
        [fname1 path1]=uigetfile(strcat(pathName,'\*.mat'),'Select the MATLAB file containing the ROIs');
        if ~fname1
            disp('User canceled');
            return
        end
        pathName=path1;    % save for later
        %load(cat(2,path1, fname1));
        load(strcat(path1,fname1));
        if exist('rois','var')
            isDone=1;
            numrois=size(rois,2);
        else
            uiwait(warndlg('ROI is not in this file. Try again','Wrong MATLAB file'));
        end
    end
end


% load files, except not

disp(' ');
if ftype==1
    disp('Select your TIFF stack');
    [fname,pname]=uigetfile(strcat(pathName,'*.*'),'Your TIFF stack!');
    info = imfinfo(cat(2,pname,fname));
    num_images = numel(info);
    first = imread(cat(2,pname,fname), 1, 'Info', info);
    size_first=size(first);  %get size of first image to preallocate space for imported images
    
    baseFileName=fname;
    pathName=pname;

else
    % Open image sequence using gui 
    disp('Select your folder containing a sequence of images');
    isDone=0;
    while ~isDone
        mydir=uigetdir(pathName,'Select a folder containing a sequence of images (non-RGB)');            %use graphical user interface to set file path in mydir
        if ~mydir
            disp('User canceled');
            return
        end
        dirlist=dir(strcat(mydir,filesep,'*.tif'));                    %save file path in structure array
        if ~isempty(dirlist)
            isDone=1;
        else
            uiwait(warndlg('No Tiff images found in this folder. Try again','Wrong folder'));
        end
    end
    first=imread(strcat(mydir,filesep,dirlist(1).name));           %pull out first image from stack
    size_first=size(first);                                        %get size of first image to preallocate space for imported images
    num_images=length(dirlist);
    
    % Save imgWidth and imgHeight
    [imgHeight,imgWidth]=size(first);

    % extract folder name (to be used as a base file name later)
    [pathName,baseFileName,ext]=fileparts(mydir);
    

end

if img1==0
    choice=questdlg('what do you want to use?','Average Image',...
    'Create average image','Choose existing image','Cancel','Create average image');
    switch choice
        case 'Create average image'
            img1=2;
        case 'Choose existing image'
            img1=1;
        case 'Cancel'
            disp('User canceled');
            return
    end
end

% Average Image Dialog
if img1==2
    % Create average intensity using all images
    imgFrom=1;
    imgTo=num_images;
    isDone=0;
    while ~isDone
        prompt = {'From','To'};
        dlg_title = 'Average images';
        num_lines = 1;
        def = {num2str(imgFrom),num2str(imgTo)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if ~isempty(answer)
            if str2double(answer{1})<str2double(answer{2})
                imgFrom=str2double(answer{1});
                imgTo=str2double(answer{2});
            else
                imgFrom=str2double(answer{2});
                imgTo=str2double(answer{1});
            end
            % Reject out-of-range errors
            if imgFrom>=1 && imgTo<=num_images
                isDone=1;
                num2avg = imgTo - imgFrom + 1;
                disp(['Averaging ',num2str(num2avg),' images (from ',num2str(imgFrom),' to ',num2str(imgTo),')']);
            end
        else

            disp('User canceled');
            return
        end
    end

    h60 = waitbar(0,'Creating average image');
    [i,j]=size(first);
    newimg=zeros(i,j);
    sumimg=zeros(i,j,'uint64');
    if ftype==1
        for ii=imgFrom:imgTo
            newimg=imread(cat(2,pname,fname),ii,'Info',info);
            sumimg=sumimg+uint64(newimg);
            waitbar(ii/num2avg,h60)
        end
    else
        for ii=imgFrom:imgTo
            newimg=imread(strcat(mydir,filesep,dirlist(ii).name));
            sumimg=sumimg+uint64(newimg);
            waitbar(ii/num2avg,h60)
        end
    end
    
    first2=zeros(i,j);
    first2=sumimg/num2avg;
    
    % check if the summation reached a ceiling
    maxSumImg=max(max(sumimg));
    if maxSumImg>=2^64
        uiwait(errordlg('Too many images were loaded. Try again.','Average Error'));
        disp('Force quit: Error due to averaging too many images');
        close(h60)
        return
    end
    close(h60)
    % Ends Hiro's import method
    

else
    % ask for exsiting image, making sure it has the same dimensions
    first2=[0 0];
    disp(' '); disp('Select the average image. The dimensions must match the film');
    while [size(first,1) size(first,2)] ~= [size(first2,1) size(first2,2)]  
        [fname2,pname2]=uigetfile(strcat(pathName,'/*.*'),'Brightfield/average image');
        info2 = imfinfo(cat(2,pname2,fname2));
        num_images2 = numel(info2);
        first2 = imread(cat(2,pname2,fname2), 1, 'Info', info2);
        if [size(first,1) size(first,2)] ~= [size(first2,1) size(first2,2)]
           disp('Your dimensions did not match, try again'); 
        end
    end
end


% ask the width in micrometer
if imgWidthMicron==0
    prompt = {['Image width is ',num2str(imgWidth),' pixels. What is the width in micrometer?']};
    dlg_title = 'Image scale';
    num_lines = 1;
    def = {num2str(imgWidthMicron)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        imgWidthMicron=str2double(answer{1});
    else
        disp('User canceled');
        return
    end
end


% Prompt for sampling frequency
if hz==0
    prompt = {'Frame rate (Hz)'};
    dlg_title = 'Enter frame rate';
    num_lines = 1;
    def = {num2str(hz)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        hz=str2double(answer{1});
        % Reject out-of-range errors
        if hz>0
            isDone=1;
        end
    else
        disp('User canceled');
        return
    end
end


%% Load Saved ROIs
intensity=[];
pixels=[];
mask2={};
maskedvec=[];

if loadrois 

    cellfig=figure;
    colormap(gray);
    imagesc(first2);
    axis image;
    if chosenRoi==0
        prompt = {'Alignment ROI'};
        dlg_title = 'Which ROI to align?';
        num_lines = 1;
        def = {num2str(chosenRoi)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if ~isempty(answer)
            chosenRoi=str2double(answer{1});
            % Reject out-of-range errors
            if chosenRoi>0
                isDone=1;
            end
        else
            disp('User canceled');
            return
        end
    end
    firstROI=images.roi.Freehand(gca,'Position',rois{chosenRoi}.Position);
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gcf,'Color',[0 0 0]);
    xlabel({'Optionally move ROI #1 to a new location.';'Double-click the ROI when done.'});
    ddd=get(gcf,'Position'); % make the figure window play nice
    wait(firstROI);
    if firstROI.Position ~= rois{chosenRoi}.Position
        correct_drift = true;
        firstP = rois{chosenRoi}.Position;
        firstP2=firstROI.Position;
        drift=[firstP2(1,1)-firstP(1,1) firstP2(1,2)-firstP(1,2)]; % change in position
    else
        correct_drift = false;
    end
    delete(firstROI);
    
    if correct_drift
        choice=questdlg('what do you want to use?','Drift Correction',...
    'Manually place each ROI','Automatic drift correction','Cancel','Manually place each ROI');
    end
    
    hh2 = waitbar(0,'Loading ROIs...');

    for i=1:numrois
        if correct_drift
            switch choice
                case 'Manually place each ROI'
                    xlabel({'Move the ROI to a new location.';'Double-click the ROI to continue.'});
                    rois{i} = images.roi.Freehand(gca,'Position',rois{i}.Position);
                    rois{i}.Label = num2str(i);
                    wait(rois{i});
                    rois{i}.Waypoints = false(size(rois{i}.Waypoints));
                    rois{i}.Color = 'w';
                case 'Automatic drift correction'
                    oldROI = rois;
                    rois{i} = images.roi.Freehand(gca,'Position',(rois{i}.Position+drift));
                    rois{i}.Label = num2str(i);
                    rois{i}.Waypoints = false(size(rois{i}.Waypoints));
                    rois{i}.Color = 'w';
                case 'Cancel'
                    disp('User canceled');
                return
            end
        else
            %rois{i}.Parent= gca;
            oldROI = rois;
            rois{i} = images.roi.Freehand(gca,'Position',(rois{i}.Position));
            rois{i}.Label = num2str(i);
            rois{i}.Waypoints = false(size(rois{i}.Waypoints));
            rois{i}.Color = 'w';
        end
        waitbar(i/numrois,hh2)
       
    end
    close(hh2) 
   

    pause(.1)
end

%% mew rois

thisisok = 0;
while thisisok == 0 % get rois until user is satisfied with them
    % measure intensity in rois       
    if ~loadrois
        cellfig=figure;
        colormap(gray);
        imagesc(first2);
        axis image;
        set(gcf,'Color',[0 0 0]);
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'YTick',[]);
        xlabel('Press Esc to undo. To finish, Shift-click inside the image, then click outside it');
        zoom reset    
        ddd=get(gcf,'Position'); % make the figure window play nice
        
        
        i=1;
        
        set(gcf,'toolbar','figure');
        % while i < (numrois+1)
        uiwait(msgbox({'Click and drag to draw the ROIs.' ...
        'Hit Esc to undo.' ...
        'When done, Shift-click inside the image.' ...
        'and then click outside the image'},'Instructions'));
    else
        xlabel('Press Esc to undo. To finish click outside of the image.');
        i = numrois+1;
    end
    maxX=get(gca,'XLim');
    maxY=get(gca,'YLim');
    quitcount=0;
    while quitcount==0
        numrois=size(rois,2);
        rois{i}=drawfreehand('Color', 'w');
        if isempty(rois{i}.Position)
            if i>1
                i=i-1;
                delete(rois{i});
                rois(end) = [];
            end
        else
            if size(rois{i}.Position) == [1,2] & ~isempty(intersect(rois{i}.Position, [maxX, maxY]))
                delete(rois{i});
                rois(end) = [];
                choice=questdlg('Done?','Tell me','Yes','No','Yes');
                switch choice
                case 'Yes'
                    quitcount=1;
                    numrois=size(rois,2);
                case 'No'
                    quitcount=0;
                case ''
                    quitcount=[];
                    disp('User canceled');
                    return
                end
                
            else
                disp(size(rois{i}.Position));
                disp(intersect(rois{i}.Position, [maxX, maxY]));
                rois{i}.Waypoints = false(size(rois{i}.Waypoints));
                rois{i}.Label = num2str(i);
                i=i+1;
            end
                
        end
    end
    set(gcf, 'Position', ddd); % make the figure normal size
    zoom out
    % New UI by Hiro (vers 2C)
    thisisok=[];
    choice=questdlg('Happy with the ROIs?','Tell me','Yes','Try again :(','Yes');
    switch choice
    case 'Yes'
        thisisok=1;
    case 'Try again :('
        thisisok=0;
    case ''
        thisisok=[];
        disp('User canceled');
        return
    end
end
    

%% Save ROIs
saverois=[];
choice=questdlg('Save the ROIs?','Save','Yes');
switch choice
	case 'Yes'
        saverois=1;
	case 'No'
        saverois=0;
        disp('ROI not saved');
    case 'Cancel'
        saverois=[];
        disp('User canceled');
        return
end

if saverois
    [~,roisPath]=uiputfile(strcat(pathName,filesep,baseFileName,'.txt'),'Save for Igor Pro');
    save(strcat(roisPath,baseFileName,'.mat'),'rois');
    disp(['ROI saved as ',strcat(baseFileName,'.mat')]);
end

%% calculate mean intensity
disp(' '); disp('Please wait');
tic
%create roi label image
labelimg=zeros(size_first(:,1),size_first(:,2));
for i=1:numrois % make logical mask for each roi
    mask2{i}=createMask(rois{i});
    pixels(i)=sum(sum(mask2{i})); %pixels in roi
    labelimg=labelimg+i*mask2{i}; %tag pixels in label image
end
labelog=logical(labelimg);
[aa bb cc]=find(labelimg); % get pixels contained in rois
labelpx=length(cc); % find number of total pixels contained

intensity=zeros(numrois,num_images); %preallocate



[dirlist,~] = sort_nat({dirlist.name}'); %dirlist changes from struct to cell
    
    h50 = waitbar(0,'Calculating...');
    for i=1:num_images
        if ftype==1
            img = imread(cat(2,pname,fname), i, 'Info', info);
        else
            img=imread(strcat(mydir,filesep,dirlist{i}));
        end
        
        maskedmat=img(labelog);
        for j=1:numrois
            intensity(j,i)=mean(maskedmat(cc==j));
        end
        waitbar(i/num_images,h50)
    end
    close(h50);

toc

%% Update ROIs figure

h3 = waitbar(0,'Updating ROI appearance...');
for i=1:numrois
    rois{i}.Color = 'y';
    rois{i}.LineWidth = 0.5;
    waitbar(i/numrois,h3)
end
close(h3)

time=(0:1:num_images-1)/hz;

pmeth = 1;

%% save roi figure

if savePlots && saverois
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,strcat(roisPath,'rois.fig'));
    saveas(gcf,strcat(roisPath,'rois.png'));
    disp(['ROI saved as ',strcat(pathName,filesep,'rois')]);
end

%% plot traces
offset=[];
    rawtraces=figure;
    % find a buffer between traces
    buffer=0.1*(max(intensity(1,:))-min(intensity(1,:)));
    offset(1)=0;
    roitext=[];
    roitext(1)=text(0,max(intensity(1,1:9)),'1 ','HorizontalAlignment','right');
    for i=2:numrois
        % find offset so traces don't overlap, plus the small buffer
        offset(i)=min(intensity(i-1,:)-max(intensity(i,:))-buffer+offset(i-1));
        %hold on, plot(time,intensity(i,:)+offset(i))
        roitext(i)=text(0,max(intensity(i,1:9)+offset(i)),cat(2,num2str(i),' '),'HorizontalAlignment','right');
    end
    % matrix of intensity values with offsets applied
    offsetraw=transpose(intensity)+transpose(repmat(transpose(offset(:,:)),1,num_images));
    hold on, plot(time, offsetraw);
    offsetraw=transpose(offsetraw);
    set(gca,'yticklabel',[]);
    axis tight
    xlabel('Time (s)')
    barf=axis;
    set(gcf,'toolbar','figure');
    
    % new plot code
    zoom reset

    % add scroll and zoom buttons on each axis
    bp = uipanel('Position',[0 0 1 .05]);
    rp = uipanel('Position',[.95 0.05 .05 .95]);
    npos2=[.01 0.01 .09 .98];
    Sreset='zoom out,xsteP=get(gca,''xlim'');for i=1:numrois,oldP=get(roitext(i),''position'');set(roitext(i),''position'',[xsteP(1) oldP(2)]),end,xsteP=get(gca,''xlim'');ysteP=get(gca,''ylim'');';
    h2=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','View all','position',npos2,'callback',Sreset);
    set(gcf,'units','pixels')
    
    xsteP=get(gca,'xlim');
    Sleft=['set(gca,''xlim'',xsteP-.1*diff(xsteP));xsteP=get(gca,''xlim'');',...
        'for i=1:numrois,oldP=get(roitext(i),''position'');set(roitext(i),''position'',[xsteP(1) oldP(2)]),end'];
    Sxin=['set(gca,''xlim'',[xsteP(1)+.1*diff(xsteP) xsteP(2)-.1*diff(xsteP)]);',...
        'xsteP=get(gca,''xlim'');','for i=1:numrois,oldP=get(roitext(i),''position'');set(roitext(i),''position'',[xsteP(1) oldP(2)]),end'];
    Sxout=['set(gca,''xlim'',[xsteP(1)-.125*diff(xsteP) xsteP(2)+.125*diff(xsteP)]);',...
        'xsteP=get(gca,''xlim'');','for i=1:numrois,oldP=get(roitext(i),''position'');set(roitext(i),''position'',[xsteP(1) oldP(2)]),end'];
    Sright=['set(gca,''xlim'',xsteP+.1*diff(xsteP));xsteP=get(gca,''xlim'');',...
        'for i=1:numrois,oldP=get(roitext(i),''position'');set(roitext(i),''position'',[xsteP(1) oldP(2)]),end'];
    
    ysteP=get(gca,'ylim');
    Sup=['set(gca,''ylim'',ysteP+.1*diff(ysteP));ysteP=get(gca,''ylim'');'];
    Syin='set(gca,''ylim'',[ysteP(1)+.1*diff(ysteP) ysteP(2)-.1*diff(ysteP)]);ysteP=get(gca,''ylim'');';
    Syout='set(gca,''ylim'',[ysteP(1)-.125*diff(ysteP) ysteP(2)+.125*diff(ysteP)]);ysteP=get(gca,''ylim'');';
    Sdown='set(gca,''ylim'',ysteP-.1*diff(ysteP));ysteP=get(gca,''ylim'');';
    
    posL=[.35 0.01 .05 .98];
    posXP=[.45 0.01 .05 .98];
    posXM=[.5 0.01 .05 .98];
    posR=[.6 0.01 .05 .98];
    hL=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','<--','position',posL,'callback',Sleft);
    hXP=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','+','position',posXP,'callback',Sxin);
    hXM=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','-','position',posXM,'callback',Sxout);
    hR=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','-->','position',posR,'callback',Sright);
    
    posU=[.01 0.6 .98 .05];
    posYP=[.01 0.5 .98 .05];
    posYM=[.01 0.45 .98 .05];
    posD=[.01 0.35 .98 .05];
    hU=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','/\','position',posU,'callback',Sup);
    hYP=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','+','position',posYP,'callback',Syin);
    hYM=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','-','position',posYM,'callback',Syout);
    hD=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','V','position',posD,'callback',Sdown);
    set(gcf,'units','pixels')
    




% New UI by Hiro (vers 2C)
dodff=[];
choice=questdlg('Convert to dF/F?','dF/F','Yes');
switch choice
    case 'Yes'
        dodff=1;
    case 'No'
        dodff=0;
        disp('Not converted to dF/F');
    case 'Cancel'
        dodff=[];
        disp('User canceled');
        return
end

%% convert to dF/F
if dodff
    
    % New UI by Hiro (vers 2C)
    uiwait(msgbox({'Your job is to select the baseline:' ...
        'Click and drag to make a rectangle over the flat part of the traces.' ...
        'A rectangle may contain many traces.' ...
        'Make as many rectangles as needed until all the traces are selected.' ...
        'Press ENTER when you are finished.' ...
        'Click OK to begin...'},'Instructions'));
    
    figure(rawtraces);
    baselines=zeros(numrois,2); % start and end points for baselines
    baselineframes=[]; % frames to calculate avg image for F0
    numrects=0;
    while sum(baselines(:,2)==0) % not all endpoints are set
        disp(cat(2,num2str(sum(baselines(:,2)==0)),' trace(s) remain(s)'));
        numrects=numrects+1;
        rectangs=imrect;
        trash=input('Adjust the rectangle if you want, then hit Enter when ready');
        rectcoords=rectangs.getPosition();
        xmin=rectcoords(1);
        xmax=rectcoords(1)+rectcoords(3);
        [a1 a2]=min(abs(time-xmin));
        xmin=a2;
        [a1 a2]=min(abs(time-xmax));
        xmax=a2;
        ymin=rectcoords(2);
        ymax=rectcoords(2)+rectcoords(4);
        
        undertop=offsetraw(:,xmin:xmax)<=ymax; % contained in rect bounds
        overbottom=offsetraw(:,xmin:xmax)>=ymin; % ditto
        rowscontained=sum(undertop.*overbottom,2)>0;
        baselines(rowscontained,1)=xmin; % set min and max values
        baselines(rowscontained,2)=xmax;  
        baselineframes=cat(2,baselineframes,xmin:xmax); % add on these frames
        
        %graphical indication of which traces have baselines?
        %delete rectangle
        delete(rectangs);
        %draw all lines in black
        hold on, plot(time, transpose(offsetraw(rowscontained,:)),'k')
        %draw baseline regions in green
        hold on, plot(time(xmin:xmax), transpose(offsetraw(rowscontained,xmin:xmax)),'g')
    end
% close(rawtraces);

    % choose dF/F method
    disp(' ');
    disp('Choose dF/F calculation method: ');
    disp('1. Standard (divide by the baseline of each ROI)');
    %disp('2. Divide by avg intensity of entire image in baseline frames');
    disp('2. There is no option 2')
    disp('3. Divide by averaging intensity of all ROIs in baseline frames')
    disp('4. Option 3 but with 30% weighting for ROI baselines');

    
    % New UI by Hiro (vers 2C)
    disp(' ');
    dffmeth=[];
    choice=questdlg('Which dF/F method?','dF/F Method','Standard','Option 3','Option 4','Standard');
    switch choice
        case 'Standard'
            dffmeth=1;
            disp('You chose the Standard dF/F method');
        case 'Option 3'
            dffmeth=3;
            disp('You chose Option 3 for the dF/F method');
        case 'Option 4'
            dffmeth=4;
            disp('You chose Option 4 for the dF/F method');
    end
    
    % remove figure
    close(rawtraces);
    
    for i = 1:numrois % avg value of each roi's baseline
       blavg(i)=mean(intensity(i,baselines(i,1):baselines(i,2)));
    end
    
    
    if dffmeth==2
        % find all baseline frames and get avg intensity in all those pixels
        blf=unique(baselineframes);
        denominator=sum(sum(sum(stack(:,:,blf))))/(size_first(:,1)*size_first(:,2)*length(blf));
    end
    if dffmeth==3 || dffmeth==4
        blpts=[];
        for i=1:numrois
            blpts=cat(2,blpts,intensity(i,baselines(i,1):baselines(i,2)));
        end
        denominator=mean(blpts);
        if dffmeth==4
            newdenom=blavg*.3+.7*denominator; 
        end
    end
    
     

    %%
    delf=[];
    for i=1:numrois
        if dffmeth==1
            delf(i,:)=100*(intensity(i,:)/blavg(i)-1);
        end
        if dffmeth==2 || dffmeth==3
            delf(i,:)=100*(intensity(i,:) - blavg(i))/denominator;
        end
        if dffmeth==4
            delf(i,:)=100*(intensity(i,:) - blavg(i))/newdenom(i);
        end
       
    end
    

    figure
    offset=[];
    buffer=0.1*(max(delf(1,:))-min(delf(1,:)));
    offset(1)=0;
    roitext=[];
    roitext(1)=text(0,max(delf(1,1:10)),'1 ','HorizontalAlignment','right');
    for i=2:numrois
        offset(i)=0-max(delf(i,:))+min(delf(i-1,:)-buffer+offset(i-1));
        roitext(i)=text(0,max(delf(i,1:10)+offset(i)),cat(2,num2str(i),' '),'HorizontalAlignment','right');
    end
    offsetdff=transpose(delf)+transpose(repmat(transpose(offset(:,:)),1,num_images));
    hold on, plot(time, offsetdff);
    set(gca,'yticklabel',[]);
    axis tight
    
    % put 10% df/F bar on the trace with the largest amplitude change
    [qq ww]=find(diff(offset)==min(diff(offset)),1);
    maxspace=abs(min(diff(offset)));
    if maxspace>90
        barsize=40;
    else
        if maxspace>50
            barsize=20;
        else
            barsize=10;
        end
    end
    
    if numrois==1
        bottom=1;
    else
        bottom=offset(ww+1);
    end
    hold on,line([time(2) time(2)],[bottom bottom+barsize],'Color','k','LineWidth',5);
    text(time(3),bottom+(barsize/2),cat(2,num2str(barsize),'% dF/F'));
    xlabel('Time (s)')  
    
    barf=axis; 
    set(gcf,'toolbar','figure');
    % new plot code, condensed
    zoom reset
 
    % add scroll and zoom buttons on each axis
    bp = uipanel('Position',[0 0 1 .05]);
    rp = uipanel('Position',[.95 0.05 .05 .95]);

    h2=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','View all','position',npos2,'callback',Sreset);
    set(gcf,'units','pixels')
    
    xsteP=get(gca,'xlim');
    ysteP=get(gca,'ylim');
    
    hL=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','<--','position',posL,'callback',Sleft);
    hXP=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','+','position',posXP,'callback',Sxin);
    hXM=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','-','position',posXM,'callback',Sxout);
    hR=uicontrol('Parent',bp,'units','normalized','style','pushbutton','string','-->','position',posR,'callback',Sright);
    
    hU=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','/\','position',posU,'callback',Sup);
    hYP=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','+','position',posYP,'callback',Syin);
    hYM=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','-','position',posYM,'callback',Syout);
    hD=uicontrol('Parent',rp,'units','normalized','style','pushbutton','string','V','position',posD,'callback',Sdown);
    set(gcf,'units','pixels')
      
end
%%


% New UI by Hiro (vers 2C)
saving=[];
choice=questdlg('Save for Igor Pro?','Igor Pro','Original only','dF/F only','Both','Both');
switch choice
    case 'Original only'
        saving=1;
    case 'dF/F only'
        saving=2;
    case 'Both'
        saving=3;
    case ''
        saving=4;
end

if saving ~=4
     disp(' ');

    
    % New UI by Hiro (vers 2C)
    [fname,tracesPath,filterIndex]=uiputfile(strcat(pathName,filesep,baseFileName,'.txt'),'Save for Igor Pro');
    if filterIndex==0
        disp('User canceled');
        return
    end
    
    fname=baseFileName;
    
    curdir=cd;
    cd(tracesPath)
    if saving==1 || saving==3
        dlmwrite(strcat(fname,'_original.txt'), intensity', 'delimiter','\t','newline','pc');
        disp([strcat(fname,'_original.txt'),' saved in ',tracesPath]);
        
        % Mat2Igor binary by Hiro (vers 3)
        HiroMat2Igor(fname,hz,0,intensity,imgWidthMicron,imgWidth,imgHeight,rois);
    end
    if saving==2 || saving==3
        dlmwrite(strcat(fname,'_dFoverF.txt'), delf', 'delimiter','\t','newline','pc');
        disp([strcat(fname,'_dFoverF.txt'),' saved in ',tracesPath]);
        
        HiroMat2Igor(fname,hz,1,delf,imgWidthMicron,imgWidth,imgHeight,rois);
    end
    cd(curdir)
end

%% save trace figure
if savePlots
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,strcat(tracesPath,'trace.fig'));
    saveas(gcf,strcat(tracesPath,'trace.png'));
    disp(['trace plot saved as ',strcat(pathName,filesep,'trace')]);
end
%% clean up
clear asdf aa bb cc maskedmat stack totalmask %clean up larger matrices
clear a1 a2 baseFileName baselineframes baselines blf bottom buffer choice ddd dirlist % and other junk
clear barsize blavg blpts bp denominator dffmeth h50 h60 img maxspace ministack
clear newdenom quitcount roitext def dlg_title ext savePlots
clear dodff first first2 newimg ftype fname2 h i img1 imgFrom imgTo info2 isDone j fname1 blah1 blah2 fcn
clear labelimg labelog labelpx loadrois mask2 maskedvec mydir num2avg
clear num_images2 numrects offset overbottom path1 pixels pmeth pname2
clear qq rawtraces rectangs rectcoords rois roipts saverois saving sumimg rowscontained
clear size_first sys texts thisisok trash undertop user ww xmax xmin ymax ymin

disp(' ');
disp('Done: Remember to save your ROI image (size window appropriately first)');