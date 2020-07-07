% ROI time measurement and dF/F tool
% Zeke Barger
% updated 25.03.14 (last change = drift correction)

disp(' '); disp(' ');
disp('   For each prompt, enter the number of your');
disp('      selection and press the enter key.');
% [user,sys] = memory;
% disp(cat(2,'         Available memory: ',num2str(sys.PhysicalMemory.Available/1000000000),' GB'));
% disp(' ');

hz=[];
while isempty(hz)
    hz=input('Enter framerate (Hz): ');
end

disp(' ');
loadrois=[];
while isempty(loadrois)
    loadrois=input('(1) Draw new ROIs or (2) load saved ROIs : ')-1;
    if ~isempty(loadrois)
        if loadrois ~= 0 && loadrois ~= 1
            loadrois=[];
        end
    end
end


disp(' '); disp('Are your data in a...');
ftype=[];
while isempty(ftype)
    ftype=input('(1) TIFF stack [one file] or (2) image sequence [many TIFF files]? ');
    if ~isempty(ftype)
        if ftype ~= 1 && ftype ~= 2
            ftype=[];
        end
    end
end

% do other prompts first
disp(' ');
img1=[];
while isempty(img1)
    img1=input('Place ROIs on 1) brightfield/other image, 2) average image : ');
    if ~isempty(img1)
        if img1 ~= 1 && img1 ~= 2
            img1=[];
        end
    end
end


if loadrois % get saved rois
    disp('Select your saved ROI .mat file');
    [fname1 path1]=uigetfile('*.*','Your saved ROIs');
    load(cat(2,path1, fname1));
    numrois=size(rois,2)-1;
end


%% load files, except not
% enoughMem=1;

disp(' ');
if ftype==1
    disp('Select your TIFF stack');
    [fname,pname]=uigetfile('*.*','Your TIFF stack!');
    info = imfinfo(cat(2,pname,fname));
    num_images = numel(info);
    first = imread(cat(2,pname,fname), 1, 'Info', info);
    size_first=size(first);  %get size of first image to preallocate space for imported images
    
%     try
%         
%     stack=zeros(size_first(:,1),size_first(:,2),num_images);
%     h = waitbar(0,'Loading images...');
%     for k = 1:num_images
%         stack(:,:,k) = imread(cat(2,pname,fname), k, 'Info', info);
%         waitbar(k/num_images,h)
%     end
%     close(h)
%     
%     catch
%         
%         disp('Not enough memory to load the sequence.');
%         disp('Performance may suffer.');
%         enoughMem=0;
%         
%     end
else
    % Open image sequence using gui 
    disp('Select the directory of your image sequence (must not be RGB)')
    mydir=uigetdir('C:\','Your image sequence!');            %use graphical user interface to set file path in mydir
    dirlist=dir(strcat(mydir,filesep,'*.tif'));                    %save file path in structure array
    first=imread(strcat(mydir,filesep,dirlist(1).name));           %pull out first image from stack
    size_first=size(first);                                        %get size of first image to preallocate space for imported images
    num_images=length(dirlist);
    
%     try
%     
%     stack=zeros(size_first(:,1),size_first(:,2),num_images);  %preallocate space for images
% 
%     % create stack for images from dirlist
%     h = waitbar(0,'Loading images...');
%     for i=(1:length(dirlist(:)))
%         stack(:,:,i)=imread(strcat(mydir,filesep,dirlist(i).name));
%         waitbar(i/num_images,h)
%     end
%     close(h)
%     
%     catch
%         
%         disp('Not enough memory to load the sequence.');
%         disp('Performance may suffer.');
%         enoughMem=0;
%         
%     end
end



if img1==2
    % find average image using first quarter of stack
    num2avg=round(num_images/4);
    if num2avg > 200 % but not more than 200 images
       num2avg=200; 
    end
    
%     if enoughMem %analyse stack
%         first2=mean(stack(:,:,1:num2avg),3);
%     else % load in the first few images
        h60 = waitbar(0,'Creating average image');
        ministack=zeros(size_first(:,1),size_first(:,2),num2avg);
        
        if ftype==1
            for i=1:num2avg
                ministack(:,:,i) = imread(cat(2,pname,fname), i, 'Info', info);
                 waitbar(i/num2avg,h60)
            end
            first2=mean(ministack(:,:,1:num2avg),3);
        else
            for i=1:num2avg
                ministack(:,:,i)=imread(strcat(mydir,filesep,dirlist(i).name));
                 waitbar(i/num2avg,h60)
            end
            first2=mean(ministack(:,:,1:num2avg),3);
        end
        close(h60)
%     end
else
    % ask for new image, making sure it has the same dimensions
    first2=[0 0];
    disp(' '); disp('Select the brightfield/other image. The dimensions must match the film');
    while [size(first,1) size(first,2)] ~= [size(first2,1) size(first2,2)]  
        [fname2,pname2]=uigetfile('*.*','Brightfield/average image');
        info2 = imfinfo(cat(2,pname2,fname2));
        num_images2 = numel(info2);
        first2 = imread(cat(2,pname2,fname2), 1, 'Info', info2);
        if [size(first,1) size(first,2)] ~= [size(first2,1) size(first2,2)]
           disp('Your dimensions did not match, try again'); 
        end
    end
end

intensity=[];
pixels=[];
mask2={};
maskedvec=[];

if loadrois % get saved rois
%     [fname1 path1]=uigetfile('*.*','Your saved ROIs');
%     load(cat(2,path1, fname1));
%     numrois=size(rois,2)-1;
    cellfig=figure;
    colormap(gray);
    imagesc(first2)
    axis image
    axis off
    
    firstP=rois{1}.getPosition();
    firstROI=impoly(gca,firstP,'closed',true);
    setVerticesDraggable(firstROI,false);
    hb1=msgbox('If necessary, move ROI 1 to correct for drift. Then press enter in the command window');
    uiwait(hb1);
    trash=input('If necessary, move ROI 1 to correct for drift. Then press enter here');
    firstP2=firstROI.getPosition();
    drift=[firstP2(1,1)-firstP(1,1) firstP2(1,2)-firstP(1,2)]; % change in position
    delete(firstROI);
    
    hh2 = waitbar(0,'Loading ROIs...');

    for i=1:numrois
        fulllist=rois{i}.getPosition();
        fulllist=fulllist+repmat(drift,size(fulllist,1),1); % correct for drift
        % aim for 40 points in each roi
        if size(fulllist,1)<40
            step1=1;
        else
            if size(fulllist,1)>300 %unless there are hella points
                step1=size(fulllist,1)/100;
            else
                step1=size(fulllist,1)/40;
            end
        end    
        rois{i}=impoly(gca,fulllist(round(1:step1:size(fulllist,1)),:),'closed',true);
        setVerticesDraggable(rois{i},false);
        waitbar(i/numrois,hh2)
       
    end
    close(hh2) 
   

    pause(.1)
%% place new rois
else
    thisisok=0;
    while thisisok == 0 % get rois until user is satisfied with them
        % measure intensity in rois       
        rois={};
        texts=[];
        numrois=[];
        cellfig=figure;
        colormap(gray);
        imagesc(first2)
        axis image
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'YTick',[]);
        xlabel('Press Esc to undo. To finish, shift+click inside the image, then click outside it');
        zoom reset
        
        ddd=get(gcf,'Position'); % make the figure window play nice
        blah1=get(gca,'XLim');
        blah2=get(gca,'YLim');
        fcn = makeConstrainToRectFcn('imfreehand',blah1,blah2);

        i=1;
        quitcount=0;
        set(gcf,'toolbar','figure');
        % while i < (numrois+1)
        hb1=msgbox('Draw ROIs. To delete the most recent ROI, press Esc. When finished, shift+click inside the image and then click outside the image');
        uiwait(hb1);

        while quitcount==0
            rois{i}=imfreehand('PositionConstraintFcn',fcn);
            if sum(size(rois{i})) > 0 % user didn't press escape
                asdf=rois{i}.getPosition;
                if isempty(asdf)
                    quitcount=1;
                else
                    texts(i)=text(asdf(1,1),asdf(1,2),num2str(i),'HorizontalAlignment',...
                        'center','BackgroundColor',[.5 .5 .5],'Margin',.01,'FontName','Arial','FontSize',10);
                    i=i+1;
                end
            else % user hit escape, delete the previous roi
                if i>1
                    i=i-1;
                    delete(texts(i));
                    delete(rois{i});
                end
                    
            end
        end
        
        set(gcf, 'Position', ddd); % make the figure normal size
        zoom out

        disp(' ');
        thisisok=[];
        while isempty(thisisok) % possibly redo the whole thing
            thisisok=abs(input('Do these ROIs look ok? 1=yes, 2=restart : ')-2);
            if ~isempty(thisisok)
                if thisisok ~= 0 && thisisok ~= 1
                    thisisok=[];
                end
            end
        end         
    end
    
    numrois=i-1;
    %possibly save rois
    disp(' ');
    saverois=[];
    while isempty(saverois)
        saverois=abs(input('Do you want to save these ROIs for later? 1=yes, 2=no : ')-2);
        if ~isempty(saverois)
            if saverois ~= 0 && saverois ~= 1
                saverois=[];
            end
        end
    end
    
    if saverois
        uisave('rois');  
    end
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
%totalmask=repmat(labelog,[1,1,num_images]); %make the mask 3d
[aa bb cc]=find(labelimg); % get pixels contained in rois
labelpx=length(cc); % find number of total pixels contained

intensity=zeros(numrois,num_images); %preallocate
%maskedvec=stack(totalmask); %apply mask to df stack (produces huge vector)
%maskedmat=vec2mat(stack(totalmask),labelpx);

% if enoughMem
% % gigantic step, gets all intensity values for all roi pixels in all frames
% maskedmat=transpose(reshape(stack(totalmask),labelpx,num_images));
% for i=1:num_images
%     for j=1:numrois
%         intensity(j,i)=mean(maskedmat(i,cc==j));
%     end
% end
% 
% else %need to load each file and process it
    h50 = waitbar(0,'Calculating...');
    for i=1:num_images
        if ftype==1
            img = imread(cat(2,pname,fname), i, 'Info', info);
        else
            img=imread(strcat(mydir,filesep,dirlist(i).name));
        end
        
        maskedmat=img(labelog);
        for j=1:numrois
            intensity(j,i)=mean(maskedmat(cc==j));
        end
        waitbar(i/num_images,h50)
    end
    close(h50);
    
% end


toc

if loadrois==0 % gotta make this work for both options???
    % update ROI figure appearance
    for i=1:numrois
        roipts=rois{i}.getPosition;
        roipts(end+1,:)=roipts(1,:);
        hold on,plot(roipts(:,1),roipts(:,2),'y','LineWidth',0.5);
        delete(rois{i});
        uistack(texts(i),'top'); %move text in front of roi
    end
else
    %just delete rois or something NO
    h3 = waitbar(0,'Updating ROI appearance...');
    for i=1:numrois
        roipts=rois{i}.getPosition;
        roipts(end+1,:)=roipts(1,:);
        hold on,plot(roipts(:,1),roipts(:,2),'y','LineWidth',0.5);
        asdf=rois{i}.getPosition;
        texts(i)=text(asdf(1,1),asdf(1,2),num2str(i),'HorizontalAlignment',...
           'center','BackgroundColor',[.5 .5 .5],'Margin',.01,'FontName','Arial','FontSize',10);   
        delete(rois{i});
        waitbar(i/numrois,h3)
    end
    close(h3)
end

time=(0:1:num_images-1)/hz;
% disp(' '); % plotting in multiple figures was not so good
% pmeth=input('Plot traces in one figure (1) or multiple windows (2)? ');
% if isempty(pmeth)
    pmeth = 1;
% end
%% plot traces
offset=[];
% if pmeth==1
    rawtraces=figure;
    % find a buffer between traces
    buffer=0.1*(max(intensity(1,:))-min(intensity(1,:)));
    offset(1)=0;
    roitext=[];
    %plot(time,intensity(1,:))
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
    
    %% new plot code
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
    

disp(' ');
dodff=[];
while isempty(dodff)
    dodff=abs(input('convert to %dF/F? 1=yes 2=no : ')-2);
    if ~isempty(dodff)
        if dodff ~= 0 && dodff ~= 1
            dodff=[];
        end
    end
end   

%% convert to dF/F
if dodff
%     disp(' ');
%     disp('Draw a rectangle and position it around a part of the trace where');
%     disp('there is no activity. The top and bottom bounds do not matter.');
    disp(' ');
    disp('INSTRUCTIONS: Choose which parts of the traces to use as baselines.');
    disp('Begin by drawing a rectangle where the traces are flat. Include as');
    disp('many traces as you want. You will be able to draw more rectangles'); 
    disp('until all traces are selected. Press ENTER when you are satisfied');
    disp('with each rectangle.');
    disp('Press enter to continue...');
    input('');
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
 
    %% choose dF/F method
    disp(' ');
    disp('Choose dF/F calculation method: ');
    disp('1. Standard (divide by each ROI`s baseline)');
    %disp('2. Divide by avg intensity of entire image in baseline frames');
    disp('There is no option 2')
    disp('3. Divide by avg intensity of all ROIs in baseline frames')
    disp('4. Option 3 but with 30% weighting for ROI baselines');
    dffmeth=[];
    while isempty(dffmeth)
        dffmeth=input('Your selection: ');
        if ~isempty(dffmeth)
            if sum(dffmeth==[1 3 4])==0
                dffmeth=[];
            end
        end
    end 
    
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
    
     
    % old Zeke method
    % denominator=sum(sum(sum(stack(:,:,xmin:xmax))))/(size_first(:,1)*size_first(:,2)*(xmax-xmin+1));
    
    % newer dF/F0 calculation, not perfect: 
    % dF = signal - avg value of baseline section
    % F0 = average intensity of entire image in all baseline frames
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
       % very old Zeke method
       % delf(i,:)=100*(intensity(i,:) - mean(intensity(i,xmin:xmax)))/denominator; 
       
    end
    
%     if pmeth==1
%        close(rawtraces); 
%     end
    figure
    offset=[];
    buffer=0.1*(max(delf(1,:))-min(delf(1,:)));
    offset(1)=0;
    roitext=[];
    roitext(1)=text(0,max(delf(1,1:10)),'1 ','HorizontalAlignment','right');
    for i=2:numrois
        offset(i)=0-max(delf(i,:))+min(delf(i-1,:)-buffer+offset(i-1));
        %hold on, plot(time,intensity(i,:)+offset(i))
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
    %% new plot code, condensed
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
disp(' ');
saving=[];
while isempty(saving)
    saving=input('(For Igor) 1) save original, 2) save dF/F, 3) save both, 4) exit : ');
    if ~isempty(saving)
        if saving ~= 1 && saving ~= 2 && saving ~= 3 && saving ~= 4
            saving=[];
        end
    end
end   

if saving ~=4
    disp(' ');
    mydir=uigetdir('','Where to save?');
    fname=[]; 
    while isempty(fname)
        fname=input('Enter base filename: ','s');
    end
    curdir=cd;
    cd(mydir)
    if saving==1 || saving==3
        dlmwrite(cat(2,fname,'_original','.txt'), intensity', 'delimiter','\t','newline','pc');
        disp(cat(2,'Saving ',fname,'_original','.txt'));
    end
    if saving==2 || saving==3
        dlmwrite(cat(2,fname,'_dFoverF','.txt'), delf', 'delimiter','\t','newline','pc');
        disp(cat(2,'Saving ',fname,'_dFoverF','.txt'));
    end
    cd(curdir)
end

clear asdf aa bb cc maskedmat stack totalmask %clean up larger matrices
clear a1 a2 baselineframes baselines blf bottom buffer ddd dirlist % and other junk
clear barsize blavg blpts bp denominator dffmeth h50 h60 img maxspace ministack
clear newdenom quitcount roitext 
clear dodff first first2 ftype fname2 h i img1 info2 j fname1 blah1 blah2 fcn
clear labelimg labelog labelpx loadrois mask2 maskedvec mydir num2avg
clear num_images2 numrects offset overbottom path1 pixels pmeth pname2
clear qq rawtraces rectangs rectcoords roipts saverois saving rowscontained
clear size_first sys texts thisisok trash undertop user ww xmax xmin ymax ymin
disp('Remember to save your ROI image (size window appropriately first)');