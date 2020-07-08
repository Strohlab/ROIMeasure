%% Igor2Mat
% Extract saved ROIs from Igor Bin file
% Saleh Altahini 


rois={};        % empty table to load ROIs into

% open the file
disp('Select your Igor .bin file');
isDone=0;
[fname, fpath]=uigetfile('*.bin','Select your Igor .bin file');
if ~fname
	disp('User canceled');
	return
end

fileID=fopen(strcat(fpath, fname),'r');
binData = fread(fileID, Inf, 'double');

% check the header
fheader = binData(1:binData(1));
if fheader(3) ~= 1.01 || fheader(2) ~= 2
    error('Unsupported Igor .bin file!!');
end

% get the ROI data
dataPackage = binData((fheader(1)+1):end);
traceLen = dataPackage(1)*dataPackage(2);
traceLen = traceLen + 3;
roiData = dataPackage(traceLen:end);

% ask for figure image
disp('Select the average image.');
[fname2, fpath2]=uigetfile(strcat(fpath,'/*.*'),'Brightfield/average image');
imageInfo = imfinfo(cat(2,fpath2,fname2));
bgImage = imread(cat(2,fpath2,fname2), 1, 'Info', imageInfo);

% prepare the figure
cellfig=figure;
colormap(gray);
imagesc(bgImage);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);

% convert ROI array to actual ROIs
n = 1;
i = 1;
while n<length(roiData)
    % convert ROI location array to n-by-2 matrix
    roiLen = roiData(n)*roiData(n+1);
    curROI = roiData((n+2):(roiLen+n+1));
    curROI = reshape(curROI,[roiData(n),roiData(n+1)]);
    
    % add the roi to figure
    rois{i} = images.roi.Freehand(gca,'Position',curROI, 'Color', 'r');
    rois{i}.Waypoints = false(size(rois{i}.Waypoints));
    rois{i}.Label = num2str(i);
    
    % counters++
    n = n+roiLen+2;
    i = i+1;
end

choice=questdlg('Save the ROIs?','Save','Yes');
switch choice
	case 'Yes'
        uisave('rois',strcat(fpath,'recoveredROIS','.mat'));
        disp(['ROI saved at ',strcat(fpath,'recoveredROIS','.mat')]);
	case 'No'
        saverois=0;
        disp('ROI not saved');
    case 'Cancel'
        saverois=[];
        disp('User canceled');
        return
end

clear bgImage binData cellfig curROI dataPackage fheader fileID fname fname2
clear fpath fpath2 i n imageInfo roiData roiLen rois traceLen