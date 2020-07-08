function HiroMat2Igor(fileName,hz,trace_type,trace,imgWidthMicron,imgWidth,imgHeight,rois)
%% Export data for Igor Pro as a binary file
% Copyright 2016 Hiro Watari
% file header length
% file version 1.01
% sampling frequency (Hz)
% type: 0=orig, 1=dF/F
% image height (pixels)
% image width (pixels)

% data header length
% number of rows
% number of cols

tic

m2i_nFileHeader=8;
m2i_nDataHeader=2;
m2i_FileVersion=1.01;
% m2i_fs=hz;          % parameter
% m2i_type=trace_type;         % parameter
% m2i_imgWidth=imgWidth;   % parameter
% m2i_imgHeight=imgHeight;  % parameter
% m2i_imgWidthMicron=imgWidthMicron;  % parameter

% File Header (vector)
m2i_FileHeader=zeros(m2i_nFileHeader,1);
m2i_FileHeader(1)=m2i_nFileHeader;
m2i_FileHeader(2)=m2i_nDataHeader;
m2i_FileHeader(3)=m2i_FileVersion;
% paramters
m2i_FileHeader(4)=hz;
m2i_FileHeader(5)=trace_type;
m2i_FileHeader(6)=imgWidth;
m2i_FileHeader(7)=imgHeight;
m2i_FileHeader(8)=imgWidthMicron;

% Data as matrix (time x ROI)
if trace_type
    % dF/F
    m2i_Data=trace';
else
    % original
    m2i_Data=trace';
end

% Data Header (vector)
m2i_DataHeader=zeros(m2i_nDataHeader,1);
m2i_nRows=size(m2i_Data,1);
m2i_nCols=size(m2i_Data,2);

m2i_DataHeader(1)=m2i_nRows;
m2i_DataHeader(2)=m2i_nCols;

% Vectorize data
m2i_Data=reshape(m2i_Data,[],1);

% Concatenate data header and data
m2i_DataPack=cat(1,m2i_DataHeader,m2i_Data);

% Loop through and copy ROI data
for i=1:m2i_nCols
    %m2i_Data=[];
    m2i_Data=rois{i}.Position;
    
    % Data Header
    m2i_DataHeader(1)=size(m2i_Data,1);
    m2i_DataHeader(2)=size(m2i_Data,2);
    
    % Vectorize
    m2i_Data=reshape(m2i_Data,[],1);
    
    % Concatenate data header and data
    m2i_DataROI=cat(1,m2i_DataHeader,m2i_Data);
    
    % Concatenate this to everything so far
    m2i_DataPack=cat(1,m2i_DataPack,m2i_DataROI);
end

% Concatenate file header with everything
m2i_Package=cat(1,m2i_FileHeader,m2i_DataPack);

if ~trace_type
    fileName=strcat(fileName,'_original.bin');
else
    fileName=strcat(fileName,'_dff.bin');
end

%Write to file
fileID=fopen(fileName,'w');
fwrite(fileID,m2i_Package,'double');
fclose(fileID);
toc

end