function [ExpImage_Image,data1,td] = TescanFrame_Load(ExpImage_Filename,Frame_Number)
%ECP_LOAD Summary of this function goes here
%   Detailed explanation goes here
data1_fields={'ImageStripSize','ViewFieldsCountX','ViewFieldsCountY','AcceleratorVoltage','Magnification','StageTilt','StageTiltY''','StageRotation','PixelSizeX','PixelSizeY'};
data1_fields{end+1}=strcat('Detector',int2str(Frame_Number-1));
ExpImage_Image=flipud(imread(ExpImage_Filename));
[data1,td]=header_read_tescan(ExpImage_Filename,data1_fields);

I2_stripsize=str2double(data1{1});
% take off the data bar
ExpImage_Image=flipud(ExpImage_Image(I2_stripsize+1:end,:));

%now extract the single frame
ECP_NumX=str2double(data1{2});
ECP_NumY=str2double(data1{3});

ECP_Wid=size(ExpImage_Image,2)/ECP_NumX;
ECP_Hig=size(ExpImage_Image,1)/ECP_NumY;

%read one of the frames
ECP_Frames=(reshape(1:ECP_NumX*ECP_NumY,ECP_NumY,ECP_NumX))';
[ECP_GridX,ECP_GridY]=meshgrid(1:ECP_NumX,1:ECP_NumY);
ECP_FrameNum=find(ECP_Frames==Frame_Number);
ECP_XExt=[1+(ECP_GridX(ECP_FrameNum)-1)*ECP_Wid:(ECP_GridX(ECP_FrameNum))*ECP_Wid];
ECP_YExt=[1+(ECP_GridY(ECP_FrameNum)-1)*ECP_Hig:(ECP_GridY(ECP_FrameNum))*ECP_Hig];

ExpImage_Image=ExpImage_Image(ECP_YExt,ECP_XExt);
end

