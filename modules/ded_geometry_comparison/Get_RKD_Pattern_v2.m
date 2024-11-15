%% Get_RKD_Pattern_v2.m
% T. B. Britton and T. Zhang 
% October 2024
% This script reads in RKD data *TFS EBSD software), process and output a
% selected RKD pattern

%%
clear;
home;
close all;

%% Input Variables

mtexLocation='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';
astroLocation='D:\AstroEBSD';


% path to data file
TFSDataLoc='D:\2024-11-04T15_22_24_Al_ben\2024-11-04T16_49_32';
InputUser.Phase_Input={'Aluminium'};
raw_loc=fullfile(TFSDataLoc,'\raw');
raw_folders=dir(raw_loc);
raw_folder=fullfile(raw_loc,raw_folders(3).name);
raw_pats_path=fullfile(raw_folder,'f_00000.pak');

% Open the file for reading in binary mode
fileID = fopen(raw_pats_path, 'rb');
rawPatsAll = fread(fileID,inf,'uint16');
fclose(fileID);

rawPatternSize = 549;
outputPatSize = 479;
%% Load MTEX

%start mtex if needed
try EBSD;
catch
    run(fullfile(mtexLocation,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(astroLocation,"start_AstroEBSD.m"));
end

RTM.Phase_Folder = fullfile(astroLocation,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% Load the EBSD data

[ebsd,ebsd_header] = bIDX_to_EBSD(TFSDataLoc);

% plotting convention - TFS
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% create a reference pattern from the Dynamics bin file
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},RTM.Phase_Folder);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Calculates the background from the stack of patterns in the raw file

% Note to Ben - I deliberately bring them out of a sunfunction to make it
% quicker, since we only need to do it once. Delete this comment upon
% publishing...

background=zeros(rawPatternSize, rawPatternSize);

s_size=256*256;
raw_data_vector = zeros(s_size*4, 1);
total_num = 1500; % limited array size, we will stop after counting this many patterns - change as needed

for n=1:total_num
    range_min=1+(n-1)*512*512;
    range_max=512*512+(n-1)*512*512;
    raw_data_vector=raw_data_vector+rawPatsAll(range_min:range_max);
end

raw_data_vector = raw_data_vector / total_num;
%seperate the chips
data_0=raw_data_vector(0*s_size+1:1*s_size); data_0=reshape(data_0,[256 256]);
data_1=raw_data_vector(1*s_size+1:2*s_size); data_1=reshape(data_1,[256 256]);
data_2=raw_data_vector(2*s_size+1:3*s_size); data_2=reshape(data_2,[256 256]);
data_3=raw_data_vector(3*s_size+1:4*s_size); data_3=reshape(data_3,[256 256]);

%arrange these in the array
background([1:256]+35,1:256)=rot90(data_2,3);
background(1:256,257:512)= rot90(data_1,2);
background([257:512]+35,[1:256]+35)= rot90(data_3,4);
background(257:512,[257:512]+35)=rot90(data_0, 2);

background = background(36:512,36:512);

% you can have a look at the background here
% figure; imagesc(BG); axis xy; axis image; colormap('gray');

%% Select a pattern to plot
x_val=16;
y_val=13;
%find the number of grid points in X and Y
x_max=max(ebsd.prop.x_in(:))+1;
y_max=max(ebsd.prop.y_in(:))+1;
pat_num=(x_val+1)+y_val*x_max;

% Read in a single pattern
[single_pattern] = rEBSD_TFS_Raw_Single(rawPatsAll, pat_num);

% correct by the background
single_pattern_bgcor = single_pattern ./ background;

% generate the "noise"
single_std=nanstd(single_pattern_bgcor(:));
single_mean=nanmean(single_pattern_bgcor(:));

noiseMask=randn(size(single_pattern_bgcor))*single_std+single_mean;

single_pattern_bgcor_noisefilled=single_pattern_bgcor;
single_pattern_bgcor_noisefilled( ...
    single_pattern==0)=noiseMask(single_pattern==0);

% deal with the uneveness my dividing by a Gaussian mask
single_pattern_bgcor_noisefilled = single_pattern_bgcor_noisefilled ./ imgaussfilt(single_pattern_bgcor_noisefilled,10);

% correct again
single_std=nanstd(single_pattern_bgcor_noisefilled(:));
single_mean=nanmean(single_pattern_bgcor_noisefilled(:));
noiseMask=randn(size(single_pattern_bgcor))*single_std+single_mean;
single_pattern_bgcor_noisefilled( ...
    single_pattern==0)=noiseMask(single_pattern==0);

figure; 
subplot(1,4,1); imagesc(single_pattern); axis xy; axis image; colormap('gray');
subplot(1,4,2); imagesc(background); axis xy; axis image; colormap('gray');
subplot(1,4,3); imagesc(single_pattern_bgcor); axis xy; axis image; colormap('gray');

% histogram stretch
single_pattern_bgcor_noisefilled = normalizeto1(single_pattern_bgcor_noisefilled);
single_pattern_bgcor_noisefilled = imadjust(single_pattern_bgcor_noisefilled,[mean(single_pattern_bgcor_noisefilled(:))-2*std(single_pattern_bgcor_noisefilled(:)), mean(single_pattern_bgcor_noisefilled(:))+2*std(single_pattern_bgcor_noisefilled(:))],[]);
subplot(1,4,4); imagesc(single_pattern_bgcor_noisefilled); axis xy; axis image; colormap('gray');

% FFT
patternFFT = log10(abs(fftshift( fft2(single_pattern_bgcor_noisefilled))).^2);

figure;
imagesc(patternFFT); axis xy; axis image; colormap('gray'); axis off; clim([2 6]);

%% Write to image files
single_pattern_16 = flipud(normalizeto16bit(single_pattern));
single_pattern_bgcor_16 = flipud(normalizeto16bit(single_pattern_bgcor));
single_pattern_bgcor_noisefilled_16 = flipud(normalizeto16bit(single_pattern_bgcor_noisefilled));

outputFormat = '.png';
patIdentifier = strcat('RKD_x',int2str(x_val),'_y',int2str(y_val),'_');
outputRawPatName = strcat(patIdentifier,'raw',outputFormat);
outputBGCorPatName = strcat(patIdentifier,'bgcor',outputFormat);
outputBGCorNoiseFilledPatName = strcat(patIdentifier,'bgcor_filled',outputFormat);

% uncomment the following to write to files
% imwrite(single_pattern_16,fullfile(TFS_DataLoc, outputRawPatName));
% imwrite(single_pattern_bgcor_16,fullfile(TFS_DataLoc, outputBGCorPatName));
% imwrite(single_pattern_bgcor_noisefilled_16,fullfile(TFS_DataLoc, outputBGCorNoiseFilledPatName));
