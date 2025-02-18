%% Reproject_To_Pattern
% Created by Tianbi Zhang in November 2024

% This script generates a gnomonic Kikuchi pattern at a known PC and orientation 
% from a diffraction sphere, 
% which can be in the form of cubic interpolation (Bruker Espirit DynamicS),
% modified Lambert projection (EMsoft), stereograms (Oxford Instruments MapSweeper/BWKD
% by Winkelmann et al.)
% Check /phases/phasefiles and phases/dynamic_templates on how to structure
% such files.

% This script accompanies the following article(s):
% "Comparison of Kikuchi Diffraction Geometries in Scanning Electron Microscope"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% (In this manuscript, the reprojections from reconstructed stereograms 
% are obtained by (1) stacking two stereograms together to get the dynamic template; 
% (2) construct the necessary files following the bwkd structure; 
% (3) reprojection [this code]. )

% Requirements: 
% (1) AstroEBSD package - this script is available as a part of the latest
% AstroEBSSD distribution.
% (2) MTEX toolbox (https://mtex-toolbox.github.io/)
% N.B. You may still need other dependencies to run the other scripts in
% AstroEBSD.

%%
clear
home
close all

%% load astro and MTEX
InputUser.Astro_loc='D:\AstroEBSD\';
InputUser.mtex_location='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';

%start mtex if needed
try EBSD;
catch
    run(fullfile(InputUser.mtex_location,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end

%% Set up the phase

%build the phases
InputUser.Phase_Folder = fullfile(InputUser.Astro_loc,'phases');
InputUser.Phase_Input  = {'Al_bwkd'}; %Si, Ferrite...
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );


%% Load an experimental pattern

path = 'C:\Users\billy\OneDrive - UBC\PhD\TKD\Kikuchi_Geometry_Comparison\Al_offaxis_tkd\';

patternName = 's15_100f_002s_bgcor.tif';
expPatPath = fullfile(path,patternName);
expPat = flipud(double(imread(expPatPath)));


[RTMexp.screensize,~] = size(expPat);

%size of the simulated pattern(s) - 
% if you want to do a difference plot, please set it to the dimension of the experimental pattern
RTMgensim.screensize = RTMexp.screensize; 

[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
%% Low level setting up stuff - you shouldn't need to change this
RTMgensim.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTMgensim.Bin_loc = fullfile(RTMgensim.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM

% [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTMgensim.screensize, RTMgensim.LPTsize );

%Define all rotation matrices needed in the code
RTMgensim.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTMgensim.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTMgensim.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

RTMexp = RTMgensim;
%% Generate a simulated pattern: PC and orientation

eangs1 = [ 139.7766 53.1139 227.3212 ] * degree;
Settings_PCin.start=[0.45425 -0.15789 0.66814];

[real_geom] = EBSP_Gnom( RTMexp,Settings_PCin.start);
[EBSD_geom_gensim ] = EBSP_Gnom( RTMgensim,Settings_PCin.start); % gnomonic coordinate - needed for band annotation

Detector_tilt=eye(3); % no detector tilt

gmatrix1=RTMgensim.Rz(eangs1(3))*RTMgensim.Rx(eangs1(2))*RTMgensim.Rz(eangs1(1));

% Simulated pattern #1
[ Pat_sim_eang1 ] = EBSP_gen( EBSD_geom_gensim,gmatrix1*Detector_tilt,screen_int); %generate the EBSP for this iteration 

%% Simulate a second pattern, e.g. from reprojection
InputUser.Phase_Input  = {'Al_bwkd'}; %Si, Ferrite
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
[ Pat_sim_eang2 ] = EBSP_gen( EBSD_geom_gensim,gmatrix1*Detector_tilt,screen_int); %generate the EBSP for this iteration 

%% Contrast normalization
Pat_sim_eang1 = Pat_sim_eang1 - mean(Pat_sim_eang1(:)); 
Pat_sim_eang1 = Pat_sim_eang1 / std(Pat_sim_eang1(:));

Pat_sim_eang2 = Pat_sim_eang2 - mean(Pat_sim_eang2(:)); 
Pat_sim_eang2 = Pat_sim_eang2 / std(Pat_sim_eang2(:));

expPat = expPat - mean(expPat(:)); 
expPat = expPat / std(expPat);

%% figure 1: exp, sim1, sim2
subplot(1,3,1); % experimental
pPattern(expPat, real_geom); axis image; axis xy; colormap('gray'); axis off;
subplot(1,3,2); % simulated #1
pPattern(Pat_sim_eang1, EBSD_geom_gensim);axis xy; axis image; colormap('gray'); axis off;
subplot(1,3,3); % simulated #2
pPattern(Pat_sim_eang2, EBSD_geom_gensim);axis xy; axis image; colormap('gray'); axis off;

%% figure 2: Experimental pattern with band overlay
figure;
s1=subplot(1,2,1); % experimental
Plot_EBSPAnnotated( Pat_sim_eang1,EBSD_geom_gensim,[],gmatrix1,Crystal_UCell{1},Crystal_Family{1},s1 );
s1=subplot(1,2,2); % simulated #1
Plot_EBSPAnnotated( expPat,real_geom,[],gmatrix1,Crystal_UCell{1},Crystal_Family{1},s1 );

%% figure 3: FFT of a pattern

figure;
subplot(1,2,1);
pPattern(expPat, real_geom); axis image; axis xy; colormap('gray');
subplot(1,2,2);
imagesc(log10(abs(fftshift( fft2(Pat_sim_eang1))).^2));colormap('gray');axis image;axis xy; axis off;colormap('gray');

%% figure 4: difference 

figure;
imagesc((expPat-Pat_sim_eang1));colormap('gray');axis image;axis xy; axis off;colormap('gray');

% To get the blue-white-red colormap, we used the following MATLAB add-on:
% Adam Auton (2024). Red Blue Colormap 
% (https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap), 
% MATLAB Central File Exchange. Retrieved November 6, 2024.
colormap(redblue(256));

%% Export a pattern
% normalize to 16 bit (or any bit depth as needed, just adjust the function
% below) so we can get a meaningful image file.
Pat_sim_eang1_16 = normalizeto16bit(Pat_sim_eang1);

% Write - set format in the file name.
% imwrite(Pat_sim_eang1_16, "output.png");

function normalized_matrix = normalizeto16bit(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
normalized_matrix = uint16(normalized_matrix * (2^16 - 1)); % Convert to 16 bit
end