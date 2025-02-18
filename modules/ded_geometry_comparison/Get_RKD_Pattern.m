%% Get_RKD_Pattern.m
% T. B. Britton and T. Zhang 
% October 2024
% This script read RKD data captured by Thermo Fisher Scientific EBSD software), process and output a
% selected RKD pattern

% This script accompanies the following article(s):
% "Comparison of Kikuchi Diffraction Geometries in Scanning Electron Microscopy"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% This script ideally works for smaller datasets with ~400 patterns. For
% much larger datasets, it is encouraged to use Get_RKD_pattern_v2.m which is structured differently.

% Requirements: 
% (1) AstroEBSD package - this script is available as a part of the latest
% AstroEBSSD distribution.
% (Optional) MTEX toolbox (https://mtex-toolbox.github.io/) - plotting the
% full dataset.
% N.B. You may still need other dependencies to run the other scripts in
% AstroEBSD.
%%
clear;
home;
close all;

%% Input Variables

mtex_location='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';
astro_location='D:\AstroEBSD';

% path to data file
TFS_DataLoc='D:\2024-11-04T15_22_24_Al_ben\2024-11-04T16_31_57';
InputUser.Phase_Input={'Al_bwkd'};


%% Load MTEX

%if MTEX loaded this doesn't matter
%if MTEX is not loaded, then this will be used to start MTEX up

%start mtex if needed
try EBSD;
catch
    run(fullfile(mtex_location,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(astro_location,"start_AstroEBSD.m"));
end

RTM.Phase_Folder = fullfile(astro_location,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% Load the EBSD data

% crystal symmetry
% CS = {... 
%   'notIndexed',...
%   crystalSymmetry('432', [3.5 3.5 3.5], 'mineral', 'Nickel', 'color', [0.53 0.81 0.98])};


% multiphase_full=fullfile(multiphase_location,'\results','multiphases.idx');

[ebsd,ebsd_header] = bIDX_to_EBSD(TFS_DataLoc);

% plotting convention - TFS
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
%% create a reference pattern from the Dynamics bin file
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},RTM.Phase_Folder);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Load a pattern from the data stack

% Select pattern position - this can be a random one
x_val=23;
y_val=19;

[single_pattern,row_patterns] = rEBSD_TFS_Cor(y_val,x_val,TFS_DataLoc); %note this is in matlab indexing
RTM.size=[479 479];
[EBSD_geom ] = EBSP_Gnom( RTM,ebsd_header); %you can change PC_in if you want    

%find the value in the original EBSD container for this
ebsd_val=find(ebsd.prop.y_in == y_val & ebsd.prop.x_in == x_val );
G_MTEX=ebsd(ebsd_val).orientations.matrix;    
% Note that the orientation data in the original dataset is indexed by the
% EBSD software and is not accurate. In the paper, we reindexed the
% patterns (both PC and orientation).

[ Pat_sim_A] = EBSP_gen( EBSD_geom,G_MTEX'*RTM.Rz(-15*pi/180),screen_int);

figure; 
subplot(2,2,1);
pPattern(single_pattern,EBSD_geom); title('Experiment')
subplot(2,2,2);
pPattern(Pat_sim_A,EBSD_geom); title('simulation')

Settings_CorX.rkd=1; %RKD geometry - fill gaps for flatfielded RKD quad chip
[ Pat_exp_cor ] = EBSP_BGCor(single_pattern,Settings_CorX); %note the flipud, this is because TFS uses Y pointing down for the pattern
Pat_exp_cor2.size=size(Pat_exp_cor);
[EBSD_geom2 ] = EBSP_Gnom( Pat_exp_cor2,ebsd_header); %you can change PC_in if you want    

subplot(2,2,3);
pPattern(Pat_exp_cor,EBSD_geom2); title('Corrected Experiment')
[ Pat_sim_B] = EBSP_gen( EBSD_geom2,G_MTEX'*RTM.Rz(-15*pi/180),screen_int);
subplot(2,2,4);
pPattern(Pat_sim_B,EBSD_geom2); title('Corrected Simulation')

%% Load the raw pattern data

raw_loc=fullfile(TFS_DataLoc,'\raw');
raw_folders=dir(raw_loc);
raw_folder=fullfile(raw_loc,raw_folders(3).name);
raw_pats=fullfile(raw_folder,'f_00000.pak');

%find the number of grid points in X and Y
x_max=max(ebsd.prop.x_in(:))+1;
y_max=max(ebsd.prop.y_in(:))+1;
pat_num=(x_val+1)+y_val*x_max;

total_num_pats = size(ebsd,1);

% Read the raw stack
[pat_raw_stack] = rEBSD_TFS_Raw(raw_pats, total_num_pats);
av_bg=sum(pat_raw_stack,3)/x_max/y_max;

% Plot patterns
pat_single=pat_raw_stack(:,:,pat_num)./av_bg;
Pat_exp_cor3.size=size(pat_single);
[EBSD_geom3 ] = EBSP_Gnom( Pat_exp_cor3,ebsd_header); %you can change PC_in if you want  

single_std=nanstd(pat_single(:));
single_mean=nanmean(pat_single(:));

EBSP2=randn(size(pat_single))*single_std+single_mean;
pat_single_cor=pat_single;
pat_single_cor(isnan(pat_single))=EBSP2(isnan(pat_single));

figure;
subplot(2,3,1);
pPattern(pat_raw_stack(:,:,pat_num),EBSD_geom3); title('Raw Exp')

subplot(2,3,2);
Pat_exp_cor3.size=size(pat_single);
pPattern(pat_single_cor,EBSD_geom3); title('Raw-Ish Exp')

subplot(2,3,3);
Pat_exp_cor3.size=size(pat_single);
[EBSD_geom3 ] = EBSP_Gnom( Pat_exp_cor3,ebsd_header); %you can change PC_in if you want 

[ Pat_sim_C] = EBSP_gen( EBSD_geom3,G_MTEX'*RTM.Rz(-15*pi/180),screen_int);
pPattern(Pat_sim_C,EBSD_geom2); title('Raw Sim')

subplot(2,3,5);
pPattern(Pat_exp_cor,EBSD_geom2); title('Corrected Experiment')
[ Pat_sim_B] = EBSP_gen( EBSD_geom2,G_MTEX'*RTM.Rz(-15*pi/180),screen_int);
subplot(2,3,6);
pPattern(Pat_sim_B,EBSD_geom2); title('Corrected Simulation')


%% The RKD dataset is a map and you can plot these data using the code below.
%% IPF-Z with unit cells
% scaling = 6; % scale the crystal shape to have a nice size
% mineral=ebsd(1).CS.mineral;
% 
% figure;
% colorKey = ipfTSLKey(ebsd.CS);
% % colorKey = ipfHSVKey(ebsd.CS);
% colorKey.inversePoleFigureDirection = zvector;
% color = colorKey.orientation2color(ebsd('indexed').orientations);
% 
% plot(ebsd,colorKey.orientation2color(ebsd.orientations),'coordinates','on')
% 
% hold on;
% 
% % grain segementation
% [grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree);
% 
% % smooth grain boundaries
% grains = smooth(grains,5);
% 
% grains=grains(grains.area > 20); %remove grains which are less than 20 um^2
% 
% % plot the unit cells
% %plot the crystal shape for this point
% 
% %here we use a hexagon, but you could use a cube for cubic
% %crystalShape.cube
% cS = crystalShape.hex(grains(mineral).CS);
% grain_centroids=grains.centroid;
% plot(grain_centroids(:,1),grain_centroids(:,2),5*scaling, grains.meanOrientation * cS * scaling,'FaceAlpha',0.2); %make the hexagons also transparent so you can see the grains behind
% xposn=ebsd.prop.x;
% yposn=ebsd.prop.y;
% xlim([min(xposn(:)) max(xposn(:))]);
% ylim([min(yposn(:)) max(yposn(:))]);
% 
% nextAxis
% colorKey.inversePoleFigureDirection = xvector;
% color = colorKey.orientation2color(ebsd('indexed').orientations);
% 
% plot(ebsd,colorKey.orientation2color(ebsd.orientations),'coordinates','on')
% 
% nextAxis
% colorKey.inversePoleFigureDirection = yvector;
% color = colorKey.orientation2color(ebsd('indexed').orientations);
% 
% plot(ebsd,colorKey.orientation2color(ebsd.orientations),'coordinates','on')

%% Plot the ODF
% odf = calcDensity(ebsd(mineral).orientations,'halfwidth',10*degree);
% 
% %
% figure;
% % plot the pole figure representation of the ODF
% plotPDF(odf,[Miller(0,0,1,ebsd(mineral).CS) Miller(0,1,0,ebsd(mineral).CS)],'projection','eangle');
% 
% mtexColorbar('location','southoutside')
% mtexColorMap blue2red
% 
% %
% %add contours
% levels=[0:33];
% hold on;
% plotPDF(odf,[Miller(0,0,1,ebsd(mineral).CS) Miller(0,1,0,ebsd(mineral).CS)],'contour',levels,'linecolor','k','projection','eangle');
% mtexColorbar('location','southoutside')

