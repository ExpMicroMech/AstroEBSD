%% Gen_Sim_Patterns.m
% This script reprojects a Kikuchi pattern for a given PC and orientation
% (Euler angles). This script does not help with pattern
% indexing/matching.
% This script can run on its own, or be called by SinglePatchEBSP.m. 

% requirement:
% AstroEBSD
% MTEX toolbox

% Tianbi Zhang and T. Ben Britton, April 2025

%% load astro and MTEX
InputUser.Astro_loc='D:\AstroEBSD_v2_standalone\';
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

opt.output = 0;

%% Analyse a single pattern

%build the phases
InputUser.Phase_Folder = fullfile(InputUser.Astro_loc,'phases');
InputUser.Phase_Input  = {'Si_bwkd'}; %Si, Ferrite
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
%% Set up the RTM

RTMgensim.screensize = 256; %size of the simulated pattern(s)
RTMgensim.Sampling_Freq=5; %Set the SO(3) sampling freq in degrees
RTMgensim.iterations = 100;%Set the number of iterations to do in the refinement step
RTMgensim.LPTsize = 128; %LPT size used in pixels

%% Low level setting up stuff - you shouldn't need to change this
RTMgensim.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTMgensim.Bin_loc = fullfile(RTMgensim.Phase_Folder,'masterpatterns'); %location of the binary files used for RTM

[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTMgensim.screensize, RTMgensim.LPTsize );

%Define all rotation matrices needed in the code
RTMgensim.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTMgensim.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTMgensim.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Generate a simulated pattern

eangs1 = [ 45 0 0 ] * degree; %
Settings_PCin.start=[ 0.4881043 0.3687062 0.4821355 ];

[EBSD_geom_gensim ] = EBSP_Gnom( RTMgensim,Settings_PCin.start); 

Detector_tilt=eye(3); %

gmatrix1=RTMgensim.Rz(eangs1(3))*RTMgensim.Rx(eangs1(2))*RTMgensim.Rz(eangs1(1));

[ Pat_sim_eang1 ] = EBSP_gen( EBSD_geom_gensim,gmatrix1*Detector_tilt,screen_int); 

Settings_Cor.MeanCentre=1;
[Pat_sim_eang1, ~] = EBSP_BGCor(Pat_sim_eang1, Settings_Cor);

figure;
pPattern(Pat_sim_eang1, EBSD_geom_gensim);axis xy; axis image; colormap('gray'); axis off; 

% pat_sim_out_16 = normalizeto16bit(Pat_sim_eang1);
% imwrite(flipud(pat_sim_out_16), "1.tif");

%%
function normalized_matrix = normalizeto16bit(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
normalized_matrix = uint16(normalized_matrix * (2^16 - 1)); % Convert to 16 bit
end