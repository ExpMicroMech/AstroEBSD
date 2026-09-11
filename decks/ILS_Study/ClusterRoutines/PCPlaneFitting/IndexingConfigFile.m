% Author: Lukas Berners
% Date: 23.07.2025
% Description: Setup configuration for the hough based indexing and further
% refinement of ebsd patterns based on initial orinetation from the
% manufacturer
%
% This script is provided "as is" without any warranty of any kind.
% The author shall not be liable for any damages arising from the use 
% of this software. Use at your own risk.
InputUser.MTEX_loc='/work/mz071159/mtex';
% InputUser.MTEX_loc='C:\Users\berners\MTEX_dev\mtex';
try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end


InputUser.Phase_Folder="/work/mz071159/AstroEBSD_v2/phases";
InputUser.Astro_loc='/work/mz071159/AstroEBSD_v2';

try
    ASTROSKRIPTS_LOC()
catch
    warning('Astroskirpts_loc not found')
end
addpath("/work/mz071159/astro_skripts")
try
    ASTROSKRIPTS_LOC()
catch
    warning('Astroskirpts_loc not found')
end


try astro_loadcheck
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end

RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
RTM.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates');
InputUser.RTM = RTM;

%% colors
InputData.cset=[1,0,0; % Colors of the Bands in rgb (1 =255)
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0]
InputData.cset = tab10;

 [cmap_div]=cbrewer('div','RdBu',30);
%  [cmap_div] = cbrewer('')
cmap_div(cmap_div<0)=0;


%% File Variables
% InputUser.Phase_Folder = Phase_Folder;
% InputUser.Astro_loc = Astro_loc;
% InputUser.OSC_folder = 'J:\Users\Berners_Lukas\clara\FeTaAl_Laves_Christina\T1_4_CG\FeTaLaves__T1.4__Area 2__EBSD+EDS Map 3\_EBSD\vt0,1'
% InputUser.OSC_folder = 'E:\sciebo\IMM\EBSD'
% InputUser.OSC_folder = input_folder
% InputUser.OSC_file = 'LaNi2_LaNi2 vibropol_Area 6_EBSD+EDS Map 1.ang'
% InputUser.OSC_folder = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\Nb_49Ni_V2_spetemberbatch'
InputUser.OSC_folder = "/hpcwork/mz071159/EBSD/__EBSDFOLDER__";
InputUser.OSC_file = '__EBSDMAP__.ang';
% InputUser.OSC_file = 'FeTaLaves_T1.4_Area 2_EBSD+EDS Map 3.ang'
InputUser.HDF5_file = InputUser.OSC_file;
InputUser.HDF5_folder = InputUser.OSC_folder;
%[ MapInfo.MapData,MicroscopeData,~,~ ] = breadh5oina(InputUser); %change input user, if the file moved
% InputUser.Phase_Input={'Nb6Ni7_em','Nb7Ni6_em','Nb57_Ni43_joubert','Nb51_Ni49_joubert'};%,'Nb57_Ni43_joubert_c1-8_c2-16'}Mg_20kv_em'}%'
InputUser.Phase_Input = {__PHASENAME__}
subFolder =__PHASEFOLDER__
if ~strcmp(subFolder,"")

    PhasesToRead = dir(fullfile(InputUser.Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
    InputUser.Phase_Input = cell(size(PhasesToRead));
    for i=1:numel(PhasesToRead)
        InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
    end
end
%'Nb518Ni482_joubert_20_2025_07_31__21_42_11_dmin-0.020_c1-4.000_c2-8.000_ehistmin-19.800_npx-1000',...
  %  'Nb518Ni482_joubert_18_2025_08_01__11_26_24_dmin-0.020_c1-5.000_c2-10.000_ehistmin-17.800_npx-1000',...
    % 'Nb518Ni482_joubert_19_2025_08_01__11_26_28_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
    % 'Nb518Ni482_joubert_21_2025_08_01__11_26_32_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.800_npx-1000',...
    % 'Nb518Ni482_joubert_19_2025_08_01__11_26_29_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.600_npx-1000',...
    % 'NbNi_dt0.400_c26.5_a4.9_20kV_2025_07_21__16_00_52_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'NbNi_dt0.400_c26.9_a4.9_20kV_2025_07_21__16_01_10_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'NbNi_dt0.400_c26.7_a4.9',...
    % 'NbNi_dt0.250_c26.5_a4.9_20kV_2025_07_22__10_59_34_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'Nb569Ni431_joubert_20kV_2025_07_24__08_56_35_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
    % %  'Nb553Ni447_joubert_20kV_2025_07_24__08_56_34_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'Nb530Ni470_joubert_20kV_2025_07_24__08_56_34_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'Nb518Ni482_joubert_20kV_2025_07_24__08_56_33_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % 'Nb496Ni504_joubert_20kV_2025_07_24__08_56_31_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
    % }
% InputUser.Phase_Input = {'La7Ni16_30kV'}
InputUser.ExportPath = fullfile(InputUser.OSC_folder,InputUser.OSC_file(1:end-4), "PatternCentrePlane");
mkdir(InputUser.ExportPath)

appendix_string=subFolder;


%% Crystal Structure Variables
InputUser.Phases = InputUser.Phase_Input;
[~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(1),InputUser.Phase_Folder );
[cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);

%
cs = crystalAlignEdax(cs);
InputUser.CS = {'notIndexed',cs};
InputUser.AlignEdax = 1;
InputUser.EdaxConversionMatrix = RTM.Ry(180*degree);
% InputUser.CS=cs
[ MapInfo.MapData,MicroscopeData,~,MapInfo.EBSPData ] = bReadOSC(InputUser); 
% MapInfo.MapData.g_num=MapInfo.g_num;
ebsd_map = MapInfo.EBSPData.ebsd_map;


%%

%% set up diffraction bits.

InputData.V_in=29E3; %voltage in eV
InputData.lattice_param_scale=1E-10;
[InputData.lambda] = calc_lambda(InputData.V_in);
InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale;

InputData.maxBNum = 6;
InputUser.maxInd = 30;
InputUser.cutoff = 0.05;
InputUser.minD = 0.2;

InputUser.AlignEdax =1 ;

% [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{1}}, InputData.lambda_p ,InputUser);

% cutoff = 0.2;
% [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,cutoff);

% m_plot = m_plot(1:10);
% I = I(1:10);
% InputData.IPlot = I;

InputUser.calibration = loadCalibrationEDAX( fullfile( InputUser.OSC_folder,InputUser.OSC_file) );
InputUser.calibration.DetPXSize = 55 % (µm, for EDAX Clarity. For other detector please change.
InputUser.calibration.DetSize = InputUser.calibration.DetPXSize*514 % 
InputUser.calibration.DetNumPixels = 514;
InputUser.PC = [InputUser.calibration.pcx,InputUser.calibration.pcy,InputUser.calibration.pcz];
[pcx,pcy, pcz] = patternCentreGeometric(InputUser.calibration,MapInfo.MapData);
MapInfo.MapData.pcx = pcx;
MapInfo.MapData.pcy = pcy;
MapInfo.MapData.pcz = pcz;
MicroscopeData.SampleTilt = InputUser.calibration.SampleTilt;
MicroscopeData.CameraTilt = InputUser.calibration.CameraElevationAngle;
MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180;

MicroscopeData.rot_edax_conversion=RTM.Ry(180*degree);

MicroscopeData.Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

%% Get an ebsd pattern to collect sizing information
pnum=1;
pat=bReadEBSP(MapInfo.EBSPData,pnum);


%% correction settings
% Settings_Cor_rtm=Settings_Cor;
Settings_Cor.gfilt=1;
Settings_Cor.gfilt_s=4;
Settings_Cor.radius=0; %use a radius mask
Settings_Cor.radius_frac=0.65;
Settings_Cor.resize=0;
Settings_Cor.Square=0; %square crop (1)
Settings_Cor.SquareCrop=1;
Settings_Cor.MeanCentre=1;
Settings_Cor.gaussian_range = 2; % try to increas a bit.
Settings_Cor.quadrants = 1;%quadrants
Settings_Cor.cross = 3
Settings_Cor.gfiltQuad=1
SettingsCor.hotpixel=1
Settings_Cor.hot_thresh=1000


pattern_info.size=size(pat);
pattern_info.size=uint32(pattern_info.size);
Settings_Cor.size=pattern_info.size;

Settings_CorDictionary = Settings_Cor;
Settings_Cor.resize=1 ;
Settings_Cor.size = uint32(pattern_info.size);

Settings_Cor.max_var_pc=0.1;
Settings_Cor.max_var_ori=5*degree;
Settings_Cor.quadrants=1;
Settings_Cor.cross=2; %make cross larger
Settings_Cor.gfiltQuad=1;
Settings_Cor.MeanCentre = 1
%%


InputUser.RTM = RTM;


RTM_info=RTM;
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
% RTM_info=struct()
RTM.screensize =double(Settings_Cor.size(1)); %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM.iterations = 6;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 500;
RTM.isHex=1;
% RTM.Rx=RTM.Rx;
% RTM.Rz=RTM.Rz;
% RTM.Ry=RTM.Ry;


RTM.PatternFlip=0;
RTM.apply_gauss=0;
%%
options.SettingsCor = Settings_Cor
% iteration options
options.MaxIter = 300;
options.StepTolerance = 1e-12;
options.OptimalityTolerance = 1e-12;
options.BGVar = 0;
options.BGVarVal = 10;
options.Diagnostics = 'off';
options.Display = 'final';
options.gradients = __USEGRADIENTS__;


gridPoints = __GRIDPOINTS__;
kernelSize =__KERNELSIZE__;
pSerachSize = __PSEARCHSIZE__;
options.fminsearch=0;

scalingFactor =1;