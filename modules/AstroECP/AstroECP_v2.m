%% AstroECP v2
%{ 

Qaiser, M. H., Berners, L., Scales, R. J., Zhang, T., Heller, M., DluhoÅ¡, J., Korte-Kerzel, S. 
& Britton, T. B. (2026). AstroECP: towards more practical electron channeling contrast imaging. 
J. Appl. Cryst. 59.  https://doi.org/10.1107/S1600576726000567  
 -------
 a digital twin of SEM within a graphical user interface to help navigate the electron channeling patterns 
 and quickly access specific ECCI conditions for subsequent analysis, based on the approaches and functions used in 
 AstroEBSD (Britton, Tong et al., 2018) and MTEX v5.10.2 (Bachmann et al., 2010; Hielscher et al., 2019). 
 AstroECP reads the dynamically calculated Kikuchi reference patterns generated from pattern simulation software 
 that include: MapSweeper EBSD dynamical patterns (higher quality Bloch Wave Kikuchi Diffraction (BWKD) patterns, 
 the open-source software EMsoft (Singh, Marc de Graefe et al., 2017) and also Bruker DynamicS. The GUI then 
 reprojects the channeling pattern based upon the crystal orientation and ECP conditions. Manipulation of the 
 simulation can be performed through virtual movement of the stage and microscope to inform movement of the sample 
 within the SEM. This helps to effectively orient the crystal for the development of crystallographic contrast with 
 a specific [uvw] direction along the optic axis of the SEM, and provide imaging conditions that are suitable for 
 the collection of ECCI micrographs.
 -------- 
%}

clear; home; close all;

%% set up the AstroEBSD and MTEX folder locations 
% Input_Data.astro_location='D:\OneDrive - UBC\Documents\GitHub\AstroEBSD_v2\';
Input_Data.astro_location='D:\OneDrive - UBC\Documents\GitHub\AstroEBSD\';
Input_Data.mtex_location='D:\OneDrive - UBC\MatLab\mtex-5.11.2\'; %working with 5.11.2

%% Example data of Si ECP. Comment out this block if you are running your own experimental data.

% Input_Data.image_folder=[Input_Data.astro_location '\modules\AstroECP'];    % example pattern location
% Input_Data.image_name='Si_SAECP_example.tif';   % example pattern name. available in https://github.com/ExpMicroMech/AstroEBSD/blob/main/modules/AstroECP/    
% addpath (Input_Data.image_folder);
% Input_Data.image_frame=1;       % frame number for Tescan data, this variable should not exist for other file types
% Input_Data.ECP_type='TESCAN';    % supported types: 'TESCAN', 'CLARA', 'TFS', 'other'. The TESCAN option has double stage Tilt 
% Input_Data.PC_in=[0.5 0.5 3.8]; % starting PC - AstroEBSD convention [PCx, PCy, DD]
% Input_Data.eangs=[177.330,16.558,-132.238]; % for the example pattern
% Input_Data.ECP_Pat_clim=[5 12]; % default settings of histogram
% Input_Data.Stage_in=[0 0 0]; % stage rotation settings, in degrees [Rx, Ry, Rz]

%% Your experimental data. Comment out this block if you want to run Example data from the above block

Input_Data.image_folder=['D:\OneDrive - UBC\Si_ECP_Collection'];          % your experimental pattern location
Input_Data.image_name='Si_ECP_WD4mm_2.tif';                               % your experimental pattern name

addpath (Input_Data.image_folder);

Input_Data.image_frame=1;       % frame number for TESCAN and CLARA data, this variable should not exist for other file types
Input_Data.ECP_type='CLARA';    % supported types: 'TESCAN', 'CLARA', 'TFS', 'other'. The TESCAN option has double stage Tilt 

Input_Data.PC_in=[0.5 0.5 2.8]; % starting pattern center - AstroEBSD convention [PCx, PCy, DD]
Input_Data.eangs=[0,0,0]; % euler angles in degrees

Input_Data.ECP_Pat_clim=[2 7]; % set the greyscale colorlim for the ECP if needed
Input_Data.Stage_in=[0 0 0];   % stage rotation settings, in degrees [Rx, Ry, Rz]

%% Crystal information

Input_Data.Phase_Input{1}='Si_20kV'; % pick the phase file - this should be placed in \AstroEBSD\phases\phasefiles
Input_Data.crystal_shape='cube'; %the crystal shape to use from MTEX - 'hex' and 'cube' coded or 'orthorhombic'
Input_Data.miller1=[0 0 1]; %miller indicies for PF plot of the crystal
Input_Data.miller2=[1 1 0]; %miller indicies for PF plot of the crystal

%% other default settings 
Input_Data.Phase_Folder = fullfile(Input_Data.astro_location,'phases'); %location of the AstroEBSD phases super-folder
Input_Data.Bin_loc = fullfile(Input_Data.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM
Input_Data.Index=1; %add the indexing button to the GUI

%% Start MTEX and AstroEBSD
try EBSD;
catch
    run(fullfile(Input_Data.mtex_location,"startup_mtex.m"));
end

try astro_loadcheck
catch
    run(fullfile(Input_Data.astro_location,"start_AstroEBSD.m"));
end

%% set some MTEX preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');

%% Load data and run the GUI

% Experimental ECP pattern load
[ECP_Pat] = AstroECP_Load(Input_Data);

%% Optional filtering 
% % background correct the ECP
% Settings_CorX.gfilt=1; %use a filter (1 or 0)
% Settings_CorX.gfilt_s=2; %low pass filter sigma
% Settings_CorX.radius=1; %use a radius mask
% Settings_CorX.radius_frac=0.65; %fraction of the pattern width to use as the mask
% ECP_Pat.pattern_BG_original=ECP_Pat.ECP_Pat_BG;
% ECP_Pat.ECP_Pat_BG=EBSP_BGCor( ECP_Pat.ECP_Pat_BG,Settings_CorX);
% ECP_Pat.pattern=ECP_Pat.ECP_Pat_BG;
% Input_Data.ECP_Pat_clim=[-0.05 0.05]; % default settings of histogram

%% Run the AstroECP GUI
[Output_Data]=fAstroECP(Input_Data,ECP_Pat);

