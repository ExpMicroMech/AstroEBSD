function [G_Refined_SO3,PH_SO3] = ECP_Index(ECP_Pat,Input_Data,PC_in)
%ECP_INDEX Summary of this function goes here
%   Detailed explanation goes here

% play with these
RTM.screensize = 128; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=4; %Set the SO(3) sampling freq in degrees - probably should be 1/2 the ECP subtended angle, if this is very small - the cost (RAM and CPU) is high
RTM.iterations = 3;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 128; %LPT size used in pixels
RTM.singleprecision= 1; %use single precision (1) or double (0)

%background correction settings
Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=5; %low pass filter sigma
Settings_CorX.resize=1; %resize correction
Settings_CorX.size=[RTM.screensize RTM.screensize]; %image height
Settings_CorX.SquareCrop=1;
%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% RTM Loading - phase + dynamic pattern, and generate the RTM filter settings
%load the crystallographic informatio
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {Input_Data.Phase_Input{1}},Input_Data.Phase_Folder); %#ok<CCAT1>
%build the look up table for the dynamical patterns
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex,'resize',512);

%set up the initial filter settings for the RTM approach
[ SettingsXCF] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );
SettingsXCF.single=RTM.singleprecision; %make the library and patterns single precision

%generate the library
[library_G,template_library] = ECP_LibraryGen(screen_int,RTM_info,RTM,SettingsXCF,PC_in);

%match
[G_Refined_SO3,PH_SO3]=ECP_LibraryMatch(ECP_Pat,template_library,library_G,screen_int,Settings_CorX,SettingsXCF,RTM,PC_in);

% %plot, for debugging
% [EBSD_geom ] = EBSP_Gnom( RTM, PC_in); %
% [ Pat_sim_SO3_refined ] = EBSP_gen( EBSD_geom,G_Refined_SO3,screen_int); %generate the EBSP for this iteration 
% figure;
% pPattern(Pat_sim_SO3_refined,EBSD_geom);
end