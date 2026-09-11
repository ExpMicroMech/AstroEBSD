function Settings_CorX = defaultStereoConfig(options)
   
    %   defaultRTMConfig return default config for Refined
    %   rtmConfig = defaultRTMConfig(
    %
    %
    %
    % for i = 1:numel(varargin)
    %     if strcmp('varargin'),'clar'
    % end
    % clarity = 0
     switch nargin
         case 0
             options = struct();
      end
    % end
    
    if ~isfield(options,'clarity')
        options.clarity = 0;
    end
    if ~isfield(options,'bgCor')
        options.bgCor = 0;
    end
    if ~isfield(options,'noMask')
        options.noMask = 0;
    end
    %% Refinement
    % RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    % RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    % RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
    % RTM.screensize = 64; %size of the library patterns and the resize of the raw EBSPs
    % RTM.Sampling_Freq=7; %Set the SO(3) sampling freq in degrees
    % RTM.iterations = 8;%Set the number of iterations to do in the refinement step
    % RTM.LPTsize = 128; %LPT size used in pixels
    % RTM.Phase_Folder = fullfile(ASTRO_LOC(),'phases'); %location of the AstroEBSD phases super-folder
    % RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates');
    % rtmConfig.RTM = RTM

    %% Image Processing
    Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
    Settings_CorX.gfilt_s=4; %low pass filter sigma
    Settings_CorX.radius=0; %use a radius mask
    Settings_CorX.radius_frac=0.65; %fraction of the pattern width to use as the mask
    Settings_CorX.hotpixel=0; %hot pixel correction
    Settings_CorX.hot_thresh=1000; %hot pixel threshold
    Settings_CorX.resize=1; %resize correction
    Settings_CorX.size=514; %image height
    Settings_CorX.RealBG=0; %use a real BG
    Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG
    Settings_CorX.SquareCrop = 0; %make square the EBSP
    Settings_CorX.SplitBG=0; %deal with a split screen
    Settings_CorX.MeanCentre = 1;
    
    % default for optimisation
    Settings_CorX.max_var_pc=0.2;
    Settings_CorX.max_var_ori=10*degree;
    
    if options.noMask==0
    % default for stereographic reprojections
        Settings_CorX.gaussianRange=5;
        Settings_CorX.gaussianRangeLarge =5 ;
        Settings_CorX.disc_sizeLarge = 0.95;
        Settings_CorX.disc_size=0.3;
    end
    % default settings for clarity
    if options.clarity == 1
        Settings_CorX.cross = 2;
        Settings_CorX.quadrants = 1;
        Settings_CorX.gfiltQuad = 1;
    end
    
    % if options.bgCor == 1
    %     Settings_CorX.gfilt = 1;
    %     Settings_CorX.gfilt_s = 4;
    % end

    rtmConfig.Settings_CorX = Settings_CorX;
end