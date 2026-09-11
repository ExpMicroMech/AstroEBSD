function rtmConfig = defaultRTMConfig(options)
   
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
    if ~isfield(options,'refine')
        options.refine = 0;
    end
    %% Refinement
    RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
    if options.refine ==1
        RTM.screensize = 256; %size of the library patterns and the resize of the raw EBSPs
        rtmConfig.Refine.num_pts=6; %number of points to use for the refinement step
        rtmConfig.Refine.ss=0.2; %initial probe volume
        rtmConfig.Refine.p=2; %order of polynomial to fit to tetrahedron
        rtmConfig.Refine.n_its=10; %number of interations
        rtmConfig.Refine.debug=1; 
        RTM.LPTsize = 256; 
        % Settings_CorX.resize = 1;
    else
        RTM.screensize = 64;
        RTM.LPTsize = 128; 
    end
    RTM.Sampling_Freq=7; %Set the SO(3) sampling freq in degrees
    RTM.iterations = 8;%Set the number of iterations to do in the refinement step
%LPT size used in pixels
    RTM.Phase_Folder = fullfile(ASTRO_LOC(),'phases'); %location of the AstroEBSD phases super-folder
    RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates');


    %% Image Processing
    Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
    Settings_CorX.gfilt_s=4; %low pass filter sigma
    Settings_CorX.radius=0; %use a radius mask
    Settings_CorX.radius_frac=0.65; %fraction of the pattern width to use as the mask
    Settings_CorX.hotpixel=0; %hot pixel correction
    Settings_CorX.hot_thresh=1000; %hot pixel threshold
    Settings_CorX.resize=1; %resize correction
    Settings_CorX.size= RTM.screensize; %image height
    Settings_CorX.RealBG=0; %use a real BG
    Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG
    Settings_CorX.SquareCrop = 1; %make square the EBSP
    Settings_CorX.SplitBG=0; %deal with a split screen
    Settings_CorX.MeanCentre = 1;
    
    % default for optimisation
    Settings_CorX.max_var_pc=0.2;
    Settings_CorX.max_var_ori=5*degree;
    Settings_CorX.screensize = Settings_CorX.size;
    % default settings for clarity
    if options.clarity == 1
        Settings_CorX.cross = 2;
        Settings_CorX.quadrants = 1;
        Settings_CorX.gfiltQuad = 1;
        RTM.EdaxConversionMatrix = RTM.Ry(180*degree); % tested for the clarity @ IMM, RWTH Aachen

    else
    
    RTM.EdaxConversionMatrix = RTM.Rx(0*degree); % unit matrix
    end
    
    rtmConfig.Settings_Cor = Settings_CorX;
    rtmConfig.generalOptions.AlignEdax = 1;
    rtmConfig.RTM = RTM;
end