%% The CODE  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Some functionallity was added by L. Berners in 2024
% There is a "grained" pca mode, which clusters the orientations and does
% not tile spatially but by grain 
%   set PCA_Setup.Grain_info=1, if you want to use this.    
%      MapInfo.quantil=PCA_Setup.quantil (here you can specify the top
%      values of imagequality to be taken into account if wanted.
%      %% how big is the orientation spread ? : MapInfo.ori_lim=PCA_Setup.ori_lim
%    %% how many paatterns per each orientation do we want to run the pca on: MapInfo.max_pat_num=PCA_Setup.max_pat_num;
%
%
%
%
%
%
%



% %% uncommment following lines if run in batch mode
% run(fullfile(InputUser.MTEX_loc,"startup.m"));
% % % % % % % % % % %cd to the HDF5 folder
% cd(InputUser.HDF5_folder)
%%

% initialise orientation matrices
RTM_setup.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM_setup.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM_setup.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
RTM_setup.Phase_Folder=InputUser.Phase_Folder;
RTM_setup.Bin_loc = fullfile(RTM_setup.Phase_Folder,'masterpatterns');
% Turn off warnings for breadHDF5
InputUser.breadhdf5_warnings_off=1;

if strcmpi(InputUser.Mode,'Isolated')
    Settings_Mode = 1;
elseif strcmpi(InputUser.Mode,'Folder')
    Settings_Mode = 2;
elseif strcmpi(InputUser.Mode,'Map_Single')
    Settings_Mode = 3;
elseif strcmpi(InputUser.Mode,'Map_All')
    Settings_Mode = 4;
else
    error(['The input setting mode ' InputUser.Mode ' is not supported']);
end

    % PCA_Setup.Mode=Settings_Mode

if ~isfield(PCA_Setup,'Grain_info')
    PCA_Setup.Grain_info=0
end

if ~isfield(PCA_Setup,'adjustTiling')
    PCA_Setup.adjustTiling = 0;
end

if Settings_Mode ==  2
	[P_list,num_patterns]=Folder_Prep(InputUser);
    PCA_Setup.P_list=P_list;
    RTM_setup.P_list=P_list;
end
PCA_Setup.Settings_Mode=Settings_Mode;
PCA_Setup.PatternFlip=InputUser.PatternFlip;
%%
%% Sets up required plugins, loads data and calculates cropping windows.
for fileloop=1:length(files)

clear tile RTM Plots EDSWindowing
%%
RTM=RTM_setup;
%%
InputUser.HDF5_file = files{fileloop};
InputUser.OSC_file = files{fileloop};  %h5 file, include the file extension
t1=clock;
if exist(InputUser.SavingFolder) == 0
    mkdir(InputUser.SavingFolder);
end

cd(InputUser.SavingFolder);

%% Creates appropriate filename and HDF5 handling parameters.
PCA_filenamesetup;

%% Load metadata
[~,~,ext]=fileparts(InputUser.HDF5_file)
if strcmp(ext,'.osc')||strcmp(ext,'.ang')

    [MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadOSC(InputUser );

elseif strcmp(ext,'.h5oina')
    [MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=breadh5oina(InputUser);
else
    [MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadHDF5(InputUser );
end
%%
[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);

if isfield(Settings_Cor, 'RealBG')
    if Settings_Cor.RealBG==1
        if ~isfield(Settings_Cor,'EBSP_bg')
            Settings_Cor.EBSP_bg = EBSP_getBackground(MapInfo.EBSPData,100);
        end
    end
end

MapInfo.MapSize1=max(MapInfo.MapData.YBeam)-min(MapInfo.MapData.YBeam)+1; % number of rows
MapInfo.MapSize2=max(MapInfo.MapData.XBeam)-min(MapInfo.MapData.XBeam)+1; % number of columns

%%
if isfield(InputUser,'SampleTilt') && isfield(InputUser,'CameraTilt')
    MicroscopeData.SampleTilt=InputUser.SampleTilt;
    MicroscopeData.CameraTilt=InputUser.CameraTilt;
end
%%
if isfield(MicroscopeData,'CameraTilt')

MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180;
else
    MicroscopeData.TotalTilt=(90-MicroscopeData.SampleTilt)*pi/180;
end
   
if isfield(MapInfo.Data_InputMap,'RadonQuality')
     Plots.Radon=MapInfo.Data_InputMap.RadonQuality;
elseif isfield(MapInfo.Data_InputMap,'ImageQualtiy')
     Plots.Radon=MapInfo.Data_InputMap.ImageQuality;
elseif isfield(MapInfo.Data_InputMap,'PatternQuality')
     Plots.Radon=MapInfo.Data_InputMap.PatternQuality;
else
     Plots.Radon=transpose(ones(MapInfo.MapSize1,MapInfo.MapSize2));
end
if isfield(MapInfo.Data_InputMap,'semsignal')
     Plots.SEMImage=MapInfo.Data_InputMap.semsignal(:,:,1);
elseif isfield(MicroscopeData,'SEMImage')
     Plots.SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));
else
     Plots.SEMImage=transpose(ones(MapInfo.MapSize2,MapInfo.MapSize1));
end
% Plots.SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));

MapInfo.MapSize1=max(MapInfo.MapData.YBeam)-min(MapInfo.MapData.YBeam)+1; % number of rows
MapInfo.MapSize2=max(MapInfo.MapData.XBeam)-min(MapInfo.MapData.XBeam)+1; % number of columns

%%
[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);
% disp('line 103')
% if exist('P_list')
  if Refine.run==1
    [MapInfo]=fRefineMap(Refine,RTM,MapInfo,MicroscopeData,InputUser,Settings_Cor,t1);
end
% if strcmp(ext,'.osc')||strcmp(ext,'.ang')
%     [Pat_inf] = ReadEBSDFile(P_list{1},InputUser.PatternFlip);
%     input_pat_size=size(Pat_inf);
%     % disp('line 106')
%     MapInfo.EBSPData.PW=input_pat_size(1);
%     MapInfo.EBSPData.PH=input_pat_size(2);
%     % disp('line 109')
%     MapInfo.Data_InputMap.PCX=zeros(size(MapInfo.Data_InputMap.XBeam_Map));
%     MapInfo.Data_InputMap.PCY=zeros(size(MapInfo.Data_InputMap.XBeam_Map));
%     MapInfo.Data_InputMap.DD=zeros(size(MapInfo.Data_InputMap.XBeam_Map));
%     MapInfo.Data_InputMap.PCX(:)=InputUser.PC_simulation(1);
%     MapInfo.Data_InputMap.PCY(:)=InputUser.PC_simulation(2);
%     MapInfo.Data_InputMap.DD(:)=InputUser.PC_simulation(3);
% else
% end

%% Refine the pattern centre if requested


%% Plot an example pattern
figure;
if Settings_Mode == 2
%     P_test=randi([1 numel(P_list)])
%     [P_r,P_c]=fPointFind(P_list(P_test),PMap);
%                 [pattern_list,num_patterns]=Folder_Prep(InputUser);
    [RefPat] = ReadEBSDFile(P_list{1},InputUser.PatternFlip);
else
    [ RefPat ] = bReadEBSP(MapInfo.EBSPData,1);
end
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
subplot(1,2,1); imagesc(RefPat); axis image; axis xy; axis tight; colormap('gray'); title('Input pattern');
subplot(1,2,2); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');

EBSD_Info.PatSizeH=size(RefPat_cor,2);
EBSD_Info.PatSizeW=size(RefPat_cor,1);
clear RefPat RefPat_cor

%% Get tiling parameters
if PCA_Setup.Grain_info==1
    % InputUser.Phase_Input={'Nb6Co7'};
    PCA_Setup.cropping=0;
    MapInfo.quantil=PCA_Setup.quantil
    MapInfo.ori_lim=PCA_Setup.ori_lim
    MapInfo.max_pat_num=PCA_Setup.max_pat_num;
    [ ~,~,~,~,~,RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
    disp(RTM_info.cif_file)
    cs=loadCIF(RTM_info.cif_file); 
    MapInfo=GrainSetup(MapInfo,cs);
else
    MapInfo=TilingSetup(MapInfo,PCA_Setup.crop_factor);
    % MapInfo = TilingSetup( MapInfo, PCA_Setup.crop_factor);
    if PCA_Setup.adjustTiling ==1
        MapInfo = adjustTilingSetup( MapInfo, PCA_Setup);
    end 
end
%% Now loop over the inputted weights

for weightloop=1:length(PCA_Setup.EBSD_weighting)%% Setup PCA, load data, generate tiles.

%Main loop over tiled areas
%Results from full AOI are stored in 'tile' structure.
%Results from RTM are stored in 'RTM' structure.
clear tile

vartol=PCA_Setup.variance_tolerance(1);
weighting=PCA_Setup.EBSD_weighting(weightloop);
disp(['weighting  ',num2str(weightloop),' of ',num2str(length(PCA_Setup.EBSD_weighting))]);
cd(InputUser.SavingFolder)

%pre-allocate empty arrays for spectra, patterns, etc for concatenation.
tile.pat=[];
tile.spec=[];
tile.av_spectra=[];
tile.weight=PCA_Setup.EBSD_weighting(weightloop);

if printing==1
    FolderSetup_PCA;
end
close all

%do a quick adjustment of cropfactor if only want to run one tile
if isnan(PCA_Setup.run_one_tile)~=1
    savetiling=PCA_Setup.crop_factor;
    PCA_Setup.crop_factor=1;
end
save(['Results_intermediate'],'MapInfo')
%% Analyse the data
if PCA_Setup.Grain_info==1
    GrainedPCA;
else
    TiledPCA;
end

% GrainedPCA;
if ~PCA_Setup.Grain_info==1
TiledPCA_Reassembly;
end
%Run the RTM if required
if PCA_Setup.RTM==1
   RTM_PCA;
   %RTM.Output=RTMOutput;
   clear RTMOutput;
end

%% Save results
%Save the results. Will overwrite as it goes through the
%weightloops/tolerances...
if printing==1
    cd(InputUser.ResultsDir)
    if PCA_Setup.RTM==1
        if PCA_Setup.PCA_EDX==1
            save(['Results'],'PCA_Setup','tile','RTM','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        else
            save(['Results'],'PCA_Setup','tile','RTM','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        end
    
    else 
        if PCA_Setup.PCA_EDX==1
            save(['Results'],'PCA_Setup','tile','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        else
            save(['Results'],'PCA_Setup','tile','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        end
    end   
end

%% Spectra saving for RC-spectra
if PCA_Setup.PCA_EDX==1
    if printing==1
        cd(InputUser.ResultsDir)
        mkdir('RCSpectra')
        cd('RCSpectra')

        %Grabs a .spx file from somewhere else, then changes channel numbers to
        %those of RC-spectra. NB operates on tile.spec if needs to change
        SPX_replace_PCA(fullfile(InputUser.Astro_loc,'modules','pca','Analysis','Bruker_Spectra_Hack'),'spec.spx',tile.spec,pwd)
        cd(InputUser.ResultsDir)

        mkdir('AverageSpectra')
        cd('AverageSpectra')
        %Grabs a .spx file from somewhere else, then changes channel numbers to
        %those of RC-spectra. NB operates on tile.spec if needs to change
        SPX_replace_PCA(fullfile(InputUser.Astro_loc,'modules','pca','Analysis','Bruker_Spectra_Hack'),'spec.spx',tile.av_spectra,pwd)
        cd(InputUser.ResultsDir)

    end

    clear ret_comps sel_comps weighting ComponentMap ElementMaps InputPats tiling_no total_ret_comps n
    %clear tile EDSWindowing
end

clear weighting vartol
end %end of weightloop
clear weightloop

cd(InputUser.SavingFolder)
close all

pTime(['completed ',num2str(fileloop), ' of ' num2str(length(files))],t1);
end %end of fileloop
clear fileloop ans
