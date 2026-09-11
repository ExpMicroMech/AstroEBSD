%% L. Berners (Beginning of 2025), This deck uses a combination of the fourier based 
% refined template matching and the real space based optimisation routine.
% You need to specify a phase and a rough estimate of the pattern centre,
% and the deck will try to refine the PC and OR at the same time.
%
%
%
% Clean up workspace
clear ; 
close all;
%%
InputUser.MTEX_loc='/work/mz071159/mtex';
% InputUser.MTEX_loc='C:\Users\berners\MTEX_dev\mtex';
try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end

InputUser.Phase_Folder="/work/mz071159/AstroEBSD_v2/phases";
InputUser.Astro_loc='/work/mz071159/AstroEBSD_v2';
try astro_loadcheck
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end



InputUser.Phase_Folder="/work/mz071159/AstroEBSD_v2/phases";
InputUser.Astro_loc='/work/mz071159/AstroEBSD_v2';
disp(InputUser)
disp(InputUser.Phase_Folder)
% warning('off','MATLAB:MKDIR:DirectoryExists')
Phase_Folder=InputUser.Phase_Folder

% color coding and FHKL plotting
cmap_div =  returnCmapDiv();





%% Set up file names
% InputUser.Phase_Input = {'Nb518Ni482_joubert_19_2025_08_01__11_26_28_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000'}
% InputUser.Phase_Input = {'Nb51_Ni49_joubert_19.5kV_dmin0.050_500npx_4.0c1_8.0c2'}
subFolder ="__PHASESUBFOLDER__"
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead)
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end
% InputUser.Phase_Input = InputUser.Phase_Input(1)
plotOnly =__PLOTONLY__;

InputUser.refineLatticeParams=__LPRMODE__;

% InputUser.Phase_Input  = {'NbNi_dt0.450_c26.7_a4.9_19.3kV_dmin0.020','NbNi_dt0.250_c26.7_a4.9_19.3kV_dmin0.020'}
% InputUser.Phase_Input = {'YNi2'}
%%
InputUser.OSC_file = '__EBSDMAP__.ang';
% InputUser.OSC_file = 'MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3.ang'
InputUser.OSC_folder =  "__EBSDFOLDER__";
% InputUser.OSC_folder = 'J:\Users\Berners_Lukas\clara\YNi2MaDinah\PCA\MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3\_EBSD\vt0,1'
% InputUser.OSC_folder = 'E:\sciebo\IMM\BandProfiles\NbNi49_ARea6'
% InputUser.PatternFolder = 'E:\sciebo\IMM\BandProfiles\NbNi49_ARea6'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 6__EBSD+EDS Map 2\_EBSD\vt0,1\Tile1\RotPCA_n44_w1,00e+00'
% InputUser.PatternFolder = "J:\Users\Berners_Lukas\clara\YNi2MaDinah\PCA\MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3\_EBSD\vt0,1"
InputUser.PatternFolder = "__PATTERNFOLDER__";
% InputUser.PatternList = ["19kV_Spot_notFiltered_BGDivGauss.tif","22kV_Spot_Filtered_THL224_BGDivGauss.tif"]%,"_Pattern_2.png"]
InputUser.ExportPath = InputUser.OSC_folder 
% InputUser.ExportPath = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 6__EBSD+EDS Map 2\PatternMatching'

InputUser.AnalysisPCHDF5File = fullfile(InputUser.ExportPath,sprintf('__PCCALIBFILE__.hdf5'));


% InputUser.fitPC = 1;
InputUser.ExportPath = fullfile(InputUser.OSC_folder,subFolder);
mkdir(InputUser.ExportPath);




% PCA_Mode = 'VMoutput' %switch between loading complete tiles or only the _pca_.mat files
% PCA_Mode = "VMoutput"
% InputUser.inputEulerMode = "commercial"
% InputUser.fitPC = 1;
% %%InputUser.ExportPath = InputUser.PatternFolder 
% % InputUser.ExportPath = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 6__EBSD+EDS Map 2\PatternMatching'
% InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('SinglePatRefinementWithLPR.hdf5'))
% InputUser.AnalysisPCHDF5File = fullfile(InputUser.ExportPath,sprintf('PatternPlots_Temp.hdf5'))
% InputUser.OldPCAFile = '_PCA_.mat'
% InputUser.fitPC = 1;



% PCA_Mode = 'VMoutput' %switch between loading complete tiles or only the _pca_.mat files
PCA_Mode = "__PCAMODE__"
if strcmp(PCA_Mode,'VMoutput')
    InputUser.OldPCAFile = '_PCA_.mat';
else
    InputUser.OldPCAFile = 'Results.mat';
end
% InputUser.fitPC = 1;
InputUser.PCA_Mode=PCA_Mode

InputUser.inputEulerMode = "__EULERMODE__";
InputUser.fitPC = __PCFITTED__;
selectPCAuto = __SELPCAUTO__;
refineFFT = __REFINEFFT__;

if InputUser.refineLatticeParams==1
    RefineLPRString = "RefineLPR";
else
    RefineLPRString = "";
end

h5exportName = sprintf('%s_SinglePatRefinement_eangs%s_pcFitted%i_PCAMode_%s_%s___PHASESUBFOLDER__.hdf5',...
          InputUser.OSC_file(1:end-4),...
          InputUser.inputEulerMode,InputUser.fitPC,PCA_Mode,RefineLPRString ...
    )
InputUser.h5name = fullfile(InputUser.ExportPath,h5exportName)

%%
for i = 1:numel(InputUser.Phase_Input)
     InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(1),InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % load ebsd map
    rtmConfig = defaultRTMConfig();
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end

%% load ebsd map

[ MapInfo.MapData,MicroscopeData,~,MapInfo.EBSPData,InputUser ] = loadEBSDEDAXClarity(InputUser);

InputUser.calibration = loadCalibrationEDAX( fullfile( InputUser.OSC_folder,InputUser.OSC_file) );
InputUser.calibration.DetPXSize = 55 % (mm, for EDAX Clarity. For other detector please change.
InputUser.calibration.DetSize = InputUser.calibration.DetPXSize*514 % 
InputUser.calibration.DetNumPixels = 514;
MicroscopeData.TotalTilt=((90-InputUser.calibration.SampleTilt)+InputUser.calibration.CameraElevationAngle)*pi/180
InputUser.PC = [InputUser.calibration.pcx,InputUser.calibration.pcy,InputUser.calibration.pcz];
[pcx, pcy, pcz] = patternCentreGeometric(InputUser.calibration,MapInfo.MapData);
MapInfo.MapData.pcx = pcx;
MapInfo.MapData.pcy = pcy;
MapInfo.MapData.pcz = pcz;
InputUser.MicroscopeData = MicroscopeData;
InputUser.MapInfo = MapInfo;
% ------------------------------ PC Refinement files ----------------------------------
% NOTE: This needs to be done after loadEBSDEDAXClarity
InputUserPCRefine = InputUser
PCRefineShift = __PCREFINESHIFT__
InputUserPCRefine.OSC_folder =  "__PCREFINEFOLDER__"
InputUserPCRefine.AnalysisPCHDF5File = fullfile(InputUserPCRefine.OSC_folder,"__PCREFINEANALYSIS__.hdf5")
InputUserPCRefine.OSC_file = "__PCREFINEOSC__.ang"




%% load tile 1
% Hard coded for now.
if strcmp(PCA_Mode, 'VMoutput')
    tileNum = 1;
    tile.rowstart = 1;
    tile.rowfin = 83;
    tile.colstart = 1;
    tile.colfin = 83;
end
% [XSearchGrid ,YSearchGrid ]= meshgrid(tile.rowstart(tileNum):1:tile.colfin(tileNum),tile.rowstart(tileNum):1:tile.rowfin(tileNum))
% [YSearchGrid ,XSearchGrid ]= meshgrid(tile.rowstart(tileNum):1:tile.rowfin(tileNum),tile.colstart(tileNum):1:tile.colfin(tileNum));
% [pnum,~]=find(MapInfo.MapData.XBeam==XSearchGrid(:)' & MapInfo.MapData.YBeam==YSearchGrid(:)');



%% check that we found the rigth stuff
% % % % plotx2east
% % % % plotzIntoPlane
% % % % figure
% % % % plot(MapInfo.EBSPData.ebsd_map(pnum),MapInfo.EBSPData.ebsd_map(pnum).prop.iq)
% % % % figure
% % % % plot(MapInfo.EBSPData.ebsd_map,MapInfo.EBSPData.ebsd_map.prop.iq)
% % % % % seems to work !!

%% load the PCA on one tile only


if strcmp(PCA_Mode, "Tiles")
    load(fullfile(InputUser.PatternFolder,InputUser.OldPCAFile),'tile','MicroscopeData');
elseif strcmp(PCA_Mode, "VMOutput")
    load(fullfile(InputUser.PatternFolder,InputUser.OldPCAFile),'VMOutput','MicroscopeData');
end
%% Try to find the best matched plane from old analysis
% DataSetName = 'Patterns/Pat_Size_514/MaxIter_40/Nb518Ni482_joubert_19_2025_08_01__11_26_29_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.600_npx-1000' 
%TODO: Needs to be automated this is now selected from the plot we did using 
% the fitPlaneToPatternCentres.m deck
try
    if (InputUser.fitPC==1)
        PatternData=readPatternDataFromH5(InputUser.AnalysisPCHDF5File,'skipPatterns')
        %%
        if selectPCAuto==-1
            warning('Setting Static PC later')
        elseif selectPCAuto==0
            PatternData = findBestMatchingPhase_PCRefinement(PatternData,514) %%
            maxStdDev = 1;% avodi magic numbers so as variable
            [PCX_mod,PCY_mod,PCZ_mod] = getPCPlaneParams(InputUser,PatternData,maxStdDev);
        else
            PatternData = PatternData(selectPCAuto);
            PatternData.MapData = MapInfo.MapData;
            maxStdDev = 1;% avodi magic numbers so as variable
            [ProjectiveTransform] = getPCPlaneParamsProjective(InputUser,PatternData,maxStdDev);
            
        end
        
        % needs improvement
        % also hard codes for now
        
        
        
    end
end
try
[MapInfoPCRefine.MapData,~,~,~] =   loadEBSDEDAXClarity(InputUserPCRefine);

if (InputUser.fitPC==1)
    % try
    PatternData=readPatternDataFromH5(InputUserPCRefine.AnalysisPCHDF5File,'skipPatterns')
    %

    if selectPCAuto==-1
        warning('Setting Static PC later')
    elseif selectPCAuto==0
        PatternData = findBestMatchingPhase_PCRefinement(PatternData,514) %%
        PatternData.MapData =MapInfoPCRefine.MapData;
        maxStdDev = 1;% avodi magic numbers so as variable
        [PCX_mod,PCY_mod,PCZ_mod] = getPCPlaneParams(InputUserPCRefine,PatternData,maxStdDev);
    else
        PatternData = PatternData(selectPCAuto)
        PatternData.MapData = MapInfoPCRefine.MapData;
        maxStdDev = 1;% avodi magic numbers so as variable
        [PCX_mod,PCY_mod,PCZ_mod] = getPCPlaneParams(InputUserPCRefine,PatternData,maxStdDev);
        
    end
    % if isfield(InputUserPCRefine,'OSCPCRefineFile')
    %
    if PCRefineShift
        InputUserPCRefine.calibration = loadCalibrationEDAX( fullfile( InputUserPCRefine.OSC_folder,InputUserPCRefine.OSC_file) );
        DeltaWD = InputUser.calibration.WD-InputUserPCRefine.calibration.WD
        % InputUser.calibration.DetSize 
        DeltaPCY = cos(InputUser.calibration.CameraElevationAngle*degree)* (DeltaWD/InputUser.calibration.DetSize )
        DeltaPCZ = sin(InputUser.calibration.CameraElevationAngle*degree)* (DeltaWD/InputUser.calibration.DetSize )
        PCZ_mod(3) = PCZ_mod(3)+DeltaPCZ;
        PCY_mod(3) = PCY_mod(3)+DeltaPCY;
        
        % end
    end

    
    
    
end
end
% % needs improvement
% % also hard codes for now



% end
% end
    % getPCPlaneParams(InputUser,PatternData,maxStdDev);
% end
    % assumes that the map is centered for now.
    % [PCX_mod,PCY_mod,PCZ_mod] = getPCPlaneParamsGeometric(InputUser,MapInfo,MicroscopeData);
% end


%% Get info for structure Parameters
FHKLoptions.Voltage = 20E3; %TODO: read from phase file.
structConfig = defaultStructureFactorConfig(FHKLoptions);     


%% Set up template matching stuff
rtmOptions.refine=1;
rtmOptions.clarity=1; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions)
rtmConfigRefine.generalOptions.AlignEdax =1 
% rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start
rtmConfigRefine.gaussBlur=0

% if we want to use pseudosymmetry. Specify as many matrices as you want.
rtmConfigRefine.options.TestPseudoSymmetry = 1;
rtmConfigRefine.options.rotPseudo = cat(3,__ROTATIONMATRIXPSEUDO__) % for tetragonal Laves phases


%%rtmConfigRefine.RTM.Rz(60*degree),...
          %%                          rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(30*degree),...
            %%    rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(60*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(120*degree)) % for orthorombic laves phase Fe Ta Al
rtmConfigRefine.options.numPseudo=size(rtmConfigRefine.options.rotPseudo,3);
%% otions for refinement % will influence the speed of the refinement by a lot.
optimisationOptions.MaxIter = __MAXITER__ ;
optimisationOptions.fminSearch = 0;
optimisationOptions.StepTolerance = 1e-12;
optimisationOptions.OptimalityTolerance = 1e-12;
optimisationOptions.BGVarVal = 00;
optimisationOptions.Diagnostics = 'off';
optimisationOptions.Display = 'final';
optimisationOptions.gradients = __GRADIENTS__;
% optimisationOptions.FunValCheck = 'off'; % should work well with NCC brakes, if the NCC is really bad will save comp. time
% optimisationOptions.minFunVal = 0.001 ; %% try to break if the funciton is too bad.
%% reprojections
InputUser.reprojectionSize = 2000;
%% GetrtmConfigRefine Interpolants
rtmConfigRefine.Settings_PCin.start =InputUser.PC %just for dummy
MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfigRefine    )
% try
%     writeStructToH5(rtmConfig, InputUser.h5name,"/rtmConfigRefine")
% catch
% warning("Could not write rtmConfig to h5")
% end
% rtmConfig.Settings_PCin.sart = [0.529145, 1-0.679621, 0.693581]

%% Start h5 writing

if ~plotOnly==1
try
    makeNewH5(InputUser.h5name)
catch
 warning("File could not be created, mak sure, you have setup git properly on your machine")
end


try
    writeStructToH5(InputUser, InputUser.h5name,"/UserInfo")
catch
 warning("Could not write Input Info to h5")
end

try
    writeStructToH5(MicroscopeData, InputUser.h5name,"/MicroscopeData")
catch
 warning("Could not write MicroscopeData Info to h5")
end
   
try
    writeStructToH5(MapInfo, InputUser.h5name,"/MapInfo")
catch
 warning("Could not write MapInfo to h5")
end

try
    writeStructToH5(rtmConfig, InputUser.h5name,"/rtmConfig")
catch
 warning("Could not write rtmConfig to h5")
end
%%
end

if strcmp(PCA_Mode, "Tiles")
    match = getTileNames(tile)
    t_size=size(match,2);

    % tile.MapInfo = 
else
    t_size =1;
end
%%
for tileNum =1:t_size
    tile.MapInfo = MapInfo;
    tile.MicroscopeData = MicroscopeData;
    tile.tileNum = tileNum;
    % tile.numPatterns = 10;
    phaseId = 1;
try
     
    rtmConfigRefine.RTM.ProjectiveTransform = ProjectiveTransform;
    catch
    try
    rtmConfigRefine.RTM.PCX_mod = PCX_mod;
    rtmConfigRefine.RTM.PCY_mod = PCY_mod;
    rtmConfigRefine.RTM.PCZ_mod = PCZ_mod;
    catch
        
        warning ("PCX mod not defined will use same PC for every coordinate")
    end
end
%%
% tileNum = 1;
% numPatterns = 10;
pnum=1
if ~plotOnly==1
% if strcmp(InputUser.inputEulerMode, "commercial")
if strcmp(PCA_Mode, "Tiles")
    tileName = match{tileNum}
    getOrisPCATile(tile,tileName,rtmConfigRefine,InputUser);
elseif strcmp(PCA_Mode, "Full_Map")
    getOrisFullMap(rtmConfigRefine,InputUser,MapInfo,MicroscopeData);
else
    getOrisPCAVMOut(VMOutput,tile,rtmConfigRefine,InputUser);
end


%% TEMPORARILY adding skipping of the rest of loop
% continue 
%%
if ~ strcmp(InputUser.inputEulerMode, "commercial")
    % if strcmp(PCA_Mode, "Tiles")
    %     tileName = match{tileNum}
    %     preparePCARefineTile(tile,tileName,rtmConfigRefine,InputUser);
    % else
    %     preparePCARefineVMOutput(VMOutput,tile,rtmConfigRefine,InputUser);
    % end
    rtmConfigRefineCoarse =rtmConfigRefine
    rtmConfigRefine.RTM.screensize=double(64);
    rtmConfigRefine.RTM.LPTsize = 2* double(rtmConfigRefine.RTM.screensize);
    % RTM.LPTsize = 64;
    rtmConfigRefine.Settings_Cor.size=128;
    rtmConfigRefine.Settings_Cor.screensize=rtmConfigRefine.Settings_Cor.size;
    pat=h5read(InputUser.h5name,"/RTM/RTMPCA/1/PCAPatterns");
    pat = pat(:,:,1);
    [patternIn,InputUser.PatternInfo] = EBSP_BGCor(pat,rtmConfigRefine.Settings_Cor);
    InputUser.EBSD_geometry = EBSP_Gnom( InputUser.PatternInfo,InputUser.PC);
    
    % [PhaseParameters,FFTParameters,EBSD_geometry] = prepareRefinementMultiPhase(InputUser,rtmConfigRefine)
    % templateMatch()
    % realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct)
    tic
    libraryMatchPatternsCoarse(InputUser,rtmConfig,MultiPhaseStruct,tileNum)
    toc
    rtmConfigRefine.Settings_Cor.size=256;
    rtmConfigRefine.Settings_Cor.screensize=rtmConfigRefine.Settings_Cor.size;
    rtmConfigRefine.RTM.Settings_Cor = rtmConfigRefine.Settings_Cor; % this way, refine5 will also internally correct the pattern.
    InputUser.MicroscopeData =MicroscopeData;
    PatternRefineFFT(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum);
    % realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct)

end
%% Patterns = VMOutput.PCA_VM_Pat_n;
if refineFFT==1
% tileNum = 1;
rtmConfigRefine.Settings_Cor.size=256;
rtmConfigRefine.Settings_Cor.screensize=rtmConfigRefine.Settings_Cor.size;
rtmConfigRefine.RTM.Settings_Cor = rtmConfigRefine.Settings_Cor; % this way, refine5 will also internally correct the pattern.
InputUser.MicroscopeData =MicroscopeData;
PatternRefineFFT(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum) ;
% isfield(RTM_setup,'SettingsCor')%%
end
% numPattePatternRefineFFTns = 2
% tileNum = 1

%% set up plotting parameters
rtmConfigRefine.Settings_Cor.size=514;
rtmConfigRefine.Settings_Cor.screensize=rtmConfigRefine.Settings_Cor.size;
rtmConfigRefine.Settings_Cor.cross=2;
rtmConfigRefine.Settings_Cor.gaussian_range = 1;
rtmConfigRefine.options.initPatSize = [514 514];
rtmConfigRefine.Settings_Cor.gaussBlur = 1;
rtmConfigRefine.Settings_Cor.quadrants=1;
rtmConfigRefine.Settings_Cor.max_var_lattice = 0.05;
% rtmConfigRefine.gaussBlur=0

%% Optimise in real space if wanted
InputUser.InputDataset ="RTMRefined"


% if InputUser.fi
if InputUser.fitPC == 0;
    optimisationOptions.refinePC = 1;
    realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions);
else
    optimisationOptions.refinePC = 0;
    realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions);
end

try
    writeStructToH5(rtmConfigRefine,InputUser.h5name,"/RTM/RTMRefined")
catch
 warning("Could not write Input Info to h5")
end

try
    writeStructToH5(optimisationOptions,InputUser.h5name,"/RTM/RTMRefined")
catch
 warning("Could not write Input Info to h5")
end





%% Refine Lattice Parameters if wanted

if ~isfield(InputUser,'refineLatticeParams')
    InputUser.refineLatticeParams=0
end
if ~(InputUser.refineLatticeParams==0)
    % InputUser.Inp
    realSpaceOptimisationLattice(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions);
end
end
if __REFINEONLY__
    continue
end
%% Plot refinement Results
plotRefinementResults(InputUser,rtmConfigRefine,structConfig,MultiPhaseStruct,tileNum)



InputUser.h5nameRepro = sprintf("%s_repro",InputUser.h5name);
%% Plot reprojections
% Inpu
% tileNum =1
plotReprojectionsFromH5(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum)
%%
%% WARNING: DO NOT RUN. All profiles will be plotted NEED debugging
try
InputUser.numHKLPLot=10;
SampleBandProfilesFromReprojections(InputUser,structConfig)
% 
catch
    warning('could not sample bandprofiles')

end

end




%% END of main deck, below we have experimental stuff
%############################################################################
%##                                                                        ##             
%##                                                                        ##
%##                                                                        ##
%##             END of main deck, below we have experimental stuff         ##
%##                                                                        ##
%##                                                                        ##
%##                                                                        ##
%############################################################################
%%% end
%% Test Stuff with the reprojected patterns To make comparison better work!
for tileNum=1:t_size
for phaseNumber = 1
    InputUser.imDir = fullfile(InputUser.ExportPath,sprintf("Tile_%03i",tileNum),InputUser.Phase_Input{phaseNumber},'ProfilesAlongBand');
        mkdir(InputUser.imDir);
    phase_name = InputUser.Phase_Input{phaseNumber}
    expPatternsRepro = h5read(InputUser.h5nameRepro,sprintf('/RTM/%s/RTMPCA/1/ExpPatternRepro',phase_name));
    simPatternsRepro = h5read(InputUser.h5nameRepro,sprintf('/RTM/%s/RTMPCA/1/SimPatternRepro',phase_name));
    %%
    expPattern1 = expPatternsRepro(:,:,2); %the first one is likely faulty
    
    simPattern1 = simPatternsRepro(:,:,2);
    %%
    figure;
    imagesc(expPattern1);
    colormap gray;
    axis('equal','off','xy');
    figure;
    imagesc(simPattern1);
    colormap gray;
    axis('equal','off','xy');
    %% from getBandProfilesSTereographic
    
    resolution=2/(size(expPattern1,1)-1);
    [gr_x,gr_y]=ndgrid(-1:resolution:1,-1:resolution:1);
    
    expPattern1 = expPattern1-mean(expPattern1(:));
    expPattern1 = expPattern1/std(expPattern1(:));
    simPattern1 = simPattern1-mean(simPattern1(:));
    simPattern1 = simPattern1/std(simPattern1(:));
    
    stereo_patternExp=griddedInterpolant(gr_x,gr_y,expPattern1,'cubic','none');
    stereo_patternSim=griddedInterpolant(gr_x,gr_y,simPattern1,'cubic','none');
    proj_ster=eangleProjection;
    
            
    %% setup structure factor bits for plotting the bands
    %% Get info for structure Parameters
    FHKLoptions.Voltage = 20E3;
    structConfig = defaultStructureFactorConfig(FHKLoptions);     
    phaseNumber = 1
    structConfig.Phase_Folder = InputUser.Phase_Folder;
    [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
    [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
    if ~isfield(InputUser,'numHKLPlot')
        InputUser.numHKLPLot=25;
    end
    [m_plot,I] =cleanFHKLDuplicates(m_plot,I,InputUser);
    structConfig.IPlot = I;
    % i=1;
    
    
    %% Set up angular resolution
    
    
    theta=(-180:0.1:180)*degree;
       f=figure;
       f.Position = [0 0 1500 1000]
       h=tiledlayout(1,1,"TileSpacing","none");
       nexttile(1)
    imagesc(stereo_patternExp.GridVectors{1},stereo_patternExp.GridVectors{2},stereo_patternExp.Values)
    colormap gray
    axis xy
    axis equal
    axis off
    clim([-4 4])
    xlim([-1 1])
    ylim([-1 1])
    hold on
    for millerNumber =1:numel(m_plot)

        m = m_plot(millerNumber);
        omega = [pi/2-asin(structConfig.lambda_p/(2* m.dspacing)) pi/2+asin(structConfig.lambda_p/(2* m.dspacing))]';
        samplingPlot=circles(vector3d(m.x,m.y,m.z),theta,omega);
        samplingPlot = samplingPlot(samplingPlot.z<0);
        [x_sample_discretePlot,y_sample_discretePlot]=proj_ster.project(samplingPlot);
        scatter(x_sample_discretePlot,y_sample_discretePlot,3,'filled','DisplayName',char(m));
    end
    legend('Location','eastoutside');
    % legend()
    exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparisonBandsLabelled_%s.png',phase_name)));
    %%
    theta=(0:0.1:180)*degree;
    BandWidth = 5*degree;
    omegaResolution = 0.05;
    omega=[pi/2-BandWidth:omegaResolution*degree:pi/2+BandWidth]; 
    omega_size=size(omega,2);
    % if you want to see an image uncomment these lines

    %%
    for millerNumber =1:numel(m_plot)
    m=m_plot(millerNumber);
    sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
    sampling = reshape(sampling,[],omega_size);
    thetaPlot = theta(any(sampling.z<0,2));
    sampling = sampling(any(sampling.z<0,2),:);

    [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
    valsExp=stereo_patternExp(y_sample_discrete,x_sample_discrete);
    valsSim=stereo_patternSim(y_sample_discrete,x_sample_discrete);
    valsExp=reshape(valsExp,[],omega_size);
    valsSim=reshape(valsSim,[],omega_size);
    valsExp = valsExp-mean(valsExp(:),'omitnan');
    valsExp = valsExp/std(valsExp(:),'omitnan');
    valsSim = valsSim-mean(valsSim(:),'omitnan');
    valsSim = valsSim/std(valsSim(:),'omitnan');
    
    %%
    %%
  
 
    %%
    % figure
    % plot(theta/degree,valsExp+omega/degree,'Color',[0.5 0.5 0.5])
    %%
    % numBins = floor(BandWidth/0.5/degree)
    binSize = 0.5;
    omegaResolution;
    binRes = floor(binSize/omegaResolution);
    numBins = floor(size(valsExp,2)/binRes);
    %
    offset = binSize;
    averageExp = zeros([size(valsExp,1),numBins]);
    averageSim = zeros([size(valsExp,1),numBins]);
    centerBin = ceil(size(valsExp,2)/2);
    omegaResolution;
    offsets = ones(1,numBins)*90;
    for binNumb = -ceil(numBins/2)+1:floor(numBins/2)-1
    
        averageExp(:,binNumb+ceil(numBins/2)+1) = mean(valsExp(:,centerBin-ceil(binRes*(-binNumb+0.5)):centerBin-floor(binRes*(-binNumb-0.5))),2);
        averageSim(:,binNumb+ceil(numBins/2)+1) = mean(valsSim(:,centerBin-ceil(binRes*(-binNumb+0.5)):centerBin-floor(binRes*(-binNumb-0.5))),2);
        offsets(binNumb+ceil(numBins/2)+1) = 90+binNumb*offset;
    end

 % Do some normalisation
 %theta
    averageSim = averageSim-mean(averageSim(:),'omitnan'); % correct for std
    averageSim = averageSim/std(averageSim(:),'omitnan');
 %
    averageExp = averageExp-mean(averageExp(:),'omitnan'); % correct for std
    averageExp = averageExp/std(averageExp(:),'omitnan');

    averageSim = averageSim/3; % correct for std
    averageSim = averageSim*binSize; % correct for bin Size
    
    averageExp = averageExp/3; % correct for std
    averageExp = averageExp*binSize; % correct for bin Size
    %%
    if max(thetaPlot/degree)>=179.99
        xlims = [90 180]
    else
        xlims = [0 90]
    end
    figure
    
    
    f=figure;
    f.Position =[0 0 1900 1500];
    h = tiledlayout(6,1,"TileSpacing","none")
    nexttile(1)
    imagesc(thetaPlot/degree,omega/degree,valsExp');
    axis equal;
    ylim([90-BandWidth/degree, 90+BandWidth/degree])
    xlim(xlims) % only half of it
    colormap 'gray'
    clim([-4 4])
    
    
    nexttile(2);
    imagesc(thetaPlot/degree,omega/degree,valsSim');
    axis equal;
    ylim([90-BandWidth/degree, 90+BandWidth/degree])
    xlim(xlims) % only half of it
    colormap 'gray'
    clim([-4 4])
    
    ax3 = nexttile(3);
    imagesc(thetaPlot/degree,omega/degree,valsSim'-valsExp');
    axis equal;
    ylim([90-BandWidth/degree, 90+BandWidth/degree])
    xlim(xlims) % only half of it
    colormap(ax3,returnCmapDiv)
    clim([-4 4])
    
    
    nexttile([3 1])
   
    plot(thetaPlot/degree,averageExp+offsets,'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
    xlim(xlims)
    xlabel("azimuth angle")
    ylabel("band opening angle (°) ")
    hold on
    plot(thetaPlot/degree,averageSim+offsets,'Color',[0 1 1],'DisplayName','Simulation')
    % legend()
    exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparison_%s_%s.png',phase_name,char(m_plot(millerNumber)))))
    %%
    f=figure;
    f.Position = [100 100 1500 1000];
    for i = 1:size(averageExp,2)
        if i == 1
        plot(thetaPlot/degree,averageExp(:,i)+offsets(i),'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
            hold on
        plot(thetaPlot/degree,averageSim(:,i)+offsets(i),'Color',[0 1 1],'DisplayName','Simulation')
        else
             plot(thetaPlot/degree,averageExp(:,i)+offsets(i),'Color',[0.5 0.5 0.5],'HandleVisibility','off')
            hold on
        plot(thetaPlot/degree,averageSim(:,i)+offsets(i),'Color',[0 1 1],'HandleVisibility','off')
        end
    end
        % en
    xlim(xlims);
    ylim([90-BandWidth/degree, 90+BandWidth/degree]);

    legend()
    xlabel("azimuth angle")
    ylabel("band opening angle")
    exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparisonOnlyBands_%s_%s.png',phase_name,char(m_plot(millerNumber)))))
    %% try to get the stuff 3d
    viewAngle = 20;
    viewAngle2 =90;
    zlimVals = 5;
    thetaLimMinus = xlims(1);
    thetaLimPlus = xlims(2);
    cmap = returnCmapDiv;
    cmap ='gray';
    [X,Y] = meshgrid(thetaPlot/degree,omega/degree);
    f=figure
    f.Position =[0 0 500 800]
    h = tiledlayout(1,3,"TileSpacing","none")
    nexttile(1)
    surf(X,Y,valsSim','EdgeColor','none')
    title('experimental')
    axis equal
    xlim([thetaLimMinus,thetaLimPlus])
    zlim([-zlimVals,zlimVals])
    clim([-zlimVals, zlimVals])
    view(viewAngle2,viewAngle)
    colormap(cmap)
    nexttile(2)
    surf(X,Y,valsExp','EdgeColor','none')
    axis equal
    xlim([thetaLimMinus,thetaLimPlus])
    zlim([-zlimVals,zlimVals])
    clim([-zlimVals, zlimVals])
    view(viewAngle2,viewAngle)
    colormap(cmap)
    title('simulated')
    
    nexttile(3)
    surf(X,Y,valsExp'-valsSim','EdgeColor','none')
    axis equal
    xlim([thetaLimMinus,thetaLimPlus])
    zlim([-zlimVals,zlimVals])
    clim([-zlimVals, zlimVals])
    view(viewAngle2,viewAngle)
    colormap(cmap)
    title('Difference')
    exportgraphics(f, fullfile(InputUser.imDir,...
                sprintf('BAndComparison3d_%s_%s.png', phase_name, char(m_plot(millerNumber)) )...
                 ))
    
    
    %%
    close all
    
    %% prepare outputs
    [m_plot_sym,indexc,~] = unique(m_plot.symmetrise(),'noSymmetry','stable');
    [~,bandGroups]=ind2sub(size(m_plot.symmetrise()),indexc);
    
    bandValues  = zeros(omega_size,numel(m_plot_sym));
    bandNames = strings(numel(m_plot_sym),1);
    for bNum = 1: numel(m_plot_sym)
        bandNames(bNum) = string(char(m_plot_sym(bNum)));
    end
    bandAngle = omega/degree;
    % bandGroups = c;
    
    %%
    % for m_number = 1:numel(m_plot_sym)
    %     % tic
    %     m = m_plot_sym(m_number);
    % 
    %     sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
    %     [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
    %     vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
    %     vals=reshape(vals,[],omega_size);
    %     sum_vals=sum(vals,1)/length(theta);
    %     bandValues(:,m_number)  = sum_vals;
    %     % toc
    % end

end
end
end
