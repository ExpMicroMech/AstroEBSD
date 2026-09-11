%% 17.10.2025 L. Berners
% In this skript we want to test the influence of potential misalignments
% on the comparison metrics between individual on-screen EBSD patterns.

% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc|rootpahtList|rootrootpath|mapFilterString|structureFilterString|additionalFilterString|rPNum|mapDirs|strucutre_dir|mapPNum|strucutre_dir|structPNum|clean_dir|h5File).
close all;
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
InputUser.Phase_Folder=Phase_Folder;
InputUser.Astro_loc=Astro_loc;
warning('off','MATLAB:MKDIR:DirectoryExists')

% exportFolder = 'C:\Users\berners\sciebo\Paper6_OccupancyAndEffectsInNbNi\CompareMisalignment'
% exportFolder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\58Nb-42Ni_v2\Area7"

%%
InputUser.h5folder = clean_dir;
% InputUser.h5folder = "J:\Users\Berners_Lukas\clara\NbCoUBC\56Nb44Co_verticalCut\Area3\Map1"
exportFolderRoot = fullfile(InputUser.h5folder,'compareNCC');
mkdir(exportFolderRoot);

InputUser.h5refinement = h5File; %$fullfile(InputUser.h5folder, "NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area_2_EBSD+EDS Map 3_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__NbCo3aILGrid_L_cleaned.hdf5");
% InputUser.h5refinement = fullfile(InputUser.h5folder, "NbCoUBC_56Nb44Co_vertical_Area 3_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__cleaned.hdf5")
% InputUser.h5refinement  = strrep(InputUser.h5refinement,'_cleaned','')
h5structure = h5info(InputUser.h5refinement,'/RTM');%,'/Patterns','min')
%%

InfoStruct = readStructFromH5(InputUser.h5refinement,'skipPatterns');


%% Set up file names
InputUser.Phase_Input = {};
countNum = 1;
InputUser.Nb = [];
% InputUser.Co = [];
% InputUser.a_lattice = [];
% InputUser.c_lattice = [];
InputUser.ILSpacing = []; 
for phaseNum = 1:numel(h5structure.Groups)
    phaseName = h5structure.Groups(phaseNum).Name;
    phaseName = phaseName(6:end); % drop the /RTM bit
    if strcmp(phaseName,'RTMPCA') || strcmp(phaseName,'tileInfo') || strcmp(phaseName,"RTMRefined") || strcmp(phaseName,"OriCorrectedRefinement")
        continue
    end
    InputUser.Phase_Input{countNum} = string(phaseName);

    
    % matches = regexp(InputUser.Phase_Input{countNum} , pattern, 'tokens');
    nb_number = [];
    % co_number = [];
    nb_number = regexp(InputUser.Phase_Input{countNum} , '3aNb-([\d\.]+)', 'tokens');
    % co_number = regexp(InputUser.Phase_Input{countNum} , 'Co([\d\.]+)', 'tokens');
    % a_lattice = regexp(InputUser.Phase_Input{countNum} , '_a([\d\.]+)', 'tokens');
    % c_lattice = regexp(InputUser.Phase_Input{countNum} , '_c([\d\.]+)', 'tokens');
    % try
    ILSpacing = regexp(InputUser.Phase_Input{countNum} , '_6c2-([\d\.]+)', 'tokens');
    nb_number = nb_number{1}{1};% if ~isempty(nb_number);
    ILSpacing = ILSpacing{1}{1};
    % co_number = co_number{1}{1};% if ~isempty(co_number);
    % a_lattice = a_lattice{1}{1};
    % c_lattice = c_lattice{1}{1};
    if ~isempty(nb_number)
        InputUser.Nb = [ InputUser.Nb; str2num(nb_number)];
        % InputUser.Co = [ InputUser.Co ; str2num(co_number) ]
        % InputUser.a_lattice = [InputUser.a_lattice; str2num(a_lattice)]
        % InputUser.c_lattice = [InputUser.c_lattice; str2num(c_lattice)]
        InputUser.ILSpacing =  [InputUser.ILSpacing; str2num(ILSpacing)];
    end
    % catch
    %     warning('no value assigned')
    % end
    % if ~isempty(co_number)
    % end
        countNum = countNum+1;


    % InputUser.h5PhaseInfo{countNum}.Name = InputUser.h5PhaseInfo{countNum}.cif(1:end-4)
end
%% setup phase Information
for i = 1:numel(InputUser.Phase_Input)
    InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
    Phases = Phase_BuilderSimplified( InputUser.Phase_Input{i},InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % cs.mineral = Phases.xtalname
    % load ebsd map
    rtmConfig = defaultRTMConfig();
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end

%% Setup Detecotr info
MicroscopeData.SampleTilt = 70;
MicroscopeData.CameraTilt = 1;
MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180;
% MicroscopeData.TotalTilt=((90-InfoStruct.UserInfo.calibration.SampleTilt)+InfoStruct.UserInfo.calibration.CameraElevationAngle)*pi/180
InputUser.Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

% InputUser.EdaxConversionMatrix = rtmConfigRefine.RTM.EdaxConversionMatrix;
InputUser.MicroscopeData = MicroscopeData;

refinementString ='RTMRefined';
tileNum = 1;
numTiles = numel(InfoStruct.RTM.tileInfo.colfin);
%%
%% Set up Corrections and refinement options
rtmOptions.refine=1;
rtmOptions.clarity=1; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions);
rtmConfigRefine.generalOptions.AlignEdax =1;
rtmConfigRefine.Settings_Cor.cross=-1;
rtmConfigRefine.Settings_Cor.quadrants=1;
rtmConfigRefine.Settings_Cor.gfiltQuad=1;
InputUser.inputEulerMode = 'commercial';
InputUser.InputDataset = "RTMRefinedLattice";
InputUser.EdaxConversionMatrix= rtmConfigRefine.RTM.EdaxConversionMatrix;
%% Correction Settings
SettingsCor = rtmConfigRefine.Settings_Cor;
SettingsCor.max_var_pc=0.05;
SettingsCor.max_var_ori=2*degree;
SettingsCor.blurf = 1.1;
SettingsCor.blur = 0;
SettingCor.gaussBlur=0;
% gaussian blur may help reduce noise
SettingsCor.SatCor = 0; % Do not do this. adds a lot of noise to optimisation landscape.!!
SettingsCor.metric = "NCC";

SettingsCor.screensize=514;
SettingsCor.size = [SettingsCor.screensize SettingsCor.screensize ];
SettingsCor.gaussian_range = 1;
SettingsCor.quadrants=1;
SettingsCor.gfiltQuad = 1;
SettingsCor.hotpixel = 0;
SettingsCor.cross=2;
rtmConfigRefine.options.initPatSize = SettingsCor.size;
%%
FHKLoptions.Voltage = double(Phases.ehistmin)*10^3; %TODO: read from phase file.
structConfig = defaultStructureFactorConfig(FHKLoptions);

rtmOptions.refine=1;
rtmOptions.clarity=1; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions);
rtmConfigRefine.generalOptions.AlignEdax =1;
% rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start


warning('off','all')
MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfigRefine    );
%%
rtmConfigRefine.gaussBlur=0;
rtmConfigRefine.gaussianRange=1;
%%
rtmConfigRefine = defaultRTMConfig();
rtmConfigRefine.Settings_Cor.screensize = 514;
rtmConfigRefine.Settings_Cor.size = 514;
PC= [0.5 0.25 0.7];
rtmConfigRefine.PC = PC;
%% Read data and build synthetic pattern
pNum = 2;
patternNumber = 1;

patOriginal = h5read(InputUser.h5refinement  ,"/RTM/RTMPCA/1/PCAPatterns");
% patOriginal = patOriginal(:,:,)
PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');
gmatrixCommercialDetectorOrig = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/G_Out',InputUser.Phase_Input{pNum},refinementString,tileNum));
% gmatrixCommercialDetectorOrig =pagemtimes( pagemtimes( pagetranspose(gmatrixRefined),InputUser.Detector_tilt) ...
                                    % , InputUser.EdaxConversionMatrix )  ;


geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber));
rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));
gmatrix = gmatrixCommercialDetectorOrig(:,:,patternNumber);

patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
patSim = EBSP_BGCor(patSim,SettingsCor);

EBSD_To_TileSim = removeCross512(patSim);
EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));

p1Name = InputUser.Phases{1};
MultiPhaseStruct(1).screen_int;
pNumCompare = 1;
metricName = 'NCC';
metrics = zeros(numel(rot_origIn),numel(InputUser.Phase_Input));
% EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);  
% Recalculation of orientation and projection center based on average of mutliple phases.
% needs to be put into a function and cleaned up later on
%
%
%
% 
% tile 
tileInfo_map = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/MapReshaped2'));
rowfin = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/rowfin'));
colfin = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/colfin'));
rowstart = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/rowstart'));
colstart = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/colstart'));
%
% compMap =
% tileInfo_map(rowstart(tileNum):rowfin(tileNum),colstart(tileNum):colfin(tileNum));
%
for tileNum =1: numTiles

% compMap = tileInfo_map(1:30,1:30); % until this is debugged.
compMap = tileInfo_map(rowstart(tileNum):rowfin(tileNum),colstart(tileNum):colfin(tileNum));
exportFolder = fullfile(exportFolderRoot, sprintf("Tile_%03i", tileNum));
mkdir(exportFolder);

PC_origIn = h5read(InputUser.h5refinement ,sprintf('/RTM/RTMPCA/%i/PCFitted',tileNum));
numPhases = numel(InputUser.Phases);
minPhaseNum =1; % exclude lit data from joubert.
numPatterns = size(PC_origIn,2);
gmatrixAll = zeros(3,3,numPatterns,numPhases);
% gmatrixFilderedAll = zeros(3,3,numPatterns);
% eangs_outAll = zeros(3,numPatterns);
eangs_outAllCombined = [];
gmatrixFilderedAllCombined = [];
allNumPixels = [];
% eangs_outAll = zeros(3,numPatterns);
nccAll = zeros(numPatterns,numPhases);
csCompare = InputUser.CSList{end};
PC_all = zeros(numPhases,numPatterns,3);

for pNum = 1:numPhases

    gmatrixCommercialDetectorOrig = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/G_Out',InputUser.Phase_Input{pNum},refinementString,tileNum));
    eangs_out = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/eangs_out',InputUser.Phase_Input{pNum},refinementString,tileNum));
    % gmatrixCommercialDetectorOrig =pagemtimes( pagemtimes( pagetranspose(gmatrixRefined),InputUser.Detector_tilt) ...
                                        % , InputUser.EdaxConversionMatrix )  ;
    % PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');
    % PC_all(pNum,:,:) = PC_origIn';
    
    % geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber));
    % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig);
    % gmatrix = gmatrixCommercialDetectorOrig(:,:,patternNumber);
    gmatrixAll(:,:,:,pNum)=gmatrixCommercialDetectorOrig;
    % eangs_outAll(:,:,pNum) = eangs_out;
      nccAll(:,pNum) = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{pNum},refinementString,tileNum));
    
end
%%
patternNumber = 1;

% meanRot = mean(gmatrixAll(:,patternNumber))
% compareMat = gmatrixAll(:,patternNumber)
% % meanRot = repmat(meanRot,size(compareMat))
% angles = [] % somehow it did not work for an array, which it technically should for mtex
% for pNum =minPhaseNum:numPhases
%     angles = [angles angle(meanRot,compareMat(pNum))/degree];
% end
%%
% angles
% figure
% plot(angles)
%%
maxPlotNum =numel(InputUser.Phase_Input);
% for pNum = minPhaseNum:patternNumber
%     if mean(NCCList(1:maxPlotNum,pNum))>0.7 & abs(std(NCCList(1:maxPlotNum,pNum)))<0.05
%     ifpKey = ipfColorKey(InputUser.CSList{pNum});
%     plotColors = ifpKey.orientation2color(oris{1});
%     plot(NCCList(1:maxPlotNum,pNum)-max(NCCList(1:maxPlotNum,pNum)),'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
%     end
% end
%% 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Getting the Best fit dependent on angle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PC_origIn = h5read(InputUser.h5refinement ,sprintf('/RTM/RTMPCA/%i/PCFitted',tileNum));

metricName = "NCC";

% occPlot = reshape(InputUser.Occupancy,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
% dtPlot = reshape(InputUser.dtSpacing,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
ILPlot = (InputUser.ILSpacing);%,numel(unique(InputUser.Nb)),numel(unique(InputUser.c_lattice)))
NbPlot = (InputUser.Nb);%(unique(InputUser.Nb)),numel(unique
[~,sortIndex] = sort(InputUser.ILSpacing);
NbPlot = reshape(NbPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)));
ILPlot = reshape(ILPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)));
% [~,newIndex] = sortrows(occPlot)
figure
hold on;
ColList =colorcube(numPatterns);
AnlgeOffset = 0;
pcOffset = 0;
% for AnlgeOffset = [0,0.05,0.1,0.2,0.5,1,2]
    % for pcOffset = [0,0.002,0.005,0.01,0.02,0.05,]
rotAxis = vector3d(0,1,0);
offset = rotation.byAxisAngle(rotAxis,AnlgeOffset*degree);
% for 
allMetricsCombined = [];
patternNums = [];
for patternNumber = 1:numPatterns
    geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber)+pcOffset);
    allMetrics = [];
    ncc_test = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{1},refinementString,tileNum));
    if ncc_test(patternNumber)<0.7
        continue
    end
    %%
    % rotIn = rotation.byMatrix(gmatrixAll(:,:,patternNumber,:));
    gmatrix =  filterRotation(gmatrixAll(:,:,patternNumber,:)); % (rot matrix rot matrix pNumber, phase Number)
    gmatrix = gmatrix*offset;
    gmatrix = gmatrix.matrix;
    % gmatrixFilderedAll(:,:,patternNumber)=gmatrix
    % We tested this in realSpaceOptimisation.m
    gTemp =  gmatrix*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(InputUser.Detector_tilt);
    gTemp_mtex= rotation.byMatrix(gTemp');
    eangsOut = [gTemp_mtex.phi1 gTemp_mtex.Phi gTemp_mtex.phi2];
    patOriginal = h5read(InputUser.h5refinement ,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
    patOriginal = patOriginal(:,:,patternNumber);
    patOriginal = EBSP_BGCor(patOriginal,SettingsCor);
    % figure
    % imagesc(patOriginal)
    %%
    for pNum = minPhaseNum:numPhases
    
            % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));

        

    
        patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
        patSim = EBSP_BGCor(patSim,SettingsCor);
        %%
        % EBSD_To_TileSim = removeCross512(patSim);
        % % EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));
        allMetrics = [allMetrics getMatchMetric(patOriginal,...
                            patSim,metricName)];

        if pNum ==minPhaseNum
            imwrite(flipud(normaliseForImwrite(patOriginal)),fullfile(exportFolder,sprintf('Pat_exp%.3fdeg_pcOFF%.3f_Flat_%i.png',AnlgeOffset,pcOffset,patternNumber)));
        end
          if patternNumber<3
               imwrite(flipud(normaliseForImwrite(patSim)),fullfile(exportFolder,sprintf('Pat_exp%.3fdeg_pcOFF%.3f_Flat_%i_sim_%s.png',AnlgeOffset,pcOffset,patternNumber,InputUser.Phase_Input{pNum})))

        end
    end
    % disp(mean(allMetrics))
    if mean(allMetrics>0.65)
        %%
        % figure
        NCCListPlot = reshape(allMetrics,numel(unique(NbPlot)),numel(unique(ILPlot)));
        % 
        % surf(occPlot(newIndex,:),dtPlot(newIndex,:),  NCCListPlot(newIndex,:)-max(NCCListPlot(:)),'FaceColor',ColList(patternNumber,:),'DisplayName',sprintf('Component %i',patternNumber))%,'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
        % hold on
        allMetricsCombined= cat(3,allMetricsCombined,NCCListPlot);
        patternNums = [patternNums, patternNumber];
        gmatrixFilderedAllCombined = cat(3,gmatrixFilderedAllCombined,gmatrix);
        eangs_outAllCombined = cat(3,eangs_outAllCombined,eangsOut);
        allNumPixels = [allNumPixels, sum(compMap(:)==patternNumber)];

    end

end
%% plotting pole figure to compare orientation data
figure
ori=orientation.byEuler(eangs_outAllCombined(:,1,:),eangs_outAllCombined(:,2,:),eangs_outAllCombined(:,3,:),InputUser.CSList{1});
plotPDF(ori,Miller({0,0,0,1},{1,0,-1,0},{1,1,-2,0},{1,1,-2,1},InputUser.CSList{1}),'MarkerSize',2,'projection','eangle');
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_PoleFigure.png',AnlgeOffset,pcOffset)));
close all
% gTemp_mtex


%%
no3d=0
figure
hold on
ColList = colorcube(size(allMetricsCombined,3)+5);


DiffValues = [];
for patNum = 1: size((allMetricsCombined),3)
    % NCCListPlot = allMetricsCombined(:,:,patNum)
    NCCPlot =  allMetricsCombined(:,:,patNum);
    NCCPlot = NCCPlot;
    DiffVals = NCCPlot-max(NCCPlot(:));

    DiffValues = cat(3,DiffValues,DiffVals);

end
    % NCCPlot = reshape(NCCPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
try
for patNum = 1: size((allMetricsCombined),3)
    DiffVals = DiffValues(:,:,patNum);
    surf(ILPlot,NbPlot,  DiffVals-min(DiffValues(:)),'FaceColor',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)));
end
catch
    for patNum = 1: size((allMetricsCombined),3)
    DiffVals = DiffValues(:,:,patNum);
    plot(ILPlot,DiffVals-min(DiffValues(:)),'Color',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)));
end
end
NCCMean = mean(DiffValues,3);
% contour()
try
contour(ILPlot,NbPlot,NCCMean,15, ...
    "FaceAlpha",0.5,'HandleVisibility','off')
catch
    no3d = 1
    warning('not plotting contour')
end
% exportTable = table()
legend('Location','southoutside','NumColumns',3)
ylabel('Occupancy, Nb')
xlabel('dt Spacing')
zlabel('NCC-maxNCC')
if no3d == 0
view(75,30)
end
clim([min(DiffValues(:)) 0 ])
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.png',AnlgeOffset,pcOffset)))
if no3d == 0
view(170,30)
end
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_rot.png',AnlgeOffset,pcOffset)))
% view(170,30)
%%

no3d=0
figure
hold on
ColList = colorcube(size(allMetricsCombined,3)+5);


DiffValues = [];
for patNum = 1: size((allMetricsCombined),3)
    % NCCListPlot = allMetricsCombined(:,:,patNum)
    NCCPlot =  allMetricsCombined(:,:,patNum);
    NCCPlot = NCCPlot;
    DiffVals = NCCPlot-max(NCCPlot(:));

    DiffValues = cat(3,DiffValues,DiffVals);

end
    % NCCPlot = reshape(NCCPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
try
for patNum = 1: size((allMetricsCombined),3)
    DiffVals = allMetricsCombined(:,:,patNum);
    surf(ILPlot,NbPlot,  DiffVals-min(DiffValues(:)),'FaceColor',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)));
end
catch
    for patNum = 1: size((allMetricsCombined),3)
    DiffVals = allMetricsCombined(:,:,patNum);
    plot(ILPlot,DiffVals-min(DiffValues(:)),'Color',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)));
end
end
NCCMean = mean(DiffValues,3);
% contour()
try
contour(ILPlot,NbPlot,NCCMean,15, ...
    "FaceAlpha",0.5,'HandleVisibility','off')
catch
    no3d = 1
    warning('not plotting contour')
end
% exportTable = table()
legend('Location','southoutside','NumColumns',3)
ylabel('Occupancy, Nb')
xlabel('dt Spacing')
zlabel('NCC-maxNCC')
if no3d == 0
view(75,30)
end
clim([min(DiffValues(:)) 0 ])
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFFNoDiff%.3fdeg_pcOFF%.3f.png',AnlgeOffset,pcOffset)))
if no3d == 0
view(170,30)
end
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFFNoDiff%.3fdeg_pcOFF%.3f_rot.png',AnlgeOffset,pcOffset)))
% view(170,30)
%%
try
figure
contourf(NbPlot',ILPlot',NCCMean',15,"ShowText",true,"LabelFormat","%0.3f", ...
    "FaceAlpha",0.75)
ylabel('Interlayerspacing')
xlabel('3a fraction Nb')
hold on
scatter(NbPlot,ILPlot,'black','filled')
colorbar
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_Flat.png',AnlgeOffset,pcOffset)))
catch
    warning("not plotting contour")
end

%%
try
figure
contourf(NbPlot',ILPlot',NCCMean',15,"ShowText",true,"LabelFormat","%0.3f", ...
    "FaceAlpha",0.75)
ylabel('Interlayerspacing')
xlabel('3a fraction Nb')
clim([-0.03 0])
hold on
scatter(NbPlot,ILPlot,'black','filled')
colorbar
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_Flat_climFixed.png',AnlgeOffset,pcOffset)))

% exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.pdf',AnlgeOffset,pcOffset)))
% saveas(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f',AnlgeOffset,pcOffset)),'svg')%,'ContentType','vector')
% figure
% plot(allMetrics)

%make an exprot table
% exportTable()
catch
    warning("not plotting contour")
end
%%
results.InputUser = InputUser;
% results.occPlot = occPlot
% results.dtPlot = dtPlot
results.ILPlot = ILPlot;
results.NbPlot = NbPlot;
results.gmatrixFilderedAllCombined = gmatrixFilderedAllCombined;
results.eangs_outAllCombined = eangs_outAllCombined;
results.allNumPixels = allNumPixels;
% results.plottingIndex = newIndex;
results.allMetricsCombined = allMetricsCombined;
results.NCCMean = NCCMean;
results.DiffValues = DiffValues;
results.patternNums = patternNums;
% results.gmatrix = gmatrix
writeStructToH5(results,fullfile(exportFolder,sprintf("OriCorrectedRefinement_oriOFF%.3fdeg_pcOFF%.3f.h5",AnlgeOffset,pcOffset)),"");
% 
%     end
% end

%%
writeStructToH5(results,InputUser.h5refinement,sprintf("RTM/OriCorrectedRefinement/%i/oriOFF%.3fdeg_pcOFF%.3f",tileNum,AnlgeOffset,pcOffset));

end

%%

for i = 1:numel(InputUser.Phase_Input)
disp(InputUser.Phase_Input{i})
end
%%
function EBSD_ToTile =  removeCross512(patOriginal)
    EBSD_ToTile = zeros(512,512);
    EBSD_ToTile(1:256,1:256) =patOriginal(1:256,1:256);
    EBSD_ToTile(1:256,257:end) =patOriginal(1:256,259:end);
    EBSD_ToTile(257:end,1:256) =patOriginal(259:end,1:256);
    EBSD_ToTile(257:end,257:end) =patOriginal(259:end,259:end);
end
function ebsdPat = fix_mean(ebsdPat)
    ebsdPat = ebsdPat-mean(ebsdPat(:));
    ebsdPat = ebsdPat./std(ebsdPat(:));
end
function rot =  filterRotation(gmatrices)

    % gmatrix = mean(rotIn);
    rotIn = rotation.byMatrix(gmatrices);
%%
    indices = 1:numel(rotIn);
    indices = indices(abs(rotIn.phi1-median(rotIn.phi1))<5*degree);
    rotIn = rotation.byMatrix(gmatrices(:,:,indices));
    indices = indices(abs(rotIn.Phi-median(rotIn.Phi))<5*degree);
    rotIn = rotation.byMatrix(gmatrices(:,:,indices));
    indices = indices(abs(rotIn.phi2-median(rotIn.phi2))<5*degree);
    rotIn = rotation.byMatrix(gmatrices(:,:,indices));
    rot = mean(rotIn);
    %%
end

% function rot =  filterRotation(gmatrices)
% 
%     % gmatrix = mean(rotIn);
%     rotIn = rotation.byMatrix(gmatrices);
% %%
%     indices = 1:numel(rotIn);
%     indices = indices(abs(rotIn.phi1-median(rotIn.phi1))<5*degree);
%     rotIn = rotation.byMatrix(gmatrices(:,:,indices));
%     indices = indices(abs(rotIn.Phi-median(rotIn.Phi))<5*degree);
%     rotIn = rotation.byMatrix(gmatrices(:,:,indices));
%     indices = indices(abs(rotIn.phi2-median(rotIn.phi2))<5*degree);
%     rotIn = rotation.byMatrix(gmatrices(:,:,indices));
%     rot = mean(rotIn);
%     %%
% end
