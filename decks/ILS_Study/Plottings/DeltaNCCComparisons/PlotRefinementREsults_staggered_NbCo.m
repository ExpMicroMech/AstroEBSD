%% 17.10.2025 L. Berners
% In this skript we want to test the influence of potential misalignments
% on the comparison metrics between individual on-screen EBSD patterns.

% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
InputUser.Phase_Folder=Phase_Folder
InputUser.Astro_loc=Astro_loc
warning('off','MATLAB:MKDIR:DirectoryExists')

% exportFolder = 'C:\Users\berners\sciebo\Paper6_OccupancyAndEffectsInNbNi\CompareMisalignment'
% exportFolder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\58Nb-42Ni_v2\Area7"

%%
InputUser.h5folder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\51Nb-49Ni_v2\Area10\Map1\NbNi3aILGrid_LPR51Nb_49Ni_2026_01\cleaned"
% InputUser.h5folder = "J:\Users\Berners_Lukas\clara\NbCoUBC\56Nb44Co_verticalCut\Area3\Map1"
exportFolder = fullfile(InputUser.h5folder,'compareNCC')
mkdir(exportFolder)
InputUser.h5name  =  fullfile(InputUser.h5folder, "Minimization_1pats_tesOpti_1Ori.hdf5");
InputUser.h5refinement = fullfile(InputUser.h5folder, "NbNi_Nb-49Ni_V2_september_Area 10_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__NbNi3aILGrid_LPR51Nb_49Ni_2026_01_cleaned.hdf5");
% InputUser.h5refinement = fullfile(InputUser.h5folder, "NbCoUBC_56Nb44Co_vertical_Area 3_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__cleaned.hdf5")
InputUser.h5refinementPatterns = InputUser.h5refinement ;%strrep(InputUser.h5refinement,'_cleaned','')
h5structure = h5info(InputUser.h5refinement,'/RTM')%,'/Patterns','min')
%%

InfoStruct = readStructFromH5(InputUser.h5refinement,'skipPatterns')


%% Set up file names
InputUser.Phase_Input = {};
countNum = 1
InputUser.Nb = [];
% InputUser.Co = [];
% InputUser.a_lattice = [];
% InputUser.c_lattice = [];
InputUser.ILSpacing = []; 
for phaseNum = 1:numel(h5structure.Groups)
    phaseName = h5structure.Groups(phaseNum).Name;
    phaseName = phaseName(6:end); % drop the /RTM bit
    if strcmp(phaseName,'RTMPCA') || strcmp(phaseName,'tileInfo') || strcmp(phaseName,"RTMRefined")
        continue
    end
    InputUser.Phase_Input{countNum} = string(phaseName);

    
    % matches = regexp(InputUser.Phase_Input{countNum} , pattern, 'tokens');
    nb_number = [];
    % co_number = [];
    nb_number = regexp(InputUser.Phase_Input{countNum} , '3aNb-([\d\.]+)', 'tokens')
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
        InputUser.Nb = [ InputUser.Nb; str2num(nb_number)]
        % InputUser.Co = [ InputUser.Co ; str2num(co_number) ]
        % InputUser.a_lattice = [InputUser.a_lattice; str2num(a_lattice)]
        % InputUser.c_lattice = [InputUser.c_lattice; str2num(c_lattice)]
        InputUser.ILSpacing =  [InputUser.ILSpacing; str2num(ILSpacing)]
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

refinementString ='RTMRefined'
tileNum = 1
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

SettingsCor.screensize=514
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
rtmConfigRefine = defaultRTMConfig(rtmOptions)
rtmConfigRefine.generalOptions.AlignEdax =1
% rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start



MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfigRefine    )
%%
rtmConfigRefine.gaussBlur=0;
rtmConfigRefine.gaussianRange=1;
%%
rtmConfigRefine = defaultRTMConfig();
rtmConfigRefine.Settings_Cor.screensize = 514;
rtmConfigRefine.Settings_Cor.size = 514;
PC= [0.5 0.25 0.7]
rtmConfigRefine.PC = PC;
%% Read data and build synthetic pattern
pNum = 2
patternNumber = 1

patOriginal = h5read(InputUser.h5refinement  ,"/RTM/RTMPCA/1/PCAPatterns");
% patOriginal = patOriginal(:,:,)
PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');
gmatrixCommercialDetectorOrig = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/G_Out',InputUser.Phase_Input{pNum},refinementString,tileNum));
% gmatrixCommercialDetectorOrig =pagemtimes( pagemtimes( pagetranspose(gmatrixRefined),InputUser.Detector_tilt) ...
                                    % , InputUser.EdaxConversionMatrix )  ;


geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber));
rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));
gmatrix = gmatrixCommercialDetectorOrig(:,:,patternNumber)

patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
patSim = EBSP_BGCor(patSim,SettingsCor);

EBSD_To_TileSim = removeCross512(patSim);
EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));

p1Name = InputUser.Phases{1}
MultiPhaseStruct(1).screen_int
pNumCompare = 1
metricName = 'NCC'
metrics = zeros(numel(rot_origIn),numel(InputUser.Phase_Input))
% EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);  
%% Recalculation of orientation and projection center based on average of mutliple phases.
% needs to be put into a function and cleaned up later on
%
%
%
 %
PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');
numPhases = numel(InputUser.Phases);
minPhaseNum =1; % exclude lit data from joubert.
numPatterns = size(PC_origIn,2)
gmatrixAll = zeros(3,3,numPatterns,numPhases);
nccAll = zeros(numPatterns,numPhases);
csCompare = InputUser.CSList{end}
PC_all = zeros(numPhases,numPatterns,3);

for pNum = 1:numPhases

    gmatrixCommercialDetectorOrig = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/G_Out',InputUser.Phase_Input{pNum},refinementString,tileNum));
    % gmatrixCommercialDetectorOrig =pagemtimes( pagemtimes( pagetranspose(gmatrixRefined),InputUser.Detector_tilt) ...
                                        % , InputUser.EdaxConversionMatrix )  ;
    % PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');
    % PC_all(pNum,:,:) = PC_origIn';
    
    % geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber));
    % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig);
    % gmatrix = gmatrixCommercialDetectorOrig(:,:,patternNumber);
    gmatrixAll(:,:,:,pNum)=gmatrixCommercialDetectorOrig;
      nccAll(:,pNum) = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{pNum},refinementString,tileNum));

end
%%
patternNumber = 1

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
maxPlotNum =numel(InputUser.Phase_Input)
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
PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');

metricName = "NCC"

% occPlot = reshape(InputUser.Occupancy,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
% dtPlot = reshape(InputUser.dtSpacing,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
ILPlot = (InputUser.ILSpacing)%,numel(unique(InputUser.Nb)),numel(unique(InputUser.c_lattice)))
NbPlot = (InputUser.Nb)%(unique(InputUser.Nb)),numel(unique
[~,sortIndex] = sort(InputUser.ILSpacing)
NbPlot = reshape(NbPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
ILPlot = reshape(ILPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
% [~,newIndex] = sortrows(occPlot)
figure
hold on
ColList =colorcube(numPatterns)
AnlgeOffset = 0
pcOffset = 0
% for AnlgeOffset = [0,0.05,0.1,0.2,0.5,1,2]
    % for pcOffset = [0,0.002,0.005,0.01,0.02,0.05,]
rotAxis = vector3d(0,1,0)
offset = rotation.byAxisAngle(rotAxis,AnlgeOffset*degree)
% for 
allMetricsCombined = []
patternNums = []
for patternNumber = 1:numPatterns
    geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber)+pcOffset);
    allMetrics = []
    ncc_test = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{1},refinementString,tileNum));
    if ncc_test(patternNumber)<0.6
        continue
    end
    %%
    % rotIn = rotation.byMatrix(gmatrixAll(:,:,patternNumber,:));
    gmatrix =  filterRotation(gmatrixAll(:,:,patternNumber,:))
    gmatrix = gmatrix*offset;
    gmatrix = gmatrix.matrix;
    patOriginal = h5read(InputUser.h5refinement ,"/RTM/RTMPCA/1/PCAPatterns");
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
            imwrite(flipud(normaliseForImwrite(patOriginal)),fullfile(exportFolder,sprintf('Pat_exp%.3fdeg_pcOFF%.3f_Flat_%i.png',AnlgeOffset,pcOffset,patternNumber)))
        end
    end
    % if mean(allMetrics>0.5)
        %%
        % figure
        NCCListPlot = reshape(allMetrics,numel(unique(NbPlot)),numel(unique(ILPlot)))
        % 
        % surf(occPlot(newIndex,:),dtPlot(newIndex,:),  NCCListPlot(newIndex,:)-max(NCCListPlot(:)),'FaceColor',ColList(patternNumber,:),'DisplayName',sprintf('Component %i',patternNumber))%,'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
        % hold on
        allMetricsCombined= cat(3,allMetricsCombined,NCCListPlot)
        patternNums = [patternNums,patternNumber]
    % end

end
%%
figure
hold on
ColList = colorcube(size(allMetricsCombined,3)+5)


DiffValues = [];
for patNum = 1: size((allMetricsCombined),3)
    % NCCListPlot = allMetricsCombined(:,:,patNum)
    NCCPlot =  allMetricsCombined(:,:,patNum)
    NCCPlot = NCCPlot
    DiffVals = NCCPlot-max(NCCPlot(:));

    DiffValues = cat(3,DiffValues,DiffVals)

end
    % NCCPlot = reshape(NCCPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
for patNum = 1: size((allMetricsCombined),3)
    DiffVals = DiffValues(:,:,patNum)
    surf(ILPlot,NbPlot,  DiffVals-min(DiffValues(:)),'FaceColor',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)))
end
NCCMean = mean(DiffValues,3);
% contour()
contour(ILPlot,NbPlot,NCCMean,15, ...
    "FaceAlpha",0.5,'HandleVisibility','off')
% exportTable = table()
legend('Location','southoutside','NumColumns',3)
ylabel('Occupancy, Nb')
xlabel('dt Spacing')
zlabel('NCC-maxNCC')
view(75,30)
clim([min(DiffValues(:)) 0 ])
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.png',AnlgeOffset,pcOffset)))

view(170,30)
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_rot.png',AnlgeOffset,pcOffset)))
% view(170,30)
%%
figure
contourf(NbPlot',ILPlot',NCCMean',15,"ShowText",true,"LabelFormat","%0.3f", ...
    "FaceAlpha",0.75)
ylabel('Interlayerspacing')
xlabel('3a fraction Nb')
hold on
scatter(NbPlot,ILPlot,'black','filled')
colorbar
exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f_Flat.png',AnlgeOffset,pcOffset)))
%%
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
%%
results.InputUser = InputUser
results.occPlot = occPlot
results.dtPlot = dtPlot
results.plottingIndex = newIndex;
results.allMetricsCombined = allMetricsCombined

writeStructToH5(results,fullfile(exportFolder,sprintf("OriCorrectedRefinement_oriOFF%.3fdeg_pcOFF%.3f.h5",AnlgeOffset,pcOffset)),"")
    % end
% end
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                   Trying to weight NCC based on differnce between
%                   simulations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

otAxis = vector3d(0,1,0)
offset = rotation.byAxisAngle(rotAxis,AnlgeOffset*degree)
% for 
allMetricsCombined = []
patternNums = []
for patternNumber = 1%:numPatterns
    geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber)+pcOffset);
    allMetrics = []
    ncc_test = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{1},refinementString,tileNum));
    if ncc_test(patternNumber)<0.6
        continue
    end
    %%
    % rotIn = rotation.byMatrix(gmatrixAll(:,:,patternNumber,:));
    gmatrix =  filterRotation(gmatrixAll(:,:,patternNumber,:))
    gmatrix = gmatrix*offset;
    gmatrix = gmatrix.matrix;
    patOriginal = h5read(InputUser.h5refinementPatterns ,"/RTM/RTMPCA/1/PCAPatterns");
    patOriginal = patOriginal(:,:,patternNumber);
    patOriginal = EBSP_BGCor(patOriginal,SettingsCor);
    patSims = zeros(numPhases-minPhaseNum,size(patOriginal,1),size(patOriginal,2));
    for pNum = minPhaseNum:numPhases
    
            % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));

        

    
        patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
        patSim = EBSP_BGCor(patSim,SettingsCor);
        patSims(pNum,:,:) = patSim;
        % EBSD_To_TileSim = removeCross512(patSim);
        % % EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));
        % allMetrics = [allMetrics getMatchMetric(patOriginal,...
                            % patSim,metricName)];
           

    
    end
    % if mean(allMetrics>0.5)
        %%

        x  = 0:0.05:2

        y = normpdf(x,1,0.4)
        figure
        plot(x,y)
        ylim([0,1])
%%
        stdPatSims = squeeze(std(patSims,[],1));
        % stdPatSimsCopy = stdPatSims;
        stdPatSimsCopy=normpdf(stdPatSims,1,0.4).*stdPatSimsCopy
        figure
        imagesc(stdPatSimsCopy)
        colorbar
        figure
        imagesc(normpdf(stdPatSims,1,0.6))
        gaussweightedstdPatSims = normpdf(stdPatSims,1,0.4)
        colorbar
        %%


%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Plotting the Joubert Phases 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PC_origIn = h5read(InputUser.h5refinement ,'/RTM/RTMPCA/1/PCFitted');

% metricName = "NCC"

% occPlot = reshape(InputUser.Occupancy,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
% dtPlot = reshape(InputUser.dtSpacing,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
% [~,newIndex] = sortrows(occPlot)

for AnlgeOffset = [0,0.05,0.1,0.2,0.5,1,2]
    for pcOffset = [0,0.002,0.005,0.01,0.02,0.05,]
figure
hold on
ColList =colorcube(numPatterns)

rotAxis = vector3d(0,1,0)
offset = rotation.byAxisAngle(rotAxis,AnlgeOffset*degree)
% for 
allMetricsCombined = []
patternNums = []
for patternNumber = 1:numPatterns
    geometry = EBSP_Gnom(SettingsCor,PC_origIn(:,patternNumber)+pcOffset);
    allMetrics = []
    ncc_test = h5read(InputUser.h5refinement ,sprintf('/RTM/%s/%s/%i/RefinedNCC',InputUser.Phase_Input{1},refinementString,tileNum));
    if ncc_test(patternNumber)<0.6
        continue
    end
    %%
    % rotIn = rotation.byMatrix(gmatrixAll(:,:,patternNumber,:));
    gmatrix =  filterRotation(gmatrixAll(:,:,patternNumber,:))
    gmatrix = gmatrix*offset;
    gmatrix = gmatrix.matrix;
    patOriginal = h5read(InputUser.h5refinementPatterns ,"/RTM/RTMPCA/1/PCAPatterns");
    patOriginal = patOriginal(:,:,patternNumber);
    patOriginal = EBSP_BGCor(patOriginal,SettingsCor);
    for pNum = 1:minPhaseNum
    
            % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));

        

    
        patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
        patSim = EBSP_BGCor(patSim,SettingsCor);
        
        % EBSD_To_TileSim = removeCross512(patSim);
        % % EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));
        allMetrics = [allMetrics getMatchMetric(patOriginal,...
                            patSim,metricName)];
           

    
    end
    % if mean(allMetrics>0.5)
        %%
        % figure
        % NCCListPlot = reshape(allMetrics,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
        % 
        plot(allMetrics(:)-max(allMetrics(:)),'Color',ColList(patternNumber,:),'DisplayName',sprintf('Component %i',patternNumber))%,'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')       % hold on
        allMetricsCombined= cat(3,allMetricsCombined,allMetrics)
        patternNums = [patternNums,patternNumber]
    % end

end
%%
figure
hold on
ColList = colorcube(size(allMetricsCombined,3)+5)
% for patNum = 1: size((allMetricsCombined),3)
%     NCCListPlot = allMetricsCombined(:,:,patNum)
%     surf(occPlot(newIndex,:),dtPlot(newIndex,:),  NCCListPlot(newIndex,:)-max(NCCListPlot(:)+min(NCCListPlot(newIndex,:)-max(NCCListPlot(:))),'FaceColor',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)))
% end
% exportTable = table()
legend('Location','southoutside','NumColumns',3)
xlabel('Occupancy, Nb')
ylabel('dt Spacing')
zlabel('NCC-maxNCC')
view(75,30)

exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.png',AnlgeOffset,pcOffset)))
% exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.pdf',AnlgeOffset,pcOffset)))
% saveas(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f',AnlgeOffset,pcOffset)),'svg')%,'ContentType','vector')
% figure
% plot(allMetrics)

%make an exprot table
% exportTable()
results.InputUser = InputUser
results.occPlot = occPlot
results.dtPlot = dtPlot
results.plottingIndex = newIndex;
results.allMetricsCombined = allMetricsCombined

writeStructToH5(results,fullfile(exportFolder,sprintf("OriCorrectedRefinement_oriOFF%.3fdeg_pcOFF%.3f.h5",AnlgeOffset,pcOffset)),"")
    end
end









        for pNum = minPhaseNum:numPhases
    %%
    SettingsCor.cross = 2
    SettingsCor.quadrants=1

            % rot_origIn = rotation.byMatrix(gmatrixCommercialDetectorOrig(:,:,patternNumber));

        

    
        patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
        patSim = EBSP_BGCor(patSim,SettingsCor);
        % patSims(pNum,:,:) = patSim;
        % EBSD_To_TileSim = removeCross512(patSim);
        % % EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));
        % allMetrics = [allMetrics getMatchMetric(patOriginal,...
                            % patSim,metricName)];
           

    % 
    % end
        f = figure
        f.Position = [100 100 1500 500]
        subplot(1,3,1)
        imagesc(patOriginal)
        clim([-3 3])
        colorbar
        axis equal xy off
        colormap gray
        hold on
        subplot(1,3,2)
        
        imagesc(patSim)
        colorbar
        axis equal xy off
        clim([-3 3])
        
        ax = subplot(1,3,3)
        imagesc(patOriginal-patSim)
        clim([-3 3])
        axis equal xy off
        % colormap(ax, returnCmapDiv())
        colorbar
        sgtitle(cleanNameLatex(InputUser.Phase_Input{pNum}))
        exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_dt_%.3f_occNb%.2f.png',dtPlot(pNum-minPhaseNum+1),occPlot(pNum-minPhaseNum+1))))
        close 
        end
        % colormap gray
        %%
        % figure
        NCCListPlot = reshape(allMetrics,numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
        % 
        % surf(occPlot(newIndex,:),dtPlot(newIndex,:),  NCCListPlot(newIndex,:)-max(NCCListPlot(:)),'FaceColor',ColList(patternNumber,:),'DisplayName',sprintf('Component %i',patternNumber))%,'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
        % hold on
        allMetricsCombined= cat(3,allMetricsCombined,NCCListPlot)
        patternNums = [patternNums,patternNumber]
    % end

end
% %%
% figure
% hold on
% ColList = colorcube(size(allMetricsCombined,3)+5)
% for patNum = 1: size((allMetricsCombined),3)
%     NCCListPlot = allMetricsCombined(:,:,patNum)
%     surf(occPlot(newIndex,:),dtPlot(newIndex,:),  NCCListPlot(newIndex,:)-max(NCCListPlot(:)),'FaceColor',ColList(patNum,:),'DisplayName',sprintf('Component %i',patternNums(patNum)))
% end
% % exportTable = table()
% legend('Location','southoutside','NumColumns',3)
% xlabel('Occupancy, Nb')
% ylabel('dt Spacing')
% zlabel('NCC-maxNCC')
% view(75,30)
% 
% exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.png',AnlgeOffset,pcOffset)))
% % exportgraphics(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f.pdf',AnlgeOffset,pcOffset)))
% % saveas(gcf,fullfile(exportFolder,sprintf('NCC_AllPhasesDT_Var_oriAveraged_oriOFF%.3fdeg_pcOFF%.3f',AnlgeOffset,pcOffset)),'svg')%,'ContentType','vector')
% % figure
% % plot(allMetrics)
% 
% %make an exprot table
% % exportTable()
% results.InputUser = InputUser
% results.occPlot = occPlot
% results.dtPlot = dtPlot
% results.plottingIndex = newIndex;
% results.allMetricsCombined = allMetricsCombined
% 
% writeStructToH5(results,fullfile(exportFolder,sprintf("OriCorrectedRefinement_oriOFF%.3fdeg_pcOFF%.3f.h5",AnlgeOffset,pcOffset)),"")





 %%



 % Remove Center portion

% figure
% imagesc(EBSD_To_Tile)
% axis equal
% axis xy
% colormap gray

tileSize =64
nTiles = floor(size(EBSD_To_Tile)/tileSize);
NCCTiles = zeros(size(EBSD_To_Tile));


        % ebsd_pat2 = EBSP_gen(EBSD_geometry,rot_origIn,MultiPhaseStruct(pNum).screen_int);
% Try to get the PC metric patchw
for row = 1:nTiles(1)
    for col = 1:nTiles(2)
        selectorx = tileSize*(row-1)+1:tileSize*(row);
        selectory = tileSize*(col-1)+1:tileSize*(col);
        metrics = getMatchMetric(EBSD_To_Tile(selectorx,selectory),...
                        EBSD_To_TileSim(selectorx,selectory),metricName);
        NCCTiles(selectorx,selectory) = metrics;
        % row
        % col
    end
end
%
figure
imagesc(NCCTiles)
colorbar
axis xy
axis equal
%%

tileSize = 128
nTiles = floor(size(EBSD_To_Tile)/tileSize);
rowIn = randi([1,nTiles(1)],12,1);
colIn = randi([1,nTiles(2)],12,1);
for pNum = [1,2]



NCCTiles = zeros(size(EBSD_To_Tile,1),size(EBSD_To_Tile,2));


patSim = EBSP_gen(geometry,gmatrix,MultiPhaseStruct(pNum).screen_int,MultiPhaseStruct(pNum).screen_int.isHex);
patSim = EBSP_BGCor(patSim,SettingsCor);

EBSD_To_TileSim = removeCross512(patSim);
EBSD_To_Tile = removeCross512(EBSP_BGCor(patOriginal(:,:,patternNumber),SettingsCor));


f =figure
f.Position = [100 100 1000 900];
% subplot(1,2,1)
h= tiledlayout(5,6,'TileSpacing','none')
for i = 1:10
    row = rowIn(i)
    col =colIn(i)
    ax = nexttile();
            selectorx = tileSize*(row-1)+1:tileSize*(row);
            selectory = tileSize*(col-1)+1:tileSize*(col);
    imagesc(EBSD_To_Tile(selectorx,selectory))
    clim([-3,3])
    axis equal
    axis xy
    axis off
    colormap(ax,'gray')
    % subplot(1,2,2)
    ax = nexttile();
    imagesc(EBSD_To_TileSim(selectorx,selectory))
    clim([-4,4])
    axis equal
    axis xy
    axis off
    colormap(ax,'gray')
    ax=nexttile()
    imagesc(fix_mean(EBSD_To_TileSim(selectorx,selectory))-fix_mean(EBSD_To_Tile(selectorx,selectory)))
    clim([-3,3])
    axis equal
    axis xy
    axis off
    colormap(ax,returnCmapDiv)
end
exportgraphics(gcf,fullfile(exportFolder,sprintf('PatchwiseComparison_%s.png',InputUser.Phase_Input{pNum})))
end
%%
selectorx = tileSize*(row-1)+1:tileSize*(row);
        selectory = tileSize*(col-1)+1:tileSize*(col);
        metrics = getMatchMetric(EBSD_To_Tile(selectorx,selectory),...
                        EBSD_To_TileSim(selectorx,selectory),metricName);
%%



        % ebsd_pat2 = EBSP_gen(EBSD_geometry,rot_origIn,MultiPhaseStruct(pNum).screen_int);
% Try to get the PC metric patchw
for row = 1:nTiles(1)
    for col = 1:nTiles(2)
        selectorx = tileSize*(row-1)+1:tileSize*(row);
        selectory = tileSize*(col-1)+1:tileSize*(col);
        selectory1 = tileSize*(col*2-1)+1:tileSize*(col*2);
        selectory2 = tileSize*(col*2)+1:tileSize*(col*2+1);
        % metrics = getMatchMetric(EBSD_To_Tile(selectorx,selectory),...
                        % EBSD_To_TileSim(selectorx,selectory),metricName);
        NCCTiles(selectorx,selectory1) = EBSD_To_Tile(selectorx,selectory);
        NCCTiles(selectorx,selectory2) = EBSD_To_TileSim(selectorx,selectory);
        % row
        % col
    end
end
%
figure
imagesc(NCCTiles)
colormap gray
axis equal
axis xy
colorbar
%%
figure
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)

for pNum = 2:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    histogram(metrics(:,pNum),[0.99:0.0005:1],'DisplayName',sprintf('%s',plName2))
    hold on
end
legend()
title(sprintf("%s vs.",plName1))
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.png',subFolder)))
%% Now try to vary the PC a bit get the "correlation without any deviations in PC and Ori
randRange =0.001
rtmConfigRefine.Settings_Cor.gaussBlur=1.5
offesets = zeros(3,numel(ori))
noiseRatios=zeros(1,numel(ori))
%%
% offesets
noise =19

for oriNum = 1:numel(ori)
    EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);
    gmatrix = ori(oriNum).matrix;
    ebsd_pat1 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNumCompare).screen_int);
    [ebsd_pat1,noiseRatio] = addNoisePoisson(ebsd_pat1, noise);
    ebsd_pat1 =EBSP_BGCor(ebsd_pat1,rtmConfigRefine.Settings_Cor);
    noiseRatios(oriNum) =noiseRatio;
    randOffset = - randRange + rand(1,3)*2*randRange;
    offesets(:,oriNum) = randOffset;
    PC_mod = PC+randOffset;
    EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC_mod);  
    
    for pNum = 1:numel( InputUser.Phase_Input)


        ebsd_pat2 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNum).screen_int);
        ebsd_pat2 = EBSP_BGCor(ebsd_pat2,rtmConfigRefine.Settings_Cor);
        %
        % figure
        % imagesc(ebsd_pat1)
        % figure
        % imagesc(ebsd_pat2)
        %
        metrics(oriNum,pNum) = getMatchMetric(ebsd_pat1,ebsd_pat2,metricName);
    end
end
% exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.png',subFolder)))
imwrite(rescale(ebsd_pat1),fullfile(exportFolder,sprintf('ExampleNoisePattern_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
imwrite(rescale(ebsd_pat2),fullfile(exportFolder,sprintf('ExampleComparePattern_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
%%
figure
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)

for pNum = 2:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    histogram(metrics(:,pNum),[0.97:0.001:1],'DisplayName',sprintf('%s',plName2))
    hold on
end
legend()
title(sprintf("%s vs.",plName1))
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
%%
figure
offsetsNormed = sqrt(offesets(1,:).^2+ offesets(2,:).^2+ offesets(3,:).^2)
[offsetsNormed_sorted,indexSorted] = sort(offsetsNormed)
for pNum = 2:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plot(offsetsNormed_sorted,metrics(indexSorted,pNum),'DisplayName',sprintf('%s',plName2))
    hold on
end
legend('Location','southwest')
title(sprintf("%s vs.",plName1))
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('DirectComparison_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
%%
f=figure
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
f.Position = [100,100,400,500]
for pNum = 3:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral);
    estimate = metrics(:,pNum)
    estimate = estimate(estimate>(estimate-0.3*std(estimate)))
    [~,fp,xfp] = kde(estimate');
    % histogram(metrics(:,pNum),[0.97:0.001:1],'DisplayName',sprintf('%s',plName2))
    scatter(xfp,fp/max(fp),'filled','DisplayName',sprintf('%s',plName2),'LineWidth',3)
    hold on
end
legend('Location','southoutside')
title(sprintf("%s vs.",plName1))
xlabel('NCC')
ylabel('probability density')
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_%s_pcVar%.4f_noise_%.2f.png',subFolder,randRange,noise)))
%%
f=figure
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
f.Position = [100,100,400,500]
comparePhase=2
for pNum = 3:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral);
    estimate = metrics(:,pNum) - metrics(:,comparePhase) 
    estimate = estimate(estimate>(estimate-0.3*std(estimate)))
    [~,fp,xfp] = kde(estimate');
    % histogram(metrics(:,pNum),[0.97:0.001:1],'DisplayName',sprintf('%s',plName2))
    scatter(xfp,fp,'filled','DisplayName',sprintf('%s',plName2),'LineWidth',3)
    xlim([-0.01,0])
    hold on
end
legend('Location','southoutside')
comparePhaseName = cleanNameLatex(InputUser.CSList{comparePhase}.mineral)
title(sprintf("%s vs. \n  %s -",plName1,comparePhaseName))
xlabel('NCC')
ylabel('probability density')
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))

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
    indices = 1:numel(rotIn)
    indices = indices(abs(rotIn.phi1-median(rotIn.phi1))<5*degree)
    rotIn = rotation.byMatrix(gmatrices(:,:,indices))
    indices = indices(abs(rotIn.Phi-median(rotIn.Phi))<5*degree)
    rotIn = rotation.byMatrix(gmatrices(:,:,indices))
    indices = indices(abs(rotIn.phi2-median(rotIn.phi2))<5*degree)
    rotIn = rotation.byMatrix(gmatrices(:,:,indices))
    rot = mean(rotIn);
    %%
end
