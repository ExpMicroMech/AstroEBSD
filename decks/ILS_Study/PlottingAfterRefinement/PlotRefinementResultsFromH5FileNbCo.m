%% L. Berners (25.08.2025),In this deck we will  load various parameters from a .h5 file
% that was produced using the RTM_PCA_COMMERCIALHOUGH_map.m script
% We will have a look at how the NCC value changes depending on the phsae
% and how this depends on the orietatioin.
%
%
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|(MTEX_loc)|(ORII5_loc)).

%
% clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|(MTEX_loc)|(ORII5_loc)).

InputUser.Phase_Folder = Phase_Folder;
% InputUser.H5Folder=  "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\51Nb-49Ni_v2\EBSD\Area6"
InputUser.H5Folder = "J:\Users\Berners_Lukas\clara\NbCoUBC\48Nb52Co_verticalCut\Area2\Map1"
% InputUser.H5File = "SinglePatRefinementWithLPR_Backup.hdf5"
InputUser.H5File = "NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area 2_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__cleaned.hdf5";
RefinementHDF5File = fullfile(InputUser.H5Folder,InputUser.H5File);
h5structure = h5info(RefinementHDF5File,'/RTM')%,'/Patterns','min')
%%

InfoStruct = readStructFromH5(RefinementHDF5File,'skipPatterns')
%% Get the phase names from the struct fields. We will need this to load the phases correctly
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
% InputUser.Phase_Input = vertcat(h5structure.Groups.Name)
%%
% text = 'NbNi_dt0.250_c26.81_a4.92_3aNb_19.3kV_dmin0.020'; % Your input text

% Regular expression to match the required patterns
% Extracting the matches

%% Load the crystal strcutures
for i = 1:numel(InputUser.Phase_Input)
    if strcmp(InputUser.Phase_Input{i},"tileInfo")
        continue
    end
     InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    h5PhaseInfo = Phase_BuilderSimplified(InputUser.Phase_Input{i},InputUser.Phase_Folder)
    charName = char(h5PhaseInfo.cif);
    h5PhaseInfo.PlotName =  charName; % change this if data not in .pha file TODO: use if else
    InputUser.h5PhaseInfo{i}=  h5PhaseInfo  
    % load ebsd map
    rtmConfig = defaultRTMConfig()
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end



%% Setup empty variables to read data from h5 file.
tileNum = 1
InputUser.Phase_Input
numPhases = numel(InputUser.Phase_Input)

phase = h5read(RefinementHDF5File,sprintf('/RTM/RTMPCA/%i/Phase',tileNum));
patternNumber = numel(phase);
NCCListLattice = zeros(numel(InputUser.Phase_Input),patternNumber);
NCCList = zeros(numel(InputUser.Phase_Input),patternNumber);
oris = cell(numel(InputUser.Phase_Input),1);
cAngle = zeros(numel(InputUser.Phase_Input),patternNumber);
NCCListFFT = zeros(numel(InputUser.Phase_Input),patternNumber);

orisFFT = cell(numel(InputUser.Phase_Input),1);

orisLattice = cell(numel(InputUser.Phase_Input),1);

cAngleLattice = zeros(numel(InputUser.Phase_Input),patternNumber);

cAngleFFT = zeros(numel(InputUser.Phase_Input),patternNumber);
%%

for phaseNum = 1:numPhases
    try
        
        NCCValsLattice = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefinedLattice/%i/RefinedNCC',InputUser.Phase_Input{phaseNum},tileNum));
        NCCListLattice(phaseNum,:) =NCCValsLattice;
        EangsLattice = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefinedLattice/%i/eangs_out',InputUser.Phase_Input{phaseNum},tileNum));
        % if size(EangsLattice(:,2) == patternNumber)
        EangsLattice=EangsLattice(:,1:3)
        % end
        oriLattice = orientation.byEuler(EangsLattice,InputUser.CSList{i});
        cAngleLattice(phaseNum,:) = angle(oriLattice,vector3d.Z)/degree;
        orisLattice{phaseNum} = oriLattice;
    catch
        warning('Lattice Refinement not read')
    end
     try
        
        NCCVals = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/RefinedNCC',InputUser.Phase_Input{phaseNum},tileNum));
        NCCList(phaseNum,:) =NCCVals;
        Eangs = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/eangs_out',InputUser.Phase_Input{phaseNum},tileNum));
        % if size(EangsLattice(:,2) == patternNumber)
        Eangs=Eangs(1:3,:)
        % end
        ori = orientation.byEuler(Eangs',InputUser.CSList{i});
        cAngle(phaseNum,:) = angle(ori,vector3d.Z)/degree;
        oris{phaseNum} = ori;
    catch
        warning('normal Refinement not read')
    end
    try
        NCCValsFFT = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMPCA/%i/Refined_PH',InputUser.Phase_Input{phaseNum},tileNum));
        NCCListFFT(phaseNum,:) = NCCValsFFT;
        EangsFFT = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMPCA/%i/eangs_out',InputUser.Phase_Input{phaseNum},tileNum));
        % if size(EangsLattice(:,2) == patternNumber)
        EangsFFT =EangsFFT(:,1:3);
        % end
        oriFFT = orientation.byEuler(EangsFFT,InputUser.CSList{i});
        cAngleFFT(phaseNum,:) = angle(oriFFT,vector3d.Z)/degree;
        orisFFT{phaseNum} = oriFFT;

     catch
        warning('FFT Refinement not read')
    end
end

%%
fNames = fieldnames(InfoStruct.RTM);
% InfoStruct.RTM.(InputUser.Phase_Input{1})
%% Plot the results from FFT based refinemend
% [cmap_blues]=cbrewer('div','PuOr',numPhases);
% cmap_blues = flipud(cmap_blues);
% cmap_blues(cmap_blues<0)=0;
% figure
% hold on
% Markers = ['s','p','^','<','>','v','o','d']
% for phaseNum = 5:numPhases-1
%     scatter(cAngle(phaseNum,:),NCCList(phaseNum,:),100,'DisplayName',InputUser.h5PhaseInfo{phaseNum}.PlotName,'MarkerFaceColor',cmap_blues(phaseNum,:),'MarkerEdgeColor','black')
% end
% ylim([0.75,0.8])
% legend('location','eastoutside')


%%
%% Plot the results from refinement of the Lattic eParameters
% [cmap_blues]=cbrewer('div','PuOr',numPhases);
% cmap_blues = flipud(cmap_blues);
% cmap_blues(cmap_blues<0)=0;
% figure
% hold on
% Markers = ['s','p','^','<','>','v','o','d']
% for phaseNum = 1:numPhases
%     % ifpKey = ipfColorKey(InputUser.CSList{1})
%     % plotColors = ifpKey.orientation2color(orisFFT{1})
%     scatter(1:patternNumber,NCCList(phaseNum,:),100,'DisplayName',cleanNameLatex(InputUser.h5PhaseInfo{phaseNum}.PlotName))%'Marker',Markers(phaseNum))%'MarkerFaceColor',cmap_blues(phaseNum,:),'MarkerEdgeColor','black')
% end
% % ylim([0.75,0.8])
% xlabel('Pattern Number')
% ylabel('NCC')
% legend('location','eastoutside')
% exportgraphics(gcf,fullfile(InputUser.H5Folder,'RefinementResults_NCCPatternNumber.png'))
% legend('location','eastoutside')
%%

%%

ipfKey = ipfColorKey(InputUser.CSList{1});
% ipfKey = rescale(ipfKey)
plotColors = ipfKey.orientation2color(oris{1});
plotColors = rescale(plotColors)
%% Only works with the interlayer spacing refinement: TODO: set up to be more autoamted
% maxPlotNum=numPhases
% 
% plot(ipfKey)
% [cmap_blues]=cbrewer('div','PuOr',numPhases);
% cmap_blues = flipud(cmap_blues);
% cmap_blues(cmap_blues<0)=0;
% figure
% hold on
% % Markers = ['s','p','^','<','>','v','o','d']
% 
% for pNum = 1:patternNumber
%     if mean(NCCList(1:maxPlotNum,pNum))>0.7 & abs(std(NCCList(1:maxPlotNum,pNum)))<0.05
%     ifpKey = ipfColorKey(InputUser.CSList{pNum});
%     plotColors = ifpKey.orientation2color(oris{1});
%     plot(InputUser.dtSpacing,  NCCList(1:maxPlotNum,pNum)-max(NCCList(1:maxPlotNum,pNum)),'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
%     end
% end
% ylim([-0.03,0])
% xlabel("interlayer spacing d_t")
% ylabel('NCC')
% exportgraphics(gcf,fullfile(InputUser.H5Folder,'RefinementResults_NCC.png'))
% legend('location','eastoutside')
%%
maxPlotNum=numPhases
minPlotNum = 6
% figure
% plot(ipfKey)
[cmap_blues]=cbrewer('div','PuOr',numPhases);
cmap_blues = flipud(cmap_blues);
cmap_blues(cmap_blues<0)=0;
figure
% hold on
% Markers = ['s','p','^','<','>','v','o','d']
% 
% lattice_Filter = 26.2
% filterIndex = InputUser.c_lattice > lattice_Filter
ILPlot = (InputUser.ILSpacing)%,numel(unique(InputUser.Nb)),numel(unique(InputUser.c_lattice)))
NbPlot = (InputUser.Nb)%(unique(InputUser.Nb)),numel(unique
[~,sortIndex] = sort(InputUser.ILSpacing)
NbPlot = reshape(NbPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
ILPlot = reshape(ILPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
% NCCPlot = 
% (InputUser.c_lattice)))
hold on
for pNum = 1:patternNumber
    % NCCListPlot = reshape(NCCList(minPlotNum:maxPlotNum,pNum),numel(unique(InputUser.Occupancy)),numel(unique(InputUser.dtSpacing)))
    if mean(NCCList(minPlotNum:maxPlotNum,pNum))>0.6 & abs(std(NCCList(minPlotNum:maxPlotNum,pNum)))<0.05
        NCCPlot = NCCList(:,pNum)
        NCCPlot = NCCPlot
        NCCPlot = reshape(NCCPlot(sortIndex), numel(unique(NbPlot)), numel(unique(ILPlot)))
    % ifpKey = ipfColorKey(InputUser.CSList{pNum});
    % plotColors = ifpKey.orientation2color(oris{1});
    surf(ILPlot, NbPlot,  NCCPlot)%-max(NCCPlot(:)))%,'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
    % scatter3(cPlot,NbPlot,NccPlot)
    end
end
% ylim([-0.03,0])
xlabel("IL, \angstroms")
ylabel('Nb')
view(30,20)
exportgraphics(gcf,fullfile(InputUser.H5Folder,'RefinementResults_NCC.png'))
% legend('location','eastoutside')

%%
% y%%
figure
plot(ipfKey)
[cmap_blues]=cbrewer('div','PuOr',numPhases);
cmap_blues = flipud(cmap_blues);
cmap_blues(cmap_blues<0)=0;
figure
hold on
Markers = ['s','p','^','<','>','v','o','d']

for pNum = 1:patternNumber
    if mean(NCCList(1:maxPlotNum,pNum))>0.7 & abs(std(NCCList1:maxPlotNum,pNum))<0.05
    ifpKey = ipfColorKey(InputUser.CSList{1})
    plotColors = ifpKey.orientation2color(orisFFT{1})
    plot(1:numPhases,NCCList(1:maxPlotNum,pNum)-max(NCCList(1:maxPlotNum,pNum)),'Marker',Markers(phaseNum),'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
    end
end
% ylim([-0.03,0])
xlabel("interlayer spacing d_t")
ylabel('NCC')
%%
% ticklabels = []
% for i = 1: numel(InputUser.h5PhaseInfo)
%     name = InputUser.h5PhaseInfo{i}.xtalname
%     name = char(name)
%    ticklabels = [ticklabels,string(name(end-3:end)) ] 
%     % xtalname: "NbNi_dt0.400_c26.7_a4.9"
% end
%     xticks([1,2,3])
% xticklabels(ticklabels')
% % exportgraphics(gcf,fullfile(InputUser.H5Folder,'RefinementResults_NCC.png'))

% legend('location','eastoutside')
%%

[cmap_blues]=cbrewer('div','PuOr',numPhases);
cmap_blues = flipud(cmap_blues);
cmap_blues(cmap_blues<0)=0;
figure
hold on
Markers = ['s','p','^','<','>','v','o','d']

for pNum = 1:patternNumber
    if mean(NCCListLattice(1:maxPlotNum,pNum))>0.7 & abs(std(NCCListLattice(1:maxPlotNum,pNum)))<0.05
    ifpKey = ipfColorKey(InputUser.CSList{1})
    plotColors = ifpKey.orientation2color(orisFFT{1})
    plot(1:numPhases,NCCListLattice(1:maxPlotNum,pNum)-max(NCCListLattice(1:maxPlotNum,pNum))+pNum/200,'Marker',Markers(phaseNum),'MarkerFaceColor',plotColors(pNum,:),'MarkerEdgeColor','black','Color','black')
    end
end
% ylim([-0.03,0])
xlabel("interlayer spacing d_t")
ylabel('NCC')
exportgraphics(gcf,fullfile(InputUser.H5Folder,'RefinementResults_NCC_Stacked.png'))

%% Build EBSD map
% for tNum = 1:4
numTiles = numel(fieldnames(InfoStruct.RTM.RTMPCA));
cropFactor = sqrt(numTiles);
% % try for now just to scale up
% for tNum = 1:numTiles
tNum=1
tile = InfoStruct.RTM.RTMPCA.(sprintf("t%i",tNum))
fullSize = size(tile.XVals)*sqrt(numTiles);
singleSize = size(tile.XVals)
fullSize = singleSize
% end

phaseList = zeros(fullSize);
eangsList = zeros(fullSize(1),fullSize(2),3);
LPRxList=  double(zeros(fullSize(1),fullSize(2),3));
PCList=  zeros(fullSize(1),fullSize(2),3);
posList = zeros(fullSize(1),fullSize(2),3);
%%

posAll = [ ]
rotAll = [ ]
phaseListAll = []
propAll = [ ]
NCCTilesAll = zeros()
NCCListAll = []
CompListAll = []
for tNum = 1:numTiles

phaseList = zeros(fullSize);
eangsList = zeros(fullSize(1),fullSize(2),3);
LPRxList=  double(zeros(fullSize(1),fullSize(2),3));
PCList=  zeros(fullSize(1),fullSize(2),3);
posList = zeros(fullSize(1),fullSize(2),3);
% CompList = zeros(fullSize);

[row,col] = ind2sub([cropFactor,cropFactor],tNum);
  % x = tNum*cropFactor + floor(mod(tNum,cropFactor))
  % y =  floor(mod(tNum/cropFactor)
  tile =  InfoStruct.RTM.RTMPCA.(sprintf("t%i",tNum));
  posList(singleSize(1)*(row-1)+1:singleSize(1)*row,singleSize(2)*(col-1)+1:singleSize(2)*col,1)=  tile.XVals;
  posList(singleSize(1)*(row-1)+1:singleSize(1)*row,singleSize(2)*(col-1)+1:singleSize(2)*col,2)=  tile.YVals;
% end

%
% for tNum = 1:numTiles
% tile = InfoStruct.RTM.RTMPCA.(sprintf("t%i",tNum))
zVals = ones(size(tile.YVals));
zVals = zVals(:);
pos = vector3d( double(tile.XVals(:)), double(tile.YVals(:)),double(zVals(:)));

%
patternNumber
phase = h5read(RefinementHDF5File,sprintf('/RTM/RTMPCA/%i/Phase',tNum));
patternNumber = numel(phase);
NCCList = zeros(numel(InputUser.Phase_Input),patternNumber);
oris = cell(numel(InputUser.Phase_Input),1);
cAngle = zeros(numel(InputUser.Phase_Input),patternNumber);
for phaseNum = 1:numPhases
  try
        
        NCCVals = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/RefinedNCC',InputUser.Phase_Input{phaseNum},tNum));
        NCCList(phaseNum,:) =NCCVals;
        Eangs = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/eangs_out',InputUser.Phase_Input{phaseNum},tNum));
        % if size(EangsLattice(:,2) == patternNumber)
        Eangs=Eangs(1:3,:);
        % end
        ori = orientation.byEuler(Eangs',InputUser.CSList{i});
        cAngle(phaseNum,:) = angle(ori,vector3d.Z)/degree;
        oris{phaseNum} = ori;
    catch
        warning('normal Refinement not read')
    end
end


[maxNCC,phaseId] = max(NCCList);
phaseId(maxNCC<0.45)=0;
NCCAllPhases = zeros(fullSize(1),fullSize(2),numPhases);
% NCCAllPhases
% oriList = orientation(eangsZero)
% REBuild the correc
% t eang.s the edax conversion matrix did not work correctly
EdaxConversionMatrix =  rtmConfig.RTM.Ry(180*degree); 
Detector_tilt =  rtmConfig.RTM.Rx(InfoStruct.MicroscopeData.TotalTilt);

%
RefinementKey = "RTMRefined";
% defaultRTMConfig
% tNum = 1
% tile.components(tile.components==0) =1% (catch indexation error)
tile.components(tile.components==0) =1;% (catch indexation error)
for patNum = 1:numel(tile.Phase);
    pId = phaseId(patNum);
    if pId == 0
        continue
    end
 
    subselection = tile.components==patNum;
    selectorsX = posList(:,:,1);
    selectorsY = posList(:,:,2);
    selection = (selectorsX(:)==tile.XVals(subselection)' & selectorsY(:)==tile.YVals(subselection)' );
    selection = subselection;
    phaseList(selection)=phaseId(patNum);
    ori = orisLattice(pId);
    % probably not super efficien
    EangsLattice = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/eangs_out',InputUser.Phase_Input{pId},tNum));
    GoutLattice = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/G_Out',InputUser.Phase_Input{pId},tNum));
    % EangsLattice = h5read(RefinementHDF5File,sprintf('/RTM/RTMPCA/%i/EangsCommercial',tileNum));
    % LPR = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/LatticeRatio',InputUser.Phase_Input{pId},tileNum));
    PC = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/PC',InputUser.Phase_Input{pId},tNum));
    % NCCVals = h5read(RefinementHDF5File,sprintf('/RTM/%s/RTMRefined/%i/RefinedNCC',InputUser.Phase_Input{phaseNum},tNum));
    % NCCListEBSD(selection) = NCCVals();

    Goutsel = GoutLattice(:,:,patNum);
    GReal =Goutsel *inv(EdaxConversionMatrix)*inv(Detector_tilt) ;
    rotIn = rotation.byMatrix(GReal');
    for i = 1:numPhases
        selectionNCCList = zeros(size(NCCAllPhases))==1;
        selectionNCCList(:,:,i)=selection;
        NCCAllPhases(selectionNCCList)=NCCList(i,patNum);
    end
    NCCListEBSD(selection) = NCCList(pId,patNum);
    CompList(selection) =patNum;
    % LPRxList(selection) = LPR(1,patNum);
    % % LPRyList(selection) = LPR(2,patNum);
    % LPRzList(selection) = LPR(3,patNum);
    PCxList(selection) = PC(1,patNum);
    PCyList(selection) = PC(2,patNum);
    PCzList(selection) = PC(3,patNum);
    eangsList(selection,1) = EangsLattice(1,patNum)';
    eangsList(selection,2) = EangsLattice(2,patNum)';
    eangsList(selection,3) = EangsLattice(3,patNum)';
    % eangsList(selection,1) = EangsLattice(1,patNum)';
    % eangsList(selection,2) = EangsLattice(2,patNum)';
    % eangsList(selection,3) = EangsLattice(3,patNum)';
    % oriList(selection) = ori.rotation;end
     eangsList(selection,1) = rotIn.phi1;
    eangsList(selection,2) = rotIn.Phi;
    eangsList(selection,3) = rotIn.phi2;
%
end

prop.NCC = NCCListEBSD;
prop.PCx = PCxList;
prop.PCy = PCyList;
prop.PCz = PCzList;
prop.Component = CompList;
rot=rotation.byEuler(eangsList(:,1),eangsList(:,2),eangsList(:,3))
% CompListAll = [CompListAll CompList];
posAll = [posAll pos];
rotAll = [rotAll rot];
phaseListAll = [phaseListAll phaseList];
propAll = [propAll prop];
NCCListAll = cat(4,NCCListAll, NCCAllPhases);
end
%%

%%
propAllClean =struct();
fn =fieldnames(propAll)
for fnum =  1:numel(fn)
    newVals = vertcat(propAll.(fn{fnum}))'
    propAllClean.(fn{fnum}) = newVals(:);
end
%
cs = crystalSymmetry("m3m")
InputUser.CSList2 = InputUser.CSList
% InputUser.CSList2{1} =cs
ebsdmap = EBSD(posAll(:),rotAll(:),phaseListAll,InputUser.CSList,propAllClean)
NCCListAllR = reshape(NCCListAll,numel(ebsdmap),numPhases);
%%
plotx2east
plotzIntoPlane
%%
exportName = char(InputUser.H5File);
exportName = exportName(1:end-6);
for i = 1:phaseNum
    phasePlot =i
    figure
    
    ipfKey = ipfHSVKey(ebsdmap(ebsdmap.phase==phasePlot).orientations);
    plot(ipfKey)
    exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%sipfKey_%i.png',exportName,i)))
    % plot the key


end
%%
grains =calcGrains(ebsdmap('indexed'))
grains = smooth(grains,5);
figure
for i = 1:phaseNum

end
 


%%

for i = 1:phaseNum
    figure
    % gbPlot = grains(ebsdmap.phase==phasePlot).boundary
    % gbPlot()% ebsdMapTemp = ebsdmap(ebsdmap.phase==phasePlot)
    % ebsdMapTemp.orientations = repmat(orientation.byEuler(0,0,0),size(ebsdMapTemp))
    % grains = calcGrains(ebsdMapTemp('indexed'))
    % grains = smooth(grains,5);
        % grainBoundaryLocal = grains(grains.phase==i).boundary
    % grainBoundaryLocal =    grainBoundaryLocal(grainBoundaryLocal.phase(:,2)==i)
    % plot(grainBoundaryLocal,'lineWidth',12,'lineColor',InputUser.CSList{i}.color)
    % hold on
    % plot(grainBoundaryLocal,'lineWidth',2,'lineColor','black')
    % hold on
    % a
    phasePlot =i
    ipfKey = ipfHSVKey(ebsdmap(ebsdmap.phase==phasePlot).orientations);
    % plotx2east
    % plotzIntoPlane
    % now generate the IPF colour maps
    ipfKey.inversePoleFigureDirection = vector3d.X;
    colors_X = ipfKey.orientation2color(ebsdmap(ebsdmap.phase==phasePlot) .orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors_Y = ipfKey.orientation2color(ebsdmap(ebsdmap.phase==phasePlot).orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors_Z = ipfKey.orientation2color(ebsdmap(ebsdmap.phase==phasePlot).orientations);
    % figure;
    plotx2east
    plotzIntoPlane
    nextAxis(1)
    plot(ebsdmap(ebsdmap.phase==phasePlot),colors_X,'micronbar','off');
    hold on

    % plot(grains.boundary,'Color',"black",'lineWidht',6)
    % hold on
    % plot(grains.boundary,'Color',InputUser.CSList{i}.color,'lineWidth',4)

    hold on
    title('IPF X')
    nextAxis(2)
    hold on
    plot(ebsdmap(ebsdmap.phase==phasePlot),colors_Y,'micronbar','off');
    hold on
    % plot(grains.boundary,'Color',"black",'lineWidht',6)
    % hold on
    % plot(grains.boundary,'Color',InputUser.CSList{i}.color,'lineWidth',4)
    hold on
    % plot(grains.boundary,'FaceAlpha',0.5)
    title('IPF Y')
    nextAxis(3)
    hold on
    plot(ebsdmap(ebsdmap.phase==phasePlot),colors_Z,'micronbar','on');
    % hold on
    % plot(grains.boundary,'Color',"black",'lineWidht',6)
    % hold on
    % plot(grains.boundary,'Color',InputUser.CSList{i}.color,'lineWidth',4)
    hold on
    % plot(grains,'FaceAlpha',0.5)
    title('IPF Z')
    % a
    hold on
    % plot(ipfKey)    
    if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%sEBSDmap_%i.png',exportName,i)))
end


% end
%%


%% Plot the LPR ratios
% 
% figure;
% plotx2east
% plotzIntoPlane
% lprRange = 0.02
% 
% figure
% % plot(ebsdmap,ebsdmap.prop.PCx)
% 
% 
% plot(ebsdmap,ebsdmap.prop.PCx,'micronbar','off');
% title('PCX')
% colorbar()
% colormap(returnCmapDiv)
% % clim([1-lprRange,1+lprRange])
% 
% nextAxis
% plot(ebsdmap,ebsdmap.prop.PCy,'micronbar','off');
% title('PCY')
% colorbar()
% colormap(returnCmapDiv)
% % clim([1-lprRange,1+lprRange])
% 
% nextAxis
% plot(ebsdmap,ebsdmap.prop.PCz,'micronbar','off');
% title('PCZ')
% colorbar()
% colormap(returnCmapDiv)
% % clim([1-lprRange,1+lprRange])
% 
% exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%s_PC.png',exportName)))
% figure;
% plotx2east
% plotzIntoPlane
% lprRange = 0.02
%%
figure
% plot(ebsdmap,ebsdmap.prop.PCx)


plot(ebsdmap,ebsdmap.prop.LPRx,'micronbar','off');
title('LPRx')
colorbar()
colormap(returnCmapDiv)
clim([1-lprRange,1+lprRange])

nextAxis
plot(ebsdmap,ebsdmap.prop.LPRy,'micronbar','off');
title('LPRy')
colorbar()
colormap(returnCmapDiv)
clim([1-lprRange,1+lprRange])

nextAxis
plot(ebsdmap,ebsdmap.prop.LPRz,'micronbar','off');
title('LPRz')
colorbar()
colormap(returnCmapDiv)
clim([1-lprRange,1+lprRange])

exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%s_LPR.png',exportName)))
%%
figure
% plot(ebsdmap,ebsdmap.prop.PCx)


plot(ebsdmap,ebsdmap.prop.NCC,'micronbar','off');
title('NCC')
clim([0.4 0.7])
colorbar()
% colormap('winter')
% clim([1-lprRange,1+lprRange])

% nextAxis
%
% plot(ebsdmap,ebsdmap.prop.LPRy,'micronbar','off');
% title('PCY')
% colorbar()
% colormap(returnCmapDiv)
% clim([1-lprRange,1+lprRange])
% 
% nextAxis
% plot(ebsdmap,ebsdmap.prop.LPRz,'micronbar','off');
% title('PCZ')
% colorbar()
% colormap(returnCmapDiv)
% clim([1-lprRange,1+lprRange])

exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%s_NCCLPR.png',exportName)))
%%
[cmap_blues]=cbrewer('seq','Blues',30);
cmap_blues = flipud(cmap_blues);
cmap_blues(cmap_blues<0)=0;
NCCAllPhases = permute(NCCListAll,[1 2 4 3]);
NCCAllPhases = reshape(NCCAllPhases,numel(ebsdmap),numPhases)
for i = 1:phaseNum 
    for j = i:phaseNum
        if i == j
            continue
        end
    figure

    plot(ebsdmap,NCCAllPhases(:,i)-NCCAllPhases(:,j),'micronbar','on');
    pName1 = InputUser.CSList{i}.mineral
    pName2 = InputUser.CSList{j}.mineral
    title(sprintf('NCC %s-\n%s',cleanNameLatex( InputUser.Phase_Input{i}),cleanNameLatex(InputUser.Phase_Input{j})))
    colorbar()

colormap(returnCmapDiv)
clim([-0.05,0.05])
    exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%s-%s_%s.png',InputUser.Phase_Input{i},InputUser.Phase_Input{j},exportName)))
        close all
    end

end
figure
plot(ebsdmap)
exportgraphics(gcf,fullfile(InputUser.H5Folder,sprintf('%sPhasemap.png',exportName)))
%%
figure

% plot(ebsdmap.phaseId)
%%
% figure
% plot(ebsdmap,ebsdmap.KAM('threshold',2.5*degree) ./ degree,'micronbar','off')
% colorbar
% setColorRange([0,0.5])
%%
figure
histogram(ebsdmap.KAM('threshold',2.5*degree)./degree)