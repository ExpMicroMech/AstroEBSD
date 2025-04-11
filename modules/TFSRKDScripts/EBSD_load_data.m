%% EBSD_load_data.m
% Example script of loading orientation data from xTalView and plot maps in MTEX
% Tianbi Zhang and T. Ben Britton, April 2025

% Requirements:
% MTEX Toolbox
% project metadata site.json
% orientation data in results/multiphase.idx
% phase info in phases (N.B. correct phase data must be loased in xTalView
% first)

%% Load MTEX
InputUser.mtex_location='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';

% start mtex if needed
try EBSD;
catch
    run(fullfile(InputUser.mtex_location,"startup.m"));
end

%% Read map metadata and extract
mapLocation = "D:\ApreoTruePixData\2025-03-31T14_04_15_Tianbi_Ebru_SiC\2025-03-31T14_23_02";
mapJson = "site.json";

EBSDMapInfo = ReadTFSJson(mapLocation, mapJson);

mapWidth = EBSDMapInfo.Info.Width; % pixels
mapHeight = EBSDMapInfo.Info.Height; % pixels
mapStepSize = EBSDMapInfo.Info.MapStepSize * 1e6; % um

% other useful info you may need
% mapkeV = EBSDMapInfo.Info.LandingEnergy / 1e3; % keV
% mapBC = EBSDMapInfo.Info.BeamCurrent * 1e9; % nA
% mapTH = EBSDMapInfo.Info.Threshold / 1e3; % keV

%% Load phase data
mapPhases = ReadTFSJson(mapLocation, "phases/active.json");

for i=1:size(mapPhases.ActivePhasesList,1)
    thisPhaseJson = mapPhases.ActivePhasesList(i).File;

    thisPhaseInfo = ReadTFSJson(fullfile(mapLocation, 'phases'), thisPhaseJson);

    cs{i+1} = crystalSymmetry(thisPhaseInfo.laue_group_symbol, ...
        [thisPhaseInfo.a thisPhaseInfo.b thisPhaseInfo.c], ...
        [ thisPhaseInfo.alpha thisPhaseInfo.alpha thisPhaseInfo.alpha] * degree, ...
        'mineral', thisPhaseInfo.name);
end

cs{1} = 'notIndexed';

% Alternative: you can use cif files
% cifFiles = {'SiC_4H.cif', 'SiC-6H_2.cif'};
% cifLocation = "D:\AstroEBSD_v2_standalone\phases\cifs";
% 
% for i=1:size(mapPhases.ActivePhasesList,1)
%     cs{i} = loadCIF(fullfile(cifLocation, cifFiles{i}));
% end
%% Read orientation data

multiphase = ReadTFSOrientationData(fullfile(mapLocation, 'results/multiphases.idx'));
multiphase.PhaseID = multiphase.PhaseID + 1; % this makes 0 = 'unindexed'.

% convert to radians for mtex
multiphase.Phi1 = multiphase.Phi1 * degree;
multiphase.PHI = multiphase.PHI * degree;
multiphase.Phi2 = multiphase.Phi2 * degree;

% can load other data if needed.

%% Plot IPF

% axis convention - confirmed by comparing with preview maps in xTalView
setMTEXpref('xAxisDirection','east'); 
setMTEXpref('zAxisDirection','outofplane');

% build the coordinate maps. To get the accurate micron bar, we multiply by
% the step size.
prop.x = double(multiphase.x) * mapStepSize;
prop.y = double(multiphase.y) * -mapStepSize;

ori = rotation('Euler', multiphase.Phi1, multiphase.PHI, multiphase.Phi2);

try
ebsd = EBSD(ori, multiphase.phaseId,cs,'Options',prop);
catch
    ebsd = EBSD(ori,  multiphase.PhaseID,cs,prop);
end

%% Plot IPF - multiphase

figure; 

for i = 1:size(mapPhases.ActivePhasesList,1)

plot(ebsd(cs{i+1}.mineral),colorKey.orientation2color(ebsd(cs{i+1}.mineral).orientations));
hold on;


end

hold off

% pot colorkey - note that this is different from the colorkey used in
% xTalView
for i = 1:size(mapPhases.ActivePhasesList,1)
ipfKey = ipfHSVKey(cs{i+1});
figure;
plot(ipfKey);
end

% alternatively you may want to try the following, to plot all 3 IPFs in
% the same figure object.
% colorKey = ipfHSVKey(cs{1+1});
% colorKey.inversePoleFigureDirection = xvector;
% color = colorKey.orientation2color(ebsd(cs{i+1}.mineral).orientations);
% plot(ebsd(cs{i+1}.mineral),color); title('IPF-X');
% nextAxis;
% colorKey.inversePoleFigureDirection = yvector;
% color = colorKey.orientation2color(ebsd(cs{i+1}.mineral).orientations);
% plot(ebsd(cs{i+1}.mineral),color); title('IPF-Y');
% nextAxis;
% colorKey.inversePoleFigureDirection = yvector;
% color = colorKey.orientation2color(ebsd(cs{i+1}.mineral).orientations);
% plot(ebsd(cs{i+1}.mineral),color); title('IPF-Z'); axis on;

%% Pole figures

odf1 = calcDensity(ebsd(cs{2}.mineral).orientations);
figure
plotPDF(odf1,[Miller(0,0,1,ebsd(cs{2}.mineral).CS) Miller(1,1,0,ebsd(cs{2}.mineral).CS)])
mtexColorbar('location','southoutside')
mtexColorMap blue2red

%% Other plots that you may like to have
% cShape = crystalShape.cube(ebsd.CS);

% [grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',2*degree);
% big_grains = grains(grains.grainSize > 2);

% figure;
% plot(ebsd,colorKey.orientation2color(ebsd.orientations));
% hold on;
% plot(grains.boundary,'linewidth',2);
% plot(ebsd('indexed'),0.4*cShape,'linewidth',2,'colored','micronbar','off');
% hold off;
% legend off;

