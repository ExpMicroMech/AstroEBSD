%% Example script of loading orientation data from Apreo

%%
InputUser.mtex_location='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';

% start mtex if needed
try EBSD;
catch
    run(fullfile(InputUser.mtex_location,"startup.m"));
end

%%
mapLocation = "D:\ApreoTruePixData\2025-04-04T10_14_42_Tianbi_Ni_redo\2025-04-04T10_20_08";
mapJson = "site.json";

EBSDMapInfo = ReadTFSJson(mapLocation, mapJson);

%% Extract some metadata

mapWidth = EBSDMapInfo.Info.Width;
mapHeight = EBSDMapInfo.Info.Height;
mapStepSize = EBSDMapInfo.Info.MapStepSize * 1e6; % um
mapkeV = EBSDMapInfo.Info.LandingEnergy / 1e3; % keV
mapBC = EBSDMapInfo.Info.BeamCurrent * 1e9; % nA
mapTH = EBSDMapInfo.Info.Threshold / 1e3; % keV

%% Load phase data
% using cif is better since mtex doesn't like json
cifFiles = {'SiC_4H.cif', 'SiC-6H_2.cif'};
cifLocation = "D:\AstroEBSD_v2_standalone\phases\cifs";

for i=1:size(mapPhases.ActivePhasesList,1)
    cs{i} = loadCIF(fullfile(cifLocation, cifFiles{i}));
end

% this doesn't look nice, consider using CIFs instead.
% mapPhases = ReadTFSJson(mapLocation, "phases/active.json");
% for i=1:size(mapPhases.ActivePhasesList,1)
%     thisPhaseJson = mapPhases.ActivePhasesList(i).File;
% 
%     thisPhaseInfo = ReadTFSJson(fullfile(mapLocation, 'phases'), thisPhaseJson);
% 
%     cs{i} = crystalSymmetry(thisPhaseInfo.laue_group_symbol);
%     cs{i}.alpha = thisPhaseInfo.alpha * degree;
%     cs{i}.beta= thisPhaseInfo.beta * degree;
%     cs{i}.gamma = thisPhaseInfo.gamma * degree;
%     cs{i}.mineral = thisPhaseInfo.name;
% end

%% Read orientation data

multiphase = ReadTFSOrientationData(fullfile(mapLocation, 'results/multiphases.idx'));
multiphase.PhaseID = multiphase.PhaseID + 1; % this makes 0 = 'unindexed'.
multiphase.Phi1 = multiphase.Phi1 * degree;
multiphase.PHI = multiphase.PHI * degree;
multiphase.Phi2 = multiphase.Phi2 * degree;

%% can load other data if needed.

%%
setMTEXpref('xAxisDirection','east'); 
setMTEXpref('zAxisDirection','intoplane');

% build the coordinate maps. To get the accurate micron bar, we multiply by
% the step size.
prop.x = double(multiphase.x) * mapStepSize;
prop.y = double(multiphase.y) * mapStepSize;

ori = rotation('Euler', multiphase.Phi1, multiphase.PHI, multiphase.Phi2);

try
ebsd = EBSD(ori, ones(size(ori)),{'notIndexed',cs{1}, cs{2}},'Options',prop);
catch
    ebsd = EBSD(ori, ones(size(ori)),{'notIndexed',cs{1}, cs{2}},prop);
end

% ebsd_rot = rotate(ebsd,10*degree);

% idk hor to figure out ipfKey for multiphases though
ipfKey = ipfHSVKey(cs{1});
figure;
plot(ipfKey);

colorKey = ipfHSVKey(cs{1});
colorKey.inversePoleFigureDirection = xvector;
color = colorKey.orientation2color(ebsd('indexed').orientations);

figure;
plot(ebsd,colorKey.orientation2color(ebsd.orientations));

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

