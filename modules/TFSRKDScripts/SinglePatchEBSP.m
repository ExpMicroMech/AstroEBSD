%% SinglePatchEBSP.m
% Example script of loading patterns from TFS xTalView projects which only
% have a single diffraction pattern.
% Used in parametric studies.

% N.B. For generating the reprojected pattern for XC, the AstroEBSD toolbox
% is needed.

% Tianbi Zhang and T. Ben Britton, April 2025

%% Directories
TFSProjectLocation = "D:\ApreoTruePixData\2025-03-31T11_47_23_Tianbi_test_5keV";
projectJson = fullfile(TFSProjectLocation, "project.json");

projectJsonText = fileread(projectJson);
projectInfo = jsondecode(projectJsonText);

rawPatFileName = "f_00000.pak";
procPatFileName = "row_00000.bin";

% genmerate a simulated pattern for XCC - requires AstroEBSD for
% reprojection
% Gen_Sim_Patterns;

%% Pre-allocate arrays for parameters

totalNumPats = size(projectInfo.Sites,1) -1;

patternSize = 256;
allPatterns = zeros(patternSize,patternSize,totalNumPats);
allProcPatterns = zeros(patternSize,patternSize,totalNumPats);
allProjects = projectInfo.Sites;

allExposure = zeros(totalNumPats,1);
allTh = zeros(totalNumPats,1);
allBC = zeros(totalNumPats,1);
allSumCount = zeros(totalNumPats,1);
allMEC = zeros(totalNumPats,1); 
% XCf = zeros(totalNumPats,1); 

for i = 2:size(projectInfo.Sites,1)
    thisProject = projectInfo.Sites(i);
    thisProjectName = thisProject.Name;
    thisProjectDir = thisProject.Dir;

    thisPatFile = fullfile(TFSProjectLocation, thisProjectDir, "raw", rawPatFileName);
    thisProcPatFile = fullfile(TFSProjectLocation, thisProjectDir, "patterns", procPatFileName);

    fileID = fopen(thisPatFile);
    rawPatVector = fread(fileID, 'int16'); % pay attention to this- might be different!!!
    fclose(fileID);

    fileID = fopen(thisProcPatFile);
    procPatVector = fread(fileID); % pay attention to this- might be different!!!
    fclose(fileID);
    
    rawPat = reshape(rawPatVector, [patternSize patternSize])';
    procPat = reshape(procPatVector, [patternSize patternSize])';
    allPatterns(:,:,i-1) = rawPat;
    allProcPatterns(:,:,i-1) = rot90(procPat,2);

    thisProjectJson = fullfile(TFSProjectLocation, thisProjectDir, "site.json");
    thisProjectJsonText = fileread(thisProjectJson);
    thisProjectInfo = jsondecode(thisProjectJsonText);

    allExposure(i-1) = (thisProjectInfo.Info.ExposureTime - 0.0003) * 1e3; % ms
    allTh(i-1) = thisProjectInfo.Info.Threshold / 1e3; % keV
    allBC(i-1) = thisProjectInfo.Info.BeamCurrent * 1e9; % nA
    allSumCount(i-1) = sum(rawPatVector);
    allMEC(i-1) = median(rawPatVector);

end

% fraction of incident electrons detected by EBSD
allDoseRatio = allSumCount ./ (allBC .* allExposure ./ (1.6 * 10^(-19)) / (10^(12)));
%% Write to image files

for selected = 1:totalNumPats
formattedTitle = sprintf("%.1fms_%.1fkeVTh_%.2fnA.tif", ...
    allExposure(selected), allTh(selected), allBC(selected));

formattedProcTitle = sprintf("%.1fms_%.1fkeVTh_%.2fnA_bgcor.png", ...
    allExposure(selected), allTh(selected), allBC(selected));

thisPattern = allPatterns(:,:,selected);
thisProcPattern = allProcPatterns(:,:,selected);

% figure;
% subplot(1,2,1);
% imagesc(thisPattern); 
% axis xy; axis image; axis off; 
% colormap('gray');
% subplot(1,2,2);
% imagesc(thisBgCorPattern); 
% axis xy; axis image; axis off; 
% colormap('gray');
% title(formattedTitle);

imgPath = fullfile(TFSProjectLocation, formattedTitle);
imgProcPath = fullfile(TFSProjectLocation, formattedProcTitle);

imwrite(normalizeto16bit(flipud(thisPattern)), imgPath);
imwrite(normalizeto16bit(flipud(thisProcPattern)), imgProcPath);
end

%% Assorted plots
figure;
scatter(allTh, allDoseRatio, 30, allBC, 'filled');
xlabel("Energy threshold (keV)");
ylabel("Median electron count");
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')

%
figure;
scatter(allTh(allBC>0.1), allMEC(allBC>0.1) , 30, allExposure(allBC>0.1), 'filled');
xlabel("Energy threshold (keV)");
ylabel("Median electron count");
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')


%% XCC

% templatePat = Pat_sim_eang1;
% 
% for i=1:totalNumPats
% 
% thisPat = allProcPatterns(:,:,i);
% 
% XCf(i) = max(normxcorr2(templatePat, thisPat), [], "all");
% end

%
figure;
imagesc(Pat_sim_eang1); axis xy; axis image; colormap('gray')
figure;
imagesc(allProcPatterns(:,:,1)); axis xy; axis image; colormap('gray')

%
cmaptemp = cool(4);
figure; hold on;
% scatter(allTh, XCf, 30, allBC, 'filled');
plot(allTh(allBC==0.4), XCf(allBC==0.4), 'Marker','o', ...
    'LineWidth',1,'Color',cmaptemp(1,:), 'MarkerFaceColor',cmaptemp(1,:), ...
    'MarkerSize',5,'DisplayName','0.4 nA');
plot(allTh(allBC==0.8), XCf(allBC==0.8), 'Marker','o', ...
    'LineWidth',1,'Color',cmaptemp(2,:), 'MarkerFaceColor',cmaptemp(2,:), ...
    'MarkerSize',5,'DisplayName','0.8 nA');
plot(allTh(allBC==1.6), XCf(allBC==1.6), 'Marker','o', ...
    'LineWidth',1,'Color',cmaptemp(3,:), 'MarkerFaceColor',cmaptemp(3,:), ...
    'MarkerSize',5,'DisplayName','1.6 nA');
plot(allTh(allBC==3.2), XCf(allBC==3.2), 'Marker','o', ...
    'LineWidth',1,'Color',cmaptemp(4,:), 'MarkerFaceColor',cmaptemp(4,:), ...
    'MarkerSize',5,'DisplayName','3.2 nA');

xlabel("Energy threshold (keV)");
ylabel("Normalized XCC")
grid on;
legend('Location','northoutside','NumColumns',4)
% a = colorbar;
% a.Label.String = 'Beam current (nA)';
% colormap('cool')

%
figure;
scatter(allTh, XCf, 30, allExposure, 'filled');
xlabel("Energy threshold (keV)");
ylabel("Normalized XCC")
grid on;
a = colorbar;
a.Label.String = 'Exposure time (ms)';
colormap('cool')


figure;
scatter(allMEC, XCf, 30, allTh, 'filled');
xlabel("Median electron count");
ylabel("Normalized XCC")
grid on;
a = colorbar;
a.Label.String = 'Energy threshold (keV)';
colormap('cool')

%
figure;
scatter(allExposure, XCf, 30, allMEC, 'filled'); hold on;
% scatter(allExposure, XCf2, 30, allMEC, 'filled'); hold on;
xlabel("Exposure time (ms)");
ylabel("Normalized XCC")
grid on;
a = colorbar;
a.Label.String = 'Median electron count';
colormap('cool')


figure;
scatter(allExposure, allMEC, 'filled'); hold on;
xlabel("Exposure time (ms)");
ylabel("Median electron count")
grid on;
colormap('cool')

%

figure;
scatter(allTh, allDoseRatio, 30, allBC, 'filled'); hold on;
xlabel("Energy threshold (keV)");
ylabel({'Fraction of' ' incident' 'electrons' 'detected' 'by EBSD'},"Rotation",0)
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')