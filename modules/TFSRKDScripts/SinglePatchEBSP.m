% TFSProjectLocation = "D:\ApreoTruePixData\2025-04-02T12_44_08_Tianbi_Si_20keV_redo_set";
TFSProjectLocation = "D:\ApreoTruePixData\2025-03-31T11_47_23_Tianbi_test_5keV";
projectJson = fullfile(TFSProjectLocation, "project.json");

projectJsonText = fileread(projectJson);
projectInfo = jsondecode(projectJsonText);

rawPatFileName = "f_00000.pak";
procPatFileName = "row_00000.bin";

Gen_Sim_Patterns;
%%
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
XCf = zeros(totalNumPats,1); 

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

allDoseRatio = allSumCount ./ (allBC .* allExposure ./ (1.6 * 10^(-19)) / (10^(12)));
%%

for selected = 1:totalNumPats
formattedTitle = sprintf("%.1fms_%.1fkeVTh_%.2fnA.tif", ...
    allExposure(selected), allTh(selected), allBC(selected));

formattedProcTitle = sprintf("%.1fms_%.1fkeVTh_%.2fnA_bgcor.png", ...
    allExposure(selected), allTh(selected), allBC(selected));

thisPattern = allPatterns(:,:,selected);
thisProcPattern = allProcPatterns(:,:,selected);
thisMask = imgaussfilt(thisPattern,15);
thisBgCorPattern = thisPattern ./ thisMask;

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

%%

% AA = [allExposure allTh allBC allSumCount];
% csvFileName = fullfile(TFSProjectLocation, 'Conditions.csv');
% csvwrite(csvFileName, AA);

% cmap = [166,206,227;
% 31,120,180;
% 178,223,138;
% 51,160,44;
% 251,154,153;
% 227,26,28;
% 253,191,111;
% 255,127,0;
% 202,178,214;
% 106,61,154;
% 255,255,153;
% 177,89,40];
%%
figure;
scatter(allTh, allDoseRatio, 30, allBC, 'filled');
xlabel("Energy threshold (keV)");
ylabel("Median electron count");
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')
%%
figure;
scatter(allTh(allBC>0.1), allMEC(allBC>0.1) , 30, allExposure(allBC>0.1), 'filled');
xlabel("Energy threshold (keV)");
ylabel("Median electron count");
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')

%%
% figure;
% scatter(allMEC, allSumCount / 65536);
% xlabel("Median electron count");
% ylabel("Average electron count"); grid on;

%%
% maxPatternInd = find(allMEC == max(allMEC(:)));

templatePat = Pat_sim_eang1;

for i=1:totalNumPats
% thisPatNameAll = patList(i,:);
% thisPatName = strip(thisPatNameAll, 'right');
% exposure(i) = str2double(thisPatName(1:end-5))

thisPat = allProcPatterns(:,:,i);

XCf(i) = max(normxcorr2(templatePat, thisPat), [], "all");
end
%%
% templatePat = allProcPatterns(:,:,totalNumPats);
% 
% for i=1:totalNumPats
% % thisPatNameAll = patList(i,:);
% % thisPatName = strip(thisPatNameAll, 'right');
% % exposure(i) = str2double(thisPatName(1:end-5))
% 
% thisPat = allProcPatterns(:,:,i);
% 
% XCf2(i) = max(normxcorr2(templatePat, thisPat), [], "all");
% end
%% 
figure;
imagesc(Pat_sim_eang1); axis xy; axis image; colormap('gray')
figure;
imagesc(allProcPatterns(:,:,1)); axis xy; axis image; colormap('gray')
%%
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

%%
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

%%
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

%%

figure;
scatter(allTh, allDoseRatio, 30, allBC, 'filled'); hold on;
% scatter(allExposure, XCf2, 30, allMEC, 'filled'); hold on;
xlabel("Energy threshold (keV)");
ylabel({'Fraction of' ' incident' 'electrons' 'detected' 'by EBSD'},"Rotation",0)
grid on;
a = colorbar;
a.Label.String = 'Beam current (nA)';
colormap('cool')