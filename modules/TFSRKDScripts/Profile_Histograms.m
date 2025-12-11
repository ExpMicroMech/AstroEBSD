pname1 = "D:\RKD_Data\multikeV\2";
% pname1 = "D:\RKD_Data\2023-11-16T12_36_34_WSe2_orig\2023-11-16T13_36_03";
pname2 = "results";
fname1 = "site.json";
fname2 = "counts.bin";
fname3 = "pq.bin";

%%
TFSJson = ReadTFSJson(pname1, fname1);

% get the MSA results
% load("Demo_Tianbi_16bin_PCA_.mat"); 

z2 = find(VMOutput.PCA_VM_num == 2);
z4 = find(VMOutput.PCA_VM_num == 4);
z10 = find(VMOutput.PCA_VM_num == 10);

allz = [z2; z4; z10];

%%
width = TFSJson.Info.Width;
height = TFSJson.Info.Height;

step = TFSJson.Info.MapStepSize * 1e6; % um

keV = TFSJson.Info.HighVoltage;
beamCurrent = TFSJson.Info.BeamCurrent; % A
exposureTime = TFSJson.Info.ExposureTime; % s

patternSize = TFSJson.Info.Geometry.PatternHeight * TFSJson.Info.Geometry.PatternWidth;

scaleFactor = beamCurrent /(1.6 * 10^(-19)) * exposureTime;

%%

mecmapnonorm = ReadTFSBinaryMapData(fullfile(pname1, pname2), fname2, TFSJson); % raw MEC
mecmap = ReadTFSBinaryMapData(fullfile(pname1, pname2), fname2, TFSJson) * patternSize ./ scaleFactor; % normalized total dose
pqmap = ReadTFSBinaryMapData(fullfile(pname1, pname2), fname3, TFSJson); % PQ

% get mean and std region by region
selectedZ = z2;
% mean(mecmapnonorm(selectedZ), "all")
% std(mecmapnonorm(selectedZ), "all")
%% MEC/PQ maps with line profile overlays

figure;
imagesc(mecmap); axis image; colorbar; axis off;
% clim([0 0.05]);
set(gca, "FontSize", 18);
figureName = strcat(num2str(keV/1e3), '.png'); 
hold on;
xx1 = [102 66];
yy1 = [135 137]; diff1 = [xx1(2) - xx1(1) yy1(2) - yy1(1)];
xx2 = [48 33];
yy2 = [56 55]; diff2 = [xx2(2) - xx2(1) yy2(2) - yy2(1)];
plot(xx1, yy1, 'LineWidth',2, 'Color', [1 1 1]);
plot(xx2, yy2, 'LineWidth',2, 'Color', [1 1 1]);

%%
aa = improfile(pqmap, xx1, yy1);

numPixel = size(aa,1);
distance = meshgrid(0:numPixel-1, 0:numPixel-1) ; % um
distance = distance(1,:) * step; % um

figure;
plot(distance, aa, 'LineWidth',2);
grid on;
xlabel("Distance (\mum)");
% ylabel("Normalized MEC");
% ylabel("Median Electron Count");
ylabel("Pattern Quality");
xlim([0 20]);

bb = improfile(pqmap, xx2, yy2);

numPixel = size(bb,1);
distance = meshgrid(0:numPixel-1, 0:numPixel-1) ; % um
distance = distance(1,:) * step; % um
figure;
plot(distance, bb, 'LineWidth',2);
grid on;
xlabel("Distance (\mum)");
% ylabel("Normalized MEC");
% ylabel("Median Electron Count");
ylabel("Pattern Quality");
xlim([0 8]);

%% Histograms of MEC and PQ for each region
cmap = [
228,26,28;
55,126,184;
77,175,74;
152,78,163;
255,127,0;
] ./ 255;

figure; 
subplot(2,1,1)
histogram(mecmapnonorm, 100, 'EdgeColor','none', 'FaceColor',cmap(4,:),'FaceAlpha', 0.5, 'DisplayName','Overall');hold on;
histogram(mecmapnonorm(allz), 100, 'EdgeColor','none', 'FaceColor',cmap(5,:),'FaceAlpha', 0.5, 'DisplayName','Flakes');hold on;
histogram(mecmapnonorm(z10), 100, 'EdgeColor','none', 'FaceColor',cmap(1,:), 'DisplayName','F1');
histogram(mecmapnonorm(z2), 100, 'EdgeColor','none', 'FaceColor',cmap(2,:), 'DisplayName','F2'); 
histogram(mecmapnonorm(z4), 100, 'EdgeColor','none', 'FaceColor',cmap(3,:), 'DisplayName','F3');
xlabel("Median Electron Count");
ylabel("Counts");
legend('Location','northeast')
subplot(2,1,2)
histogram(pqmap, 100, 'EdgeColor','none', 'FaceColor',cmap(4,:),'FaceAlpha', 0.5, 'DisplayName','Overall');hold on;
histogram(pqmap(allz), 100, 'EdgeColor','none', 'FaceColor',cmap(5,:),'FaceAlpha', 0.5, 'DisplayName','Flakes');hold on;
histogram(pqmap(z10), 100, 'EdgeColor','none', 'FaceColor',cmap(1,:), 'DisplayName','F1');
histogram(pqmap(z2), 100, 'EdgeColor','none', 'FaceColor',cmap(2,:), 'DisplayName','F2'); 
histogram(pqmap(z4), 100, 'EdgeColor','none', 'FaceColor',cmap(3,:), 'DisplayName','F3');
xlabel("Pattern Quality");
ylabel("Counts");


