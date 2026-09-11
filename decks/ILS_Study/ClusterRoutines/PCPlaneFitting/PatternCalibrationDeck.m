clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
% close all;
% % InputUser.MTEX_loc=MTEX_loc
% % InputUser.MTEX_loc='C:\Users\berners\mtex-5.7.0'
% try EBSD;
% catch
%     run(fullfile(InputUser.MTEX_loc,"startup.m"));
% end
% 
%
%% sourcing out stuff to a config .m file


run __INDEXCONFIGFILE__.m
exportpath = InputUser.ExportPath;
%optional :
% run plotEBSDDataIndexing.m
%%
% Settings_Cor

%%
%%
if isfield(MapInfo,'g_num')
    MapInfo.MapData.g_num=MapInfo.g_num;
end

%%

%%
InputUser.RTMBinningFactor = 8;

pat=bReadEBSP(MapInfo.EBSPData,pnum);

[patternIn,InputUser.PatternInfo] = EBSP_BGCor(pat,Settings_Cor);
InputUser.PatternInfoRTM = InputUser.PatternInfo;
Settings_Cor.gaussBlur = 0; %reduces noise a bit
Settings_Cor.SettingsCorRTM = Settings_Cor;
InputUser.SettingsCor = Settings_Cor;
InputUser.SettingsCorRTM.size =  double(floor(InputUser.SettingsCor.size/InputUser.RTMBinningFactor));
InputUser.PatternInfoRTM.size = double(floor(InputUser.SettingsCor.size/InputUser.RTMBinningFactor));
RTM.screensize =double(InputUser.PatternInfoRTM.size(1));
RTM.LPTsize = 2* double(RTM.screensize);

%%
[PhaseParameters,FFTParameters,EBSD_geometry] = prepareRefinementMultiPhase(InputUser,RTM)
%%
InputUser.SettingsCorRTM.cross = 2;

% figure
% imagesc(patternIn1)
%%
% profile on 
pnum = 1
% for i = 1:10
[G_Refined, PeakHight] = IndexPatternFFT(patternIn,FFTParameters,RTM,PhaseParameters,InputUser,pnum)
%% start the hough based refinement Input parameters in analogy to input deck single, but we can take some rom the loaded data


%% Test EDAX and Pattern alignment
% [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
% [EBSD_simulation ] = EBSP_Gnom( RTM_info,InputUser.PC); %
% [ EBSD_geometry ] = EBSP_Gnom(pattern_info,Settings_PCin.start);

%% Construct a grid to extract patterns for refinement
InputUser.gridPoints = gridPoints;
InputUser.kernelSize = kernelSize;

%%

InputUser.EBSD_geometry = EBSD_geometry;
%% Extract the patterns we want to refine

kernelSize = InputUser.kernelSize;
PatternSearchSize = pSerachSize;
gridRes = uint32((MicroscopeData.NCOLS -2*kernelSize-5)/ InputUser.gridPoints);
[XGrid,YGrid] = meshgrid(1 + kernelSize + PatternSearchSize:gridRes:MicroscopeData.NCOLS - kernelSize - PatternSearchSize, ....
                1+kernelSize+PatternSearchSize:gridRes:MicroscopeData.NROWS-kernelSize-PatternSearchSize);
for gridNumber = 1:size(XGrid(:))
    allX = XGrid(gridNumber)-PatternSearchSize:XGrid(gridNumber)+PatternSearchSize;
    allY = YGrid(gridNumber)-PatternSearchSize:YGrid(gridNumber)+PatternSearchSize;
    [XSearchGrid,YSearchGrid] = meshgrid(allX(:), allY(:));

    [pnum,~]=find(MapInfo.MapData.XBeam==XSearchGrid(:)' & MapInfo.MapData.YBeam==YSearchGrid(:)');
    [maxCI,ciIndex] = max(ebsd_map(pnum).ci.*ebsd_map(pnum).iq); %filter for patterns within a ragne that have maximum CI, from this we now, that the indexing is likeyl correct
    % CI is not the best metric though!!!
    % gridNumber
    maxCI;
    fprintf("CI Index: %i\n", ciIndex);
    fprintf("Xsearch grid at CIIndex %i\n", XSearchGrid(ciIndex));
    fprintf("X grid at grid nunmber: %i\n",XGrid(gridNumber));
    % disp(XSearchGrid(ciIndex))
    % disp( XGrid(gridNumber) )
    % disp(XGrid)
    XGrid(gridNumber) = XSearchGrid(ciIndex);
    YGrid(gridNumber) = YSearchGrid(ciIndex);
end

%
for i = 1:size(XGrid(:),1)
    ebsd_pat = kernelAveraging(XGrid(i),YGrid(i),kernelSize,MapInfo);
    XGrid(i),YGrid(i);
    patternsRefine(:,:,i) = ebsd_pat;%,Settings_Cor;
end
%%
Settings_Cor.cross=2
figure
imagesc(EBSP_BGCor(patternsRefine(:,:,i),Settings_Cor))
colormap gray
clim([-3,3])
axis equal
axis xy

%%
[patternNumber,~]=find(MapInfo.MapData.XBeam==XGrid(:)' & MapInfo.MapData.YBeam==YGrid(:)');
%%
% pat = bReadEBSP(MapInfo.EBSPData,105727);
% imagesc(pat)
% colormap gray;
% axis equal;
% axis xy;
%%
% pnum=1
% regroup in the future
InputUser.XGrid = XGrid;
InputUser.YGrid = YGrid;
InputUser.pnum = pnum;
if ~isfield(options,'gradients')
    options.gradients=0;
end
if options.gradients==1
    InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('__EBSDMAP___PCRefinement___JOBID__gradients.hdf5'))
else
    InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('__EBSDMAP___PCRefinement___JOBID__.hdf5'))
end
try
    makeNewH5(InputUser.h5name)
catch
 warning("File could not be created, mak sure, you have setup git properly on your machine")
end
%
try
    writeStructToH5(InputUser, InputUser.h5name,"/UserInfo")
catch
 warning("Could not write Input Info to h5")
end
try
    writeStructToH5(RTM, InputUser.h5name,"/UserInfo/RTMInfo")
catch
 warning("Could not write Input Info to h5")
end
try
    writeStructToH5(MicroscopeData, InputUser.h5name,"/MicroscopeData")
catch
 warning("Could not write Input Info to h5")
end

try
    writeStructToH5(MapInfo, InputUser.h5name,"/MapInfo")
catch
 warning("Could not write Input Info to h5")
end
%%

options.SettingsCor = Settings_Cor
% iteration options
% options.MaxIter = 20;
% options.StepTolerance = 1e-8;
% options.OptimalityTolerance = 1e-8;
% options.BGVar = 0;
% options.BGVarVal = 10;
% options.Diagnostics = 'off';

% for StepTolerance = [1e-4,1e-6,1e-8]
% for BGVar = [0,1]
% for patSize =514%[64,128,256,514]% [64,128,256,514]
%     options.StepTolerance = StepTolerance;
%     options.BGVar = BGVar;
% 
%     Settings_Cor.size = [patSize, patSize]
%     InputUser.ExportPath = fullfile(exportpath,sprintf("Size%i",patSize),sprintf("Tolerance%.2e_BgVar%i",StepTolerance,BGVar))
% end
% end
% end
%%
for phaseNum = 1: size(InputUser.Phases,2)
screen_int = PhaseParameters.screen_int(phaseNum);
for maxIter = [300] % 1200 would be "defualt for 6 varirable 200*numVar)
for patSize =[514]% [64,128,256,514]
    % if patSize >=256a
    %     SettingsCor.gaussBlur = 2; %reduces noise a bit
    % end
    Settings_Cor.size = [patSize, patSize]
    options.Display='final';
    radPerPix = atan(1/(2*InputUser.PC(3)))* (1/patSize);
    pcPerPix = (1/patSize);
    options.scalingFactor =scalingFactor;
    options.oriScaling = options.scalingFactor * 1/radPerPix;
    options.pcScaling = options.scalingFactor * 1/pcPerPix;
    % options.BGVar = BGVar;
    % options.StepTolerance = StepTolerance;
    options.SettingsCor = Settings_Cor
    options.MaxIter = maxIter
    InputUser.ExportPath = fullfile(exportpath,sprintf("/PatSize_%i/MaxIter_%i/%s",patSize,maxIter,InputUser.Phases{phaseNum}));
    mkdir(InputUser.ExportPath)
    % InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('PatternPlots.hdf5'))
    InputUser.hdf5DataSet = sprintf("/Patterns/PatSize_%i/MaxIter_%i/%s",patSize,maxIter,InputUser.Phases{phaseNum});
    tic    
    [PCBest,eangs_out,phbest,eangs] = fitPatternCentreMultiPoint(MapInfo,MicroscopeData,InputUser,Settings_Cor,screen_int,patternsRefine,patternNumber,options)
    disp("----------- Time Passed for fitting step -----------")
    toc
    disp("------------- continue plotting --------------------")
    disp("")
    disp("")
    save( fullfile(InputUser.ExportPath,sprintf('ResultsPatternFit_%i',patSize)),"PCBest","phbest","eangs_out","options","InputUser","Settings_Cor","MicroscopeData","MapInfo","patternsRefine","pnum","eangs","XGrid" ,"YGrid")
% %%
%     f = figure
%     subplot(3,2,1)
%     plt_lim = Settings_Cor.max_var_pc*double(Settings_Cor.size(1));
%     scatter(double(XGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(1)-PCBest(:,1))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     colorbar
%     xlabel('X , µm')
%     ylabel('\Delta PCX, px')
%     % ylim([0.45,0.55])
%     % ylim([-plt_lim,plt_lim])
%     title('PCX and x coord')
    
    
%     subplot(3,2,2)
     % figure
%     scatter(double(YGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(1)-PCBest(:,1))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     % ylim([-plt_lim,plt_lim])
%     xlabel('Y , µm')
%     ylabel('\Delta PCY, px')
%     title('PCy and y coord')
%     colorbar
%     %
%     % subplot(3,3,2)
%     % % figure
%     % scatter(eangs_out(:,2)/degree,(InputUser.PC(1)-PCBest(:,1))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     % % ylim([-plt_lim,plt_lim])
%     % xlabel('PHI , °')
%     % ylabel('\Delta PCX, px')
%     % title('PCx and x coord')
%     % colorbar
%     %
%     plt_lim = Settings_Cor.max_var_pc*double(Settings_Cor.size(1));
%     scatter(double(XGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(2)-PCBest(:,2))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     colorbar
%     xlabel('X , µm')
%     ylabel('\Delta PCY, px')
%     % ylim([0.45,0.55])
%     % ylim([-plt_lim,plt_lim])
%     title('PCy and x coord')
    
    
%     subplot(3,2,4)
%     % figure
%     scatter(double(YGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(2)-PCBest(:,2))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     % ylim([-plt_lim,plt_lim])
%     xlabel('Y , µm')
%     ylabel('\Delta PCY, px')
%     title('PCY and y coord')
%     colorbar
    
    
%     subplot(3,2,5)
%     % figure
%     scatter(double(XGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(3)-PCBest(:,3))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     % ylim([-plt_lim,plt_lim])
%     xlabel('X, µm')
%     ylabel('\Delta PCZ, px')
%     title('PCZ and x coord')
    
    
%     colorbar
    
    
%     subplot(3,2,6)
%     % figure
%     scatter(double(YGrid(:)*InputUser.calibration.XSTEP),(InputUser.PC(3)-PCBest(:,3))*double(Settings_Cor.size(1)),30,phbest,'filled')
%     % ylim([-plt_lim,plt_lim])/
%     ylabel('\Delta PCZ, px')
%     title('PCZ and y coord')
%     colorbar
%     f.Position = [100 100 1000 1000]
%     exportgraphics(gcf,fullfile(InputUser.ExportPath,sprintf('PatternFit_%i.png',patSize)))
%     % title('PCX and x coord')
% % %% plot as imagessac
% %     f = figure
% %     subplot(1,3,1)
% %     % imagesc(plt_lim)
% %     XPlot = double(XGrid(:)*InputUser.calibration.XSTEP);
% %     YPlot = double(YGrid(:)*InputUser.calibration.XSTEP);
% %     scatter(XPlot,YPlot,100,PCBest(:,1),'filled')
% %     colormap(cmap_div)
% %     clim([mean(PCBest(:,1))-0.02,mean(PCBest(:,1))+0.02])
% %     colorbar
% %     axis equal
% %     ylim([min(YPlot(:)),max(YPlot(:))])
% %     xlim([min(XPlot(:)),max(XPlot(:))])
% %     f.Position = [100, 100, 1400,400]
% %     % f = figure
% %     subplot(1,3,2)
% %     % imagesc(plt_lim)
% %     XPlot = double(XGrid(:)*InputUser.calibration.XSTEP);
% %     YPlot = double(YGrid(:)*InputUser.calibration.XSTEP);
% %     scatter(XPlot,YPlot,100,PCBest(:,2),'filled')
% %     colormap(cmap_div)
% %     clim([mean(PCBest(:,2))-0.02,mean(PCBest(:,2))+0.02])
% %     colorbar
% %     axis equal
% %     ylim([min(YPlot(:)),max(YPlot(:))])
% %     xlim([min(XPlot(:)),max(XPlot(:))])
% %     f.Position = [100, 100, 1400,400]
% %     exportgraphics(gcf,fullfile(InputUser.ExportPath,sprintf('PatternFit_%i.png',patSize)))
% %%
% XPlot = double(XGrid*InputUser.calibration.XSTEP);
% YPlot = double(YGrid*InputUser.calibration.XSTEP);
%  f =    figure
%  titles = ["\DeltaPCX (px)", "\DeltaPCY (px)", "\DeltaPCZ (px)"]   ;
%  for pcVar = 1:3
    
%     PCPixels = (InputUser.PC(pcVar)-PCBest(:,pcVar))*double(Settings_Cor.size(1))
%     pcColor = reshape(PCPixels,size(XGrid));
%     subplot(1,3,pcVar)
%     s =pcolor(XPlot,YPlot,pcColor);
%     colormap(cmap_div)
%     colorbar
    
%     clim([mean(PCPixels(:))-2*std(PCPixels(:)),mean(PCPixels(:))+2*std(PCPixels(:))])
%     axis equal
%     ylim([min(YPlot(:)),max(YPlot(:))])
%     xlim([min(XPlot(:)),max(XPlot(:))])
%     xlabel('x, µm')
%     ylabel('y,µm')
%     title(titles(pcVar))
%     s.FaceColor = 'interp';
%         s.Marker = "o";
%         % s.MarkerSize = "10"
%         % s.EdgeColor = 'interp'
%         s.MarkerFaceColor='flat';
% end
% f.Position = [100, 100, 1400,400]
% exportgraphics(gcf,fullfile(InputUser.ExportPath,sprintf('DeltaPC_pixels%i.png',patSize)))
% %%
% XPlot = double(XGrid*InputUser.calibration.XSTEP);
% YPlot = double(YGrid*InputUser.calibration.XSTEP);
%  f =    figure
%  titles = ["PCX", "PCY", "PCZ"]   
%  for pcVar = 1:4
%      ax= subplot(1,4,pcVar)
%      if pcVar<4
%         PCPixels = PCBest(:,pcVar)%)*double(Settings_Cor.size(1))
%          pcColor = reshape(PCPixels,size(XGrid));
 
%         % s.Mark
%         s =pcolor(XPlot,YPlot,pcColor);
%         s.FaceColor = 'interp';
%         s.Marker = "o";
%         % s.MarkerSize = "10"
%         % s.EdgeColor = 'interp'
%         s.MarkerFaceColor='flat';
%         colormap(cmap_div);
%         colorbar
%         % clim([mean(PCBest(:,pcVar))-0.02,mean(PCBest(:,pcVar))+0.02])
     
%         title(titles(pcVar))
%         clim([mean(PCPixels(:))-2*std(PCPixels(:)),mean(PCPixels(:))+2*std(PCPixels(:))])
%      else
         
%        pcColor = reshape(phbest,size(XGrid));
        
%         s=pcolor(XPlot,YPlot,pcColor);
%         s.FaceColor = 'interp';
%         s.Marker = "o";
%         % s.MarkerSize = "10"
%         % s.EdgeColor = 'interp'
%         s.MarkerFaceColor='flat';
%         title('NCC')
%         colormap(ax,flipud(summer))
%         colorbar
%         clim([mean(pcColor(:))-1*std(pcColor(:)),mean(pcColor(:))+1*std(pcColor(:))])
%      end
%         axis equal

%     ylim([min(YPlot(:)),max(YPlot(:))])
%     xlim([min(XPlot(:)),max(XPlot(:))])
%     xlabel('x, µm')
%     if pcVar ==1
%         ylabel('y,µm')
%     end
    

% end
% f.Position = [100, 100, 1700,400];
% exportgraphics(gcf,fullfile(InputUser.ExportPath,sprintf('PCVariation_Relative%i.png',patSize)))
% %%
% XPlot = double(XGrid*InputUser.calibration.XSTEP);
% YPlot = double(YGrid*InputUser.calibration.XSTEP);
%  f =    figure;
%  titles = ["\Deltaphi1 (°)", "\DeltaPhi (°)", "\Deltaphi2 (°)"]   ;
%  eangs_inConv = [eangs.phi1(:),eangs.Phi(:),eangs.phi2(:)]
%  eangs_outConv = [eangs_out.phi1(:),eangs_out.Phi(:),eangs_out.phi2(:)]
%  for pcVar = 1:3
    
%     PCPixels = (eangs_inConv(:,pcVar)-eangs_outConv(:,pcVar))/degree;
%     pcColor = reshape(PCPixels,size(XGrid));
%     subplot(1,3,pcVar)
%     s=pcolor(XPlot,YPlot,pcColor);
%     s.FaceColor = 'interp';
%         s.Marker = "o";
%         % s.MarkerSize = "10"
%         % s.EdgeColor = 'interp'
%         s.MarkerFaceColor='flat';
%     colormap(cmap_div)
%     colorbar
    
%     clim([mean(PCPixels(:))-2*std(PCPixels(:)),mean(PCPixels(:))+2*std(PCPixels(:))])
%     axis equal
%     ylim([min(YPlot(:)),max(YPlot(:))])
%     xlim([min(XPlot(:)),max(XPlot(:))])
%     xlabel('x, µm')
%     ylabel('y,µm')
%     title(titles(pcVar))

% end
% f.Position = [100, 100, 1400,400]
% exportgraphics(gcf,fullfile(InputUser.ExportPath,sprintf('EangVariation%i.png',patSize)))

end
end
end
% end