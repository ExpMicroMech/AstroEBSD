%% L. Berners (Beginning of 2025), This deck uses a combination of the fourier based 
% refined template matching and the real space based optimisation routine.
% You need to specify a phase and a rough estimate of the pattern centre,
% and the deck will try to refine the PC and OR at the same time.
%
%
%
% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;

InputUser.Phase_Folder=Phase_Folder
InputUser.Astro_loc=Astro_loc
%% color coding
cmap_div =  returnCmapDiv()

InputData.cset=[1,0,0; % Colors of the Bands in rgb (1 =255)
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0];
InputData.cset=[InputData.cset;InputData.cset*0.5];
structConfig.cset = InputData.cset
%% File names


% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\YNi2MaDinah\YNi2_NS_D_repolish\Area1\_EBSD\vt0,1\Tile3\RotPCA_n12_w1,00e+00'
% InputUser.PatternFolder ="I:\001Clara\FeTaAl_Laves_Christina\Ta4-Fe2-Al6_repolish\Area2\_EBSD\vt0,1\Tile1\RotPCA_n14_w1,00e+00"
% InputUser.PatternFolder = 'E:\sciebo\IMM\ElectronHolography\EnergySeries(Holography)\BGCorrecction'
InputUser.PatternFolder = "J:\Users\Berners_Lukas\clara\C15_2_comparisonExperiment\Area1\_EBSD\vt0,1\Tile1\patternsForMatching"
InputUser.PatternList = ["_Pattern_1.png","_Pattern_2.png","_Pattern_3.png","_Pattern_4.png","_Pattern_5.png"]%,"_Pattern_2.png"]

% InputUser.ExportPath = 'E:\sciebo\IMM\EBSD\YNi2\YNi_NS_D\Area1\_EBSD_noHotpixCorr\vt0,1\Tile1\RotPCA_n17_w1,00e+00\PatternsForMatching'
InputUser.ExportPath = "J:\Users\Berners_Lukas\clara\C15_2_comparisonExperiment\Area1\_EBSD\vt0,1\Tile1\patternsForMatching"

InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('SinglePatRefinement.hdf5'))

InputUser.PCA_Mode = ""
InputUser.inputEulerMode = ""
%% Read Phases

% PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/Y*Ni*20kV*.pha'))'
% <<<<<<< HEAD
% 
% =======
% % 
% >>>>>>> b1a126541c31bace566ae213cb0959e25929d7f1
% InputUser.Phase_Input = cell(size(PhasesToRead));
% for i=1:numel(PhasesToRead)
%     InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
% end
% InputUser.Phase_Input = {'Y95Ni2_19.0kV_dmin0.020'}%'Nb6Co7_em','Nb518Ni482_joubert_20_2025_07_25__09_07_26_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.800_npx-1000'}
% subFolder = 'NbNiLPR_51Nb49Ni_IL'
% subFolder = "NbCo3aILGrid_LPR_High_2026_01_Replot"
% subFolder = "YNi_29.3kV"
% subFolder = "YNi2_20kVVoltageVar"
% subFolder = "YNi2_dbwCompare"
subFolder = "C15Al2Compare"
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
InputUser.Phase_Input = cell(size(PhasesToRead)-1);
for i=1:numel(PhasesToRead) % -1 for debugging
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end

% cutoff = 0.2       

% extractions  = struct();


for i = 1:numel(InputUser.Phase_Input)
    InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % load ebsd map
    rtmConfig = defaultRTMConfig()
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end


hList = {{0 0 4}',{4 0 4}', {4 4 4}', {2 2 6}', {2 2 4}',{6 6 10}'}


hType = 'HKL';
hnames = strings(numel(hList),numel(InputUser.Phase_Input))
InputUser.mPlotNames = hList
InputUser.htype = "HKL";

%% Get info for structure Parameters
FHKLoptions.Voltage = 19.5*10^3; % little bit below primary beam in kV
structConfig = defaultStructureFactorConfig(FHKLoptions);     
structConfig = defaultStructureFactorConfig(FHKLoptions);
structConfig.cutoff=0.2;
structConfig.minD = 0.5;
structConfig.maxInd = 10
%% Set up template matching stuff
% rtmConfig for coarse (Dictionary style pattern matching)
rtmConfig = defaultRTMConfig()
rtmConfig.generalOptions.AlignEdax =1 
rtmConfig.Settings_PCin.start = [0.529201, 1-0.683268, 0.693859]
InputUser.PC = rtmConfig.Settings_PCin.start 

%% Correction settings:
% rmtConfigRefine for FFT based and real space matching

rtmOptions.refine=1;
rtmOptions.clarity=1; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions)
rtmConfigRefine.generalOptions.AlignEdax =1 
rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start

Settings_Cor_plot= rtmConfigRefine.Settings_Cor,
%
Settings_Cor_plot.size=514
Settings_Cor_plot.cross=2
Settings_Cor_plot.quadrants=0
Settings_Cor_plot.gaussian_range = 1
%% Experiment Setup 
% For the map data stuff we will need to out source this.
MicroscopeData.SampleTilt = 70
MicroscopeData.CameraTilt = 1
% MicroscopeData.SampleTilt=InputUser.SampleTilt
% MicroscopeData.CameraTilt=InputUser.CameraTilt
% MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180;
MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180;
% MicroscopeData.TotalTilt=0
% Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);







%% Start h5 writing


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

MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfig)
%%
setupRTMSinglePats(rtmConfig,InputUser)%,MapInfo,MicroscopeData)


%%
%%
tic
libraryMatchPatternsCoarse(InputUser,rtmConfig,MultiPhaseStruct,1)
toc

%%

try
    writeStructToH5(rtmConfig, InputUser.h5name,"/rtmConfigRefine")
catch
 warning("Could not write rtmConfig to h5")
end




%%
InputUser.MicroscopeData = MicroscopeData
tic
PatternRefineFFT(InputUser,rtmConfigRefine,MultiPhaseStruct,1) 
toc
%%
numPatterns = numel(InputUser.PatternList)

optimisationOptions.Display = 'iter'
InputDataset = 'RTMRefined'
rtmConfigRefine.options = struct()
rtmConfigRefine.options.initPatSize = [514 514]
%
optimisationOptions.MaxIter = 200
realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct,1,optimisationOptions)
%%
InputUser.plotHoughAfterRefinement = 1;
InputUser.fitPC = 0
% InputUser.CSList = 
rtmConfigRefine.Settings_Cor = Settings_Cor_plot
plotRefinementResults(InputUser,rtmConfigRefine,structConfig,MultiPhaseStruct,1)
    %%
    InputUser.averagedOrisAfterRefinement = 0

    InputUser.h5nameRepro = sprintf("%s_repro",InputUser.h5name);
    % %% Plot reprojections
    % % Inpu
    % % tileNum =1
    plotReprojectionsFromH5(InputUser,rtmConfigRefine,MultiPhaseStruct,1)
    % %%
    % %% WARNING: DO NOT RUN. All profiles will be plotted NEED debugging
    % try
    %%
    InputUser.numHKLPLot=10;
    SampleBandProfilesFromReprojections(InputUser,structConfig)

%%
% PatternInfo_new.screensize = 514
% for i = 1
%     % for i = 1
%         %%
%         pc_best
%         [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,pc_best);
%         % [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,PC_best);
%         [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,g_refined,MultiPhaseStruct(phaseNumber).screen_int);
%         % [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_Refined,screen_int);
%         Pat_sim=EBSP_BGCor(Pat_sim,Settings_Cor_plot);
%         % [ Pat_sim_pseudo ] = EBSP_gen(EBSD_geom_plot,RTM.Rz(60*degree)*g_refined,screen_int);
%         % Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));Settings_Cor_plot
%           Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
%         % figure
%         Pat_exp_cor=EBSP_BGCor(Pat_exp,Settings_Cor_plot);
%         % Pat_exp=Pat_exp-mean(Pat_exp(:));
%         % Pat_exp=Pat_exp/std(Pat_exp(:));
%         %%
%         % imwrite(flipud(rescale(Pat_sim)),fullfile(export_path,['refined_pattern.tif']))
%         structConfig.Phase_Folder = Phase_Folder
%         [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
% 
% %%
% 
%         [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff)
%         m_plot = m_plot(1:10)
%         if m_plot.CS.numSym>=24 %cubic :D
%             I = I(1:5)
%             m_plot = m_plot(1:5)
%         else
%             I = I(1:10)
%         end
%         InputData.IPlot = I;
%         %%Pat_exp
%         InputUser.Pattern = Pat_exp;
%         InputUser.geometry = EBSD_geom_plot
%         InputUser.gmatrix_refined = conv_EA_to_G(eangs_out) 
% 
% %%
%         % [gridded_interpolant, InputUser] = PatternToInterpolant(InputUser,'ApplyMask',InputUser.PatternMask,'bgCor');
%         InputUser.pattern_info = struct();
%         InputUser.pattern_info.size = size( InputUser.Pattern );
%         InputUser.PatternMask = ones(size(InputUser.Pattern));
%         InputUser.PC_in = pc_best;
%         colorbar()
%         %%
%         try
%         InputUser= rmfield(InputUser,'PatternMask')
%         end
%         InputUser.cs = MultiPhaseStruct(phaseNumber).cs_phase;
%         InputUserSim = InputUser;
%         InputUserSim.Pattern = Pat_sim;
%         [gridded_interpolant, InputUser] = PatternToInterpolant(InputUser,'gaussian',1.5,'bgCor');
%         InputUserSim.Pattern = Pat_sim;
% 
% %
%         [gridded_interpolantSim, InputUserSim] = PatternToInterpolant(InputUserSim,'gaussian',1.5,'bgCor');
%       %
%         figure
%         imagesc(InputUser.Pattern_masked)
%         axis xy
%         axis equal
%         clim([-3.3,3.3])
%         colormap gray
% %% 
%         figure
%         imagesc(InputUserSim.Pattern_masked)
%         colorbar
%         clim([-3.3,3.3])
%            axis xy
%         axis equal
%         colormap gray
%        %%
%         img_out=reprojections_discrete(InputUser,gridded_interpolant,'img_size',2000,'fundamental');
%         imwrite(img_out,fullfile(export_path,[sprintf('resampledPattern_exp%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
%         tic
%         img_outSim=reprojections_discrete(InputUserSim,gridded_interpolantSim,'img_size',2000,'fundamental');
%         imwrite(img_outSim,fullfile(export_path,[sprintf('resampledPattern_Sim%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
% toc
%         %%
%         img_outSim = img_outSim -mean(img_outSim(:));
%         img_outSim = img_outSim/std(img_outSim(:));
%         img_out = img_out -mean(img_out(:));
%         img_out = img_out/std(img_out(:));
%         figure
%         imagesc(img_out)
%         axis('equal','xy')
%         colormap gray
%         clim([-3 3])
%         % %%
%         % img_out = ( img_out-mean(img_out(:)) );
%         % img_out = img_out / std( img_out(:));
%         %%
%         f = figure
%         f.Position = [100 100 1100 1100]
%         f.Units = 'pixels'
%         imagesc(img_outSim-img_out)
%         axis equal
%         axis xy
%         axis off
%         colormap(cmap_div)
%         %
%         % figure
%         % imagesc(diffImg)
%         exportgraphics(gcf,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',600)
%         % diffImg = img_outSim-img_out;%,cmap_div
%         % max(diffImg(:))
%         % min(diffImg(:))
%         %  figure
%         % imagesc(diffImg)
%         % % exportgraphics(gcf,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',600)
%         % rescale(diffImg)
%         % %%
%         % diffImg =rescale(diffImg)
%         % %%
%         % imwrite(diffImg,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]))
%         %%
%         % figure
%         % imshow(diffImg)
%         %% Plot the refinement results
%      f=   figure
%         % f.Position = [100 100 1500 600]
%         f.Units = 'pixels'
%         % figure
%         subplot(2,3,1)
%         imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim);
%         axis equal;
%         axis xy;
%         axis off
%         title('Simulated Pattern')
%         clim([-3.5,3.5])
%         ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         % figure
%         % ax5=subplot(2,3,4)
%         text(0,-0.07, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_out(1)/degree,eangs_out(2)/degree,eangs_out(3)/degree),'units','normalized','fontsize',14)
%         % text(0, -0.14, sprintf('LPR= %.4f , %.4f, %.4f', LPR(1),LPR(2),LPR(3)),'units','normalized','fontsize',14)
%         text(0, -0.21, sprintf('PC= %.4f , %.4f, %.4f', pc_best(1),pc_best(2),pc_best(3)),'units','normalized','fontsize',14)
%         text(0, -0.28, sprintf('NCC= %.4f', phbest),'units','normalized','fontsize',14)
%         % imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_sim_pseudo)
% 
%         % axis equal;
%         % axis xy;
%         % colormap(ax2,'gray');
%         % ylim([EBSD_geometry.y_gn_min,EBSD_geometry.y_gn_max])
%         % % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         % axis off
%         % title('Pseudosymmetry')
%         % clim([-3.5,3.5])
%         % m = Miller({1,1,-2,0},cs_phase)
%         % m=Miller({0,0,0,1},{2,0,-2,0},{2,2,-4,0},cs_phase)
%         % for j=1:size(ori,2)
%                 % r=ori(j).matrix;
%                 % euler=ori(j).Euler();
%                 % r = 
%         pattern =Pat_exp_cor;
%         % fname=fullfile(export_path, sprintf( '%s_gnomonic_%.1f_%.1f_%.1f_%i.png' ,InputUser.Phase_Input{i}, euler(1),euler(2),euler(3),j));
%         % imwrite(rescale(flipud(pattern)),fname)
%         % fname = 'I:\001Clara\EnergyFiltering\Si\RTM\bandsOnPattern.png'
%         % title(sprintf('Rz = %.0f °, Rx =%.0f°, Rz2 = %.0f°',angles(j),angles(i),0))
%         % label_bands(pattern,EBSD_geom_plot,m,inv(rotation.byMatrix(g_refined)),InputData)
%         % fname2=fullfile(export_path, sprintf( '%s_gnomonic_%.1f_%.1f_%.1f_%i_with_bands.png' ,InputUser.Phase_Input{i}, euler(1),euler(2),euler(3),j));
%         % exportgraphics(gcf,fname2,'resolution',300)
% 
%                 ax2=subplot(2,3,4)
%         structConfig.IPlot = InputData.IPlot
%         imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim-0.2)
%         label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),structConfig)
%         % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
%         axis equal;
%         axis xy;
%         colormap(ax2,'gray');
%         ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         axis off
%         title('Simulated. Pattern','fontsize',14)
%         clim([-3.5,3.5])
%         % axis(app.UIAxes, 'tight')
%         % set(gcf,'Resize','off')
% 
% 
% 
% 
%         ax2=subplot(2,3,2)
% 
%         imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
%         % label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
%         axis equal;
%         axis xy;
%         colormap(ax2,'gray');
%         ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         axis off
%         title('Exp. Pattern','fontsize',14)
%         clim([-3.5,3.5])
%         % axis(app.UIAxes, 'tight')
%         % set(gcf,'Resize','off')
% 
% 
%         ax2=subplot(2,3,5)
% 
%         imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
%         label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),structConfig)
%         % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
%         axis equal;
%         axis xy;
%         colormap(ax2,'gray');
%         ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         axis off
%         title('Exp. Pattern','fontsize',14)
%         clim([-3.5,3.5])
%         % axis(app.UIAxes, 'tight')
%         % set(gcf,'Resize','off')
% 
% 
%         ax=subplot(2,3,3)%subplot(1,3,3)
%         % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
%         imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
%         axis equal;
%         axis xy;
%         colormap gray;
%         ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         clim([-3,3])
%         % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         axis off
%         title('Difference')
%         colormap(ax,cmap_div)
%         colorbar('position',[0.91,0.63,0.01,0.25])        
% 
%         % set(gcf, 'Position' , [100 000 1500 900])
%         %%
%         % ax=subplot(2,3,6)%subplot(1,3,3)
%         % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,normalize(Pat_exp_cor)./normalize(Pat_sim))
%         % axis equal;
%         % axis xy;
%         % colormap gray;
%         % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
%         % clim([0.01,10])
%         % % caxis(log([0,10]));
%         % % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
%         % axis off
%         % title('Division')
%         % colormap(ax,cmap_div)
%         % colorbar('position',[0.91,0.13,0.01,0.25])        
%         % % set(ax6,'colorscale','log')
%         set(gcf, 'Position' , [100 000 1500 900])
%         f.Units = 'pixels'
%         exportgraphics(gcf, fullfile(export_path,[sprintf('refined_pattern_Pattern%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
%         % error('end')
%         close all
%         %%
%         numBins = 3
% binSize = 100/numBins
% bins = prctile(Pat_sim(:),[0:binSize:100])
% binsExp = prctile(Pat_exp_cor(:),[0:binSize:100])
% % f=figure
% f= figure
% f.Position=[20,20,1900,1200]
% h = tiledlayout(3,numBins, 'TileSpacing', 'none', 'Padding', 'none');
% for i = 1:numBins
% % subplot(5,3,i*3-2)
% greyAx = nexttile(i)
% % if bins(i)<0
%     % setValue = 
% plotImage = (imgaussfilt(Pat_sim,1.5));
% plotImage(plotImage<bins(i))=0;
% plotImage(plotImage>bins(i+1))=0;
% imagesc(plotImage)
% clim([bins(i),bins(i+1)])
% axis equal
% axis xy
% axis off
% colormap(greyAx,'gray')
% 
% % subplot(5,3,i*3-1)
% greyAx= nexttile(i+numBins)
% plotImageExp = (imgaussfilt(Pat_exp_cor,1.5));
% plotImageExp = imgaussfilt(plotImageExp);
% plotImageExp(plotImageExp<binsExp(i))=0;
% plotImageExp(plotImageExp>binsExp(i+1))=0;
% imagesc(plotImageExp)
% clim([binsExp(i),binsExp(i+1)])
% axis equal
% axis off
% axis xy
% colormap(greyAx,'gray')
% ax = nexttile(i+numBins*2)
% % plotImageExp = (imgaussfilt(Pat_exp_cor,1.5));
% % plotImageExp = imgaussfilt(plotImageExp);
% % plotImageExp(plotImageExp<bins(i))=-1000;
% % plotImageExp(plotImageExp>bins(i+1))=-1000;
% % fusedImage = imfuse(plotImageExp,plotImage);
% imagesc(plotImageExp-plotImage)
% % imagesc(fusedImage)
% clim([-3 3])
% colormap(ax,cmap_div)
% colorbar
% axis equal
% axis off
% axis xy
% end
% f.Position=[20,20,1900,1200]
% exportgraphics(gcf, fullfile(export_path,[sprintf('refined_pattern_Pattern_Thresholded%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
% close all
%     end
% % end
% %%
% %%% end
% %%
% 
% 
% % figure
% % surf(Pat_sim_Nb6Ni7-Pat_sim_Nb7Ni6,edgecolor='None')
% % 
% % colormap(cmap_div)
% % clim([-1.5,1.5])
% % view([20 80])
% % title('difference_simulation')
% % colormap(cmap_div)
% % 
% % %%
% % figure
% % blurr=imgaussfilt(Pat_sim_Nb6Ni7-Pat_sim_Nb7Ni6,5);
% % surf(blurr,edgecolor='None')
% % 
% % colormap(cmap_div)
% % clim([-1.5,1.5])
% % view([20 80])
% % title('difference_simulation')
% % colormap(cmap_div)
% % 
% % %%
% % figure
% % 
% % surf(Pat_sim_Nb6Ni7,edgecolor='None')
% % colormap(cmap_div)
% % clim([-3,3])
% % view([20 80])
% % title('Nb6Ni7')
% % axis off
% % %%
% % figure
% % surf(Pat_sim_Nb7Ni6,edgecolor='None')
% % colormap(cmap_div)
% % clim([-3,3])
% % view([20 80])
% % title('Nb7Ni6')
% % axis off
% % %%
% % figure
% % surf(Pat_exp,edgecolor='None')
% % colormap(cmap_div)
% % clim([-3,3])
% % view([20 80])
% % axis off
% % %%
% % figure
% % blurr=imgaussfilt(Pat_sim_Nb7Ni6-Pat_exp,5);
% % surf(blurr,edgecolor='None')
% % surf(blurr,edgecolor='None')
% % colormap(cmap_div)
% % clim([-1,1])
% % view([20 80])
% % title('Nb7Ni6-Pat_exp')
% % %%
% % figure
% % surf(Pat_Pat_exp,edgecolor='None')
% % colormap(cmap_div)
% % clim([-3,3])
% % view([20 80])
% % %%
% % Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));
% % figure
% % Pat_exp=Pat_exp-mean(Pat_exp(:));
% % Pat_exp=Pat_exp/std(Pat_exp(:));
% % imagesc(Pat_exp)
% % axis equal;
% % axis xy;
% % colormap gray;
% % %% Save the image based on the same thresholds
% % [ Pat_sim_pseudo ] = EBSP_gen(EBSD_geom_plot,RTM.Rz(60*degree)*g_refined,screen_int);
% % file_name=fullfile(export_path,'simulated_pattern.png')
% % Pat_sim_export=Pat_sim;
% % Pat_sim_export(Pat_sim_export>3.5)=3.5;
% % Pat_sim_export(Pat_sim_export<-3.5)=-3.5;
% % imwrite(flipud(rescale(Pat_sim_export)),file_name)
% % 
% % file_name=fullfile(export_path,'experimental_pattern.png')
% % Pat_exp_export=Pat_exp;
% % Pat_exp_export(Pat_exp_export>3.5)=3.5;
% % Pat_exp_export(Pat_exp_export<-3.5)=-3.5;
% % imwrite(flipud(rescale(Pat_exp_export)),file_name)
% % 
% % file_name=fullfile(export_path,'simulated_pattern_pseudo.png')
% % Pat_sim_export=Pat_sim_pseudo;
% % Pat_sim_export(Pat_sim_export>3.5)=3.5;
% % Pat_sim_export(Pat_sim_export<-3.5)=-3.5;
% % imwrite(flipud(rescale(Pat_sim_export)),file_name)
% % 
% % %%
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % figure
% % 
% % imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_sim);
% % clim([-3.5,3.5])
% % axis equal;
% % axis xy;
% % % axis off
% % colormap gray
% % 
% % 
% % %%
% % figure
% % 
% % imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_exp);
% % clim([-3.5,3.5])
% % axis equal;
% % axis xy;
% % % axis off
% % colormap gray
% 
% 
% 
% 
% % % function label_bands_alt(pattern,EBSD_geometry,m,r,InputData)
% % %     theta=(0:1:360)*degree;
% % %     m_sym=m.symmetrise('unique','antipodal');
% % % 
% % %     figure
% % %     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern)
% % %     colormap gray
% % %     axis xy;
% % %     axis off;
% % %     axis equal
% % %     hold on
% % %     fix that we plot too many bands, and transparncy gets gone
% % %     for i=1:size(m,2)
% % %         %
% % %         m_sym=m(i).symmetrise('unique','antipodal');
% % %         %
% % %         col_count=1;
% % %         for j=1:size(m_sym,1)
% % %             col_plot=InputData.cset(i,:);
% % %             col_plot=[col_plot 1];
% % %             if j==1
% % %                 plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
% % %             end
% % % 
% % %             for k=1:3
% % %                 % depending on case, plot band center or band edges
% % %                 if k==1
% % %                 omega=pi/2;
% % %                 elseif k==2
% % %                     omega = pi/2-asin(InputData.lambda_p/(2*m_sym(j).dspacing));
% % %                 else
% % %                     omega = pi/2+asin(InputData.lambda_p/(2*m_sym(j).dspacing));
% % %                 end
% % %                 plane=circels(r*m_sym(j)./norm(r*m_sym(j)),theta,omega);
% % % 
% % %             hold on
% % %             m_sym(j)
% % %             gn_coords=[plane.x./plane.z, plane.y./plane.z];
% % %             only select coordinates that are in front
% % %             gn_coords=gn_coords(plane.z>0.001,:);
% % %             get rid of everything that is bigger than the pattern
% % %             gn_coords=find_limits(gn_coords,EBSD_geometry);
% % %             plot the bands
% % %             if ~isempty(gn_coords)
% % %                 col_plot=InputData.cset(i,:)*(1-(col_count/size(m_sym,1)));
% % %                 lw=col_count*2
% % %                 col_plot=[col_plot 0.5];
% % %                 if k==1
% % %                    plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
% % %                 elseif k==2
% % %                       l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m_sym(j)));
% % %                 else
% % %                     l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
% % %                 end
% % % 
% % %                 l=patchline()
% % %                 l.lh=[0,0,0,0.5]
% % % 
% % %             end
% % %                hold on
% % %         a
% % %             end
% % %             if ~isempty(gn_coords)
% % %                 col_count=col_count+1;
% % %             end
% % %         disp([i j])
% % %         end
% % %     end
% % %     legend()
% % %         xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
% % %         ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
% % %         legend('FontSize',8)
% % %         axis off
% % % end
% % % 
% % % function label_bands(pattern,EBSD_geometry,m,r,InputData)
% % %     theta=(0:1:360)*degree;
% % %     % m_sym=m.symmetrise('unique','antipodal');
% % % 
% % %     % figure
% % %     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern)
% % %     clim([-3,3])
% % %     colormap gray
% % %     axis xy;
% % %     axis off;
% % %     axis equal
% % %     hold on
% % %     % fix that we plot too many bands, and transparncy gets gone
% % %     for i=1:size(m,2)
% % %         %%
% %         m_sym=m(i).symmetrise('unique','antipodal');
% %         %%
% %         col_count=1;
% %         % disp_name = char(m)
% %         control =1;
% %         for j=1:size(m_sym,1)
% % 
% %             col_plot=InputData.cset(i,:);
% %             col_plot=[col_plot 1];
% %             % if j==1
% %             %     plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
% %             % end
% % 
% %             for k=1:3
% %                 %% depending on case, plot band center or band edges
% %                 if k==1
% %                 omega=pi/2;
% %                 elseif k==2
% %                     omega = pi/2-asin(InputData.lambda_p/(2*m_sym(j).dspacing));
% %                 else
% %                     omega = pi/2+asin(InputData.lambda_p/(2*m_sym(j).dspacing));
% %                 end
% %                 plane=circels(r*m_sym(j)./norm(r*m_sym(j)),theta,omega);
% % 
% %             hold on
% %             % m_sym(j)
% %             gn_coords=[plane.x./plane.z, plane.y./plane.z];
% %             % only select coordinates that are in front
% %             gn_coords=gn_coords(plane.z>0.001,:);
% %             % get rid of everything that is bigger than the pattern
% %             gn_coords=find_limits(gn_coords,EBSD_geometry);
% %             %plot the bands
% %             if ~isempty(gn_coords)
% % 
% %                 col_plot=InputData.cset(i,:);%*(1-(col_count/size(m_sym,1)));
% %                 lw=2
% %                 col_plot=[col_plot 0.5];
% %                 if k==1
% %                    plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
% %                 elseif k==2
% %                     if control == 1
% %                         control =0;
% %                         l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m(i)));
% %                     else
% %                         l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
% %                     end
% %                 else
% %                     l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
% %                 end
% % 
% %                 % l=patchline()
% %                 % l.lh=[0,0,0,0.5]
% % 
% %             end
% %                hold on
% %         % a
% %             end
% %             if ~isempty(gn_coords)
% %                 col_count=col_count+1;
% %             end
% %     %     disp([i j])
% %         end
% %     end
% %     % legend()
% %         xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
% %         ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
% %         legend('FontSize',8)
% %         % axis off
% % end
% 
% % function gn_coords=find_limits(gn_coords,EBSD_geometry)
% % %%
% %     y_width=EBSD_geometry.y_gn_max-EBSD_geometry.y_gn_min;
% %     x_width=EBSD_geometry.x_gn_max-EBSD_geometry.x_gn_min;
% %     gn_coords(abs(gn_coords)<1e-4)=0;
% %     condition=gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.2) & gn_coords(:,2)<(EBSD_geometry.y_gn_max+y_width*0.2) & gn_coords(:,1)>(EBSD_geometry.x_gn_min-x_width*0.2) & gn_coords(:,1)<(EBSD_geometry.x_gn_max+x_width*0.2);
% %     %%
% %     gn_coords=gn_coords(condition,:);
% % %     gn_coords=gn_coords(find(gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.1) & gn_coords(:,2)<(EBSD_geometry.y_gn_min+y_width*0.1)));
% %    %% 
% % end
% % function c= circels(normals,theta,omega)
% % %     normals=normals(:);
% % %     c_all=zeros([length(normals),length(theta),length(omega)])
% % c_all=[];
% %     for i = 1:length(normals)
% % %%
% %       n = normals.subSet(i);
% % 
% %       % generate circles
% %       c = axis2quat(n,theta) * axis2quat(orth(n),omega)*n;
% % %       c=c';
% % %       c_all=[c_all c];
% %       % plot circles
% % %       h = [h,line(c,varargin{:},'hold','doNotDraw')]; %#ok<AGROW>
% % 
% %     end
% %     %%
% % end
% % function [lambda_out] = calc_lambda(V_in)
% % %CALC_LAMBDA calculate the wavelength of a particle, with relativistic correction
% % %
% % % V_in = voltage in eV
% % %
% % % lambda_out = output wavelength in m
% % %
% % % e.g. 
% % % [lambda_out] = calc_lambda(V_in)
% % % where:
% % % V_in=20E3; %voltage of the microscrope in eV
% % 
% % %constants
% % h=6.626E-34; %Planck constant in Js
% % me=9.1E-31; %mass of electron in kg
% % e_charge=1.602E-19; %charge of an electron in C
% % c=2.998E8; %speed of light in m/s
% % 
% % %equation with relativistic correction
% % lambda_out = ( h / sqrt(2 * me * e_charge * V_in) ) / (sqrt(1 + e_charge * V_in /2/me/c^2)); % in m
% % 
% % end
% % 
% % 
% % %%
% 
% % function p = patchline(xs,ys,varargin)
% % %https://de.mathworks.com/matlabcentral/fileexchange/36953-patchline
% % % Plot lines as patches (efficiently)
% % %
% % % SYNTAX:
% % %     patchline(xs,ys)
% % %     patchline(xs,ys,zs,...)
% % %     patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
% % %     p = patchline(...)
% % %
% % % PROPERTIES: 
% % %     Accepts all parameter-values accepted by PATCH.
% % % 
% % % DESCRIPTION:
% % %     p = patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
% % %         Takes a vector of x-values (xs) and a same-sized
% % %         vector of y-values (ys). z-values (zs) are
% % %         supported, but optional; if specified, zs must
% % %         occupy the third input position. Takes all P-V
% % %         pairs supported by PATCH. Returns in p the handle
% % %         to the resulting patch object.
% % %         
% % % NOTES:
% % %     Note that we are drawing 0-thickness patches here,
% % %     represented only by their edges. FACE PROPERTIES WILL
% % %     NOT NOTICEABLY AFFECT THESE OBJECTS! (Modify the
% % %     properties of the edges instead.)
% % %
% % %     LINUX (UNIX) USERS: One test-user found that this code
% % %     worked well on his Windows machine, but crashed his
% % %     Linux box. We traced the problem to an openGL issue;
% % %     the problem can be fixed by calling 'opengl software'
% % %     in your <http://www.mathworks.com/help/techdoc/ref/startup.html startup.m>.
% % %     (That command is valid at startup, but not at runtime,
% % %     on a unix machine.)
% % %
% % % EXAMPLES:
% % %%% Example 1:
% % %
% % % n = 10;
% % % xs = rand(n,1);
% % % ys = rand(n,1);
% % % zs = rand(n,1)*3;
% % % plot3(xs,ys,zs,'r.')
% % % xlabel('x');ylabel('y');zlabel('z');
% % % p  = patchline(xs,ys,zs,'linestyle','--','edgecolor','g',...
% % %     'linewidth',3,'edgealpha',0.2);
% % %
% % %%% Example 2: (Note "hold on" not necessary here!)
% % %
% % % t = 0:pi/64:4*pi;
% % % p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
% % % p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
% % % l = legend('sine(t)','cosine(t)');
% % % tmp = sort(findobj(l,'type','patch'));
% % % for ii = 1:numel(tmp)
% % %     set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
% % % end
% % %
% % %%% Example 3 (requires Image Processing Toolbox):
% % %%%   (NOTE that this is NOT the same as showing a transparent image on 
% % %%%         of the existing image. (That functionality is
% % %%%         available using showMaskAsOverlay or imoverlay).
% % %%%         Instead, patchline plots transparent lines over
% % %%%         the image.)
% % %
% % % img = imread('rice.png');
% % % imshow(img)
% % % img = imtophat(img,strel('disk',15));
% % % grains = im2bw(img,graythresh(img));
% % % grains = bwareaopen(grains,10);
% % % edges = edge(grains,'canny');
% % % boundaries = bwboundaries(edges,'noholes');
% % % cmap = jet(numel(boundaries));
% % % ind = randperm(numel(boundaries));
% % % for ii = 1:numel(boundaries)
% % % patchline(boundaries{ii}(:,2),boundaries{ii}(:,1),...
% % %     'edgealpha',0.2,'edgecolor',cmap(ind(ii),:),'linewidth',3);
% % % end
% % %
% % % Written by Brett Shoelson, PhD
% % % brett.shoelson@mathworks.com
% % % 5/31/2012
% % % 
% % % Revisions:
% % % 6/26 Improved rice.png example, modified FEX image.
% % % 11/30/30/2020 Fixed an issue on line 113 if varargin is empty.
% % %
% % % Copyright 2012 MathWorks, Inc.
% % %
% % % See also: patch, line, plot
% % [zs,PVs] = parseInputs(varargin{:});
% % if rem(numel(PVs),2) ~= 0
% %     % Odd number of inputs!
% %     error('patchline: Parameter-Values must be entered in valid pairs')
% % end
% % % Facecolor = 'k' is (essentially) ignored here, but syntactically necessary
% % if isempty(zs)
% %     p = patch([xs(:);NaN],[ys(:);NaN],'k');
% % else
% %     p = patch([xs(:);NaN],[ys(:);NaN],[zs(:);NaN],'k');
% % end
% % % Apply PV pairs
% % for ii = 1:2:numel(PVs)
% %     set(p,PVs{ii},PVs{ii+1})
% % end
% % if nargout == 0
% %     clear p
% % end
% % end
% % function [zs,PVs] = parseInputs(varargin)
% % if ~isempty(varargin) && isnumeric(varargin{1})
% %     zs = varargin{1};
% %     PVs = varargin(2:end);
% % else
% %     PVs = varargin;
% %     zs = [];
% % end
% % end

% function rot_mat=get_h_rot(h)