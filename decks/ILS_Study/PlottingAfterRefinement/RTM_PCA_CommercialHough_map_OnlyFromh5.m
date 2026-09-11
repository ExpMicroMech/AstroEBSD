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
warning('off','MATLAB:MKDIR:DirectoryExists')


%% color coding and FHKL plotting
cmap_div =  returnCmapDiv()

%% only plot. Do not refine.
plotOnly = 1;
InputUser.refineLatticeParams=0;
plotBandsOnly = 0;
doReprojection = 0;
plotBandsAlongOnly=0;
%% Set up file names
% InputUser.Phase_Input = {'Nb518Ni482_joubert_19_2025_08_01__11_26_28_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000'}
% subFolder =
% subFolder = 'NbNiLPR_58Nb42Ni_IL'
% subFolder = 'NbNiLPR_51Nb49Ni_IL'
subFolder = "C15Al2Compare"
% subFolder = 'NbNiLPR_51Nb49Ni_IL'
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
InputUser.Phase_Input = cell(size(PhasesToRead)-1);
for i=1:numel(PhasesToRead) % -1 for debugging
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end

% InputUser.OSC_folder = "E:\sciebo\IMM\BandProfiles\PlottingFromkH5Refinement";
% InputUser.OSC_folder = "J:\Users\Berners_Lukas\clara\YNi2MaDinah\YNi2_NS_D_repolish\Area4\Map3";
InputUser.OSC_folder = "J:\Users\Berners_Lukas\clara\YNi2MaDinah\YNi2_NS_D_repolish\Area1\Map2"
InputUser.OSC_folder = "J:\Users\Berners_Lukas\clara\C15_2_comparisonExperiment\Area1"
% InputUser.OSC_folder "= "J:\Users\Berners_Lukas\clara\NbCoUBC\48Nb52Co_verticalCut\Area2\Map1"
% InputUser.OSC_folder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\51Nb-49Ni_v2\Area10\Map1";
% InputUser.OSC_folder = "C:\Users\berners\sciebo\BandProfiles\PlottingFromkH5Refinement"

InputUser.h5Folder = fullfile(InputUser.OSC_folder,subFolder)
% InputUser.h5Folder = fullfile(InputUser.OSC_folder)
InputUser.ExportPath = InputUser.h5Folder

% InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('2026-02-13T19-23-20.551388_ANALYZED_DATA_berners_NbNi_58Nb-Ni_V2_september_Area 7_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__NbNi3aILGrid_LPR58Nb_42Ni_2026_01_cleaned.hdf5'))
% InputUser.h5name = fullfile(InputUser.ExportPath,sprintf('NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area 2_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted1_PCAMode_Tiles__NbCo3aILGrid_LPR_Low_2026_01_cleaned.hdf5'))
InputUser.h5name = fullfile(InputUser.ExportPath,sprintf("C15-MgAlCa_C15_2_Area 1_EBSD+EDS Map 1_SinglePatRefinement_eangscommercial_pcFitted-1_PCAMode_Tiles__C15Al2Compare.hdf5"))
InputUser.h5nameRepro = fullfile(InputUser.ExportPath,sprintf('SinglePatRefinementWithLPR2a.hdf5_repro'))



% PCA_Mode = 'VMoutput' , Full_Map, Tiles %switch between loading complete tiles or only the _pca_.mat files
% Tiles
% Full_Map
PCA_Mode = "Tiles";
% PCA_Mode = "Full_Map";
InputUser.PCA_Mode = PCA_Mode;
InputUser.inputEulerMode = "commercial";%"commercial"; % commercial alternatively.
InputUser.fitPC =0;
InputUser.averagedOrisAfterRefinement = 0;
refineFFT = 0;
getMapBG = 0;
realSapceOptimisation = 0;
%%
for i = 1:numel(InputUser.Phase_Input)
    InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % load ebsd map
    rtmConfig = defaultRTMConfig()
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end
if plotOnly==1
end
% InputUser=setUpPhasesFromH5(InputUser,h5info(InputUser.h5name,'/RTM'))
% load ebsd mapf

%% Get info for structure Parameters
% this assumes all phases have the same voltage
Phases = Phase_BuilderSimplified( InputUser.Phase_Input{1},InputUser.Phase_Folder );
%
try
    FHKLoptions.Voltage = double(Phases.ehistmin)*10^3; %TODO: read from phase file.
catch
    FHKLoptions.Voltage = 19.5*10^3; %20 is default typically, and we are slightly below this.
end
structConfig = defaultStructureFactorConfig(FHKLoptions);
structConfig.cutoff=0.2;
structConfig.minD = 0.5;
structConfig.maxInd = 10; % This will be slow initially but we need a max of 15 for the mu phsae % is much quicker for the plotting !
%% Set up template matching stuff
rtmOptions.refine=1;
rtmOptions.clarity=0; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions)
rtmConfigRefine.generalOptions.AlignEdax =1
rtmConfigRefine.RTM.Sampling_Freq=5;
% rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start
rtmConfigRefine.gaussBlur=1.5
rtmConfigRefine.Settings_Cor.max_var_pc = 0.1;
% rtmConfigRefine.Settings_Cor.max_var_pc = 0.1

% if we want to use pseudosymmetry. Specify as many matrices as you want.
rtmConfigRefine.options.TestPseudoSymmetry = 0;
rtmConfigRefine.options.rotPseudo = cat(3,rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree)) % for tetragonal Laves phases


%%rtmConfigRefine.RTM.Rz(60*degree),...
%%                          rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(30*degree),...
%%    rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(60*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(120*degree)) % for orthorombic laves phase Fe Ta Al
rtmConfigRefine.options.numPseudo=3;
%% otions for refinement % will influence the speed of the refinement by a lot.
optimisationOptions.MaxIter = 20;
optimisationOptions.StepTolerance = 1e-8;
optimisationOptions.OptimalityTolerance = 1e-8;
optimisationOptions.BGVarVal = 00;
optimisationOptions.Diagnostics = 'off';
optimisationOptions.Display = 'iter';
%% reprojections and Plotting
InputUser.numPatterns ="";
InputUser.reprojectionSize = 2000;
%% GetrtmConfigRefine Interpolants
% rtmConfigRefine.Settings_PCin.start =InputUser.PC %just for dummy
MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfigRefine    )
% try
%     ncc_testctToH5(rtmConfig, InputUser.h5name,"/rtmConfigRefine")
% catch
%  warning("Could not write rtmConfig to h5")
%
%% Plot refinement Results
% InputUser.averagedOrisAfterRefinement
InputUser.InputDataset ="RTMRefined";
InputUser.PatMatchMode=''
h = figure;
set(h, 'WindowStyle', 'Docked');
rtmConfigRefine.Settings_Cor.gaussian_range=1;
rtmConfigRefine.Settings_Cor.hotpixel=1
rtmConfigRefine.Settings_Cor.radius_frac = 0.69;
rtmConfigRefine.Settings_Cor.radius = 0;
% rtmConfigRefine.Settings_Cor.radius = 1;
% rtmConfigRefine.Settings_Cor.gaussBlur = 0;
rtmConfigRefine.Settings_Cor.quadrants=1; % otherwise cross does not get applied !!!!

rtmConfigRefine.Settings_Cor.cross = 2;
rtmConfigRefine.Settings_Cor.screensize = 514;
rtmConfigRefine.Settings_Cor.size = 514;
% tmConfigRefine.Settings_Cor.cross
%%
InputUser.h5nameRepro = sprintf("%s_repro",InputUser.h5name);
if plotBandsOnly~=1
    % InputUser.numPatterns = 1
    try
        InputUser = rmfield(InputUser,'numPatterns')
    end
    tileNum = 1
    plotRefinementResults(InputUser,rtmConfigRefine,structConfig,MultiPhaseStruct,tileNum)
    %
    %
end
%%

% %% Plot reprojections
% % Inpu
tileNum =0
if doReprojection==1
    plotReprojectionsFromH5(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum)
    % %%
    % %% WARNING: DO NOT RUN. All profiles will be plotted NEED debugging
    % try
end
%%
if ~plotBandsAlongOnly == 1
    InputUser.numHKLPLot=10;
    SampleBandProfilesFromReprojections(InputUser,structConfig)
end
%
% catch
%     warning('could not sample bandprofiles')
%
% end
%
% end
%
error('end')



%% END of main deck, below we have experimental stuff
%############################################################################
%##                                                                        ##
%##                                                                        ##
%##                                                                        ##
%##             END of main deck, below we have experimental stuff         ##
%##                                                                        ##
%##                                                                        ##
%##                                                                        ##
%############################################################################
%%% end
doPlots = 0
%% Test Stuff with the reprojected patterns To make comparison better work!
for phaseNumber = 1:numel(InputUser.Phase_Input)
    %
    InputUser.imDir = fullfile(InputUser.ExportPath,sprintf("Tile_%03i",tileNum),InputUser.Phase_Input{phaseNumber},'ProfilesAlongBand');
    mkdir(InputUser.imDir);
    phase_name = InputUser.Phase_Input{phaseNumber}
    expPatternsRepro = h5read(InputUser.h5nameRepro,sprintf('/RTM/%s/RTMPCA/1/ExpPatternRepro',phase_name));
    simPatternsRepro = h5read(InputUser.h5nameRepro,sprintf('/RTM/%s/RTMPCA/1/SimPatternRepro',phase_name));
    %%



    FHKLoptions.Voltage = double(Phases.ehistmin)*10^3; %TODO: read from phase file.
    structConfig = defaultStructureFactorConfig(FHKLoptions);
    FHKLoptions.Voltage = 20E3;
    structConfig = defaultStructureFactorConfig(FHKLoptions);
    % phaseNumber = 1;
    structConfig.Phase_Folder = InputUser.Phase_Folder;
    [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
    [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
    if ~isfield(InputUser,'numHKLPlot')
        InputUser.numHKLPLot=25;
    end
    [m_plot,I] =cleanFHKLDuplicates(m_plot,I,InputUser);
    structConfig.IPlot = I;



    for patNum = 1:size(expPatternsRepro,3) %only do first 5 otherwise takes too long.
        disp(sprintf("\n \n \n \n \n \n \n \n ################################################################ \n \n \n \n"))
        disp( sprintf(" Running Phase: %i,    Pattern : %i \n\n\n\n",phaseNumber,patNum))
        disp(sprintf("################################################################ \n \n \n \n"))
        
        expPattern1 = expPatternsRepro(:,:,1); %the first one is likely faulty
    
    simPattern1 = simPatternsRepro(:,:,1);
    %%
    figure
    imagesc(expPattern1)
    colormap gray
    axis('equal','off','xy')
    figure
    imagesc(simPattern1)
    colormap gray
    axis('equal','off','xy')
    %% from getBandProfilesSTereographic
    
    resolution=2/(size(expPattern1,1)-1);
    [gr_x,gr_y]=ndgrid(-1:resolution:1,-1:resolution:1);
    
    expPattern1 = expPattern1-mean(expPattern1(:));
    expPattern1 = expPattern1/std(expPattern1(:));
    simPattern1 = simPattern1-mean(simPattern1(:));
    simPattern1 = simPattern1/std(simPattern1(:));
    
    stereo_patternExp=griddedInterpolant(gr_x,gr_y,expPattern1,'cubic','none');
    stereo_patternSim=griddedInterpolant(gr_x,gr_y,simPattern1,'cubic','none');
    proj_ster=eangleProjection;
    
    
    %% setup structure factor bits for plotting the bands
    %% Get info for structure Parameters
    
    
    % i=1;
    
    
    %% Set up angular resolution
    if ~ doPlots == 0 
    
    theta=(-180:0.1:180)*degree;

    f=figure;
    f.Position = [0 0 1500 1000];
    h=tiledlayout(1,1,"TileSpacing","none");
    nexttile(1)
    imagesc(stereo_patternExp.GridVectors{1},stereo_patternExp.GridVectors{2},stereo_patternExp.Values)
    colormap gray
    axis xy
    axis equal
    axis off
    clim([-4 4])
    xlim([-1 1])
    ylim([-1 1])
    hold on

    for millerNumber =1:numel(m_plot)
        
        m = m_plot(millerNumber)
        omega = [pi/2-asin(structConfig.lambda_p/(2* m.dspacing)) pi/2+asin(structConfig.lambda_p/(2* m.dspacing))]'
        samplingPlot=circles(vector3d(m.x,m.y,m.z),theta,omega);
        samplingPlot = samplingPlot(samplingPlot.z<0)
        [x_sample_discretePlot,y_sample_discretePlot]=proj_ster.project(samplingPlot);
        scatter(x_sample_discretePlot,y_sample_discretePlot,3,'filled','DisplayName',char(m));
    end
    legend('Location','eastoutside')
    % legend()
    exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparisonBandsLabelled_%s.png',phase_name)))
    end
    %%
    theta=(0:0.1:180)*degree;
    BandWidth = 5*degree;
    omegaResolution = 0.05;
    omega=[pi/2-BandWidth:omegaResolution*degree:pi/2+BandWidth];
    omega_size=size(omega,2);
    % if you want to see an image uncomment these lines
    
    %%
    for millerNumber =1:numel(m_plot)
        m=m_plot(millerNumber);
        sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
        sampling = reshape(sampling,[],omega_size);
        thetaPlot = theta(any(sampling.z<0,2));
        sampling = sampling(any(sampling.z<0,2),:);
        
        [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
        valsExp=stereo_patternExp(y_sample_discrete,x_sample_discrete);
        valsSim=stereo_patternSim(y_sample_discrete,x_sample_discrete);
        valsExp=reshape(valsExp,[],omega_size);
        valsSim=reshape(valsSim,[],omega_size);
        valsExp = valsExp-mean(valsExp(:),'omitnan');
        valsExp = valsExp/std(valsExp(:),'omitnan');
        valsSim = valsSim-mean(valsSim(:),'omitnan');
        valsSim = valsSim/std(valsSim(:),'omitnan');
        
        %%
        %%
        
        
        %%
        % figure
        % plot(theta/degree,valsExp+omega/degree,'Color',[0.5 0.5 0.5])
        %%
        % numBins = floor(BandWidth/0.5/degree)
        binSize = 0.5;
        omegaResolution;
        binRes = floor(binSize/omegaResolution);
        numBins = floor(size(valsExp,2)/binRes);
        %
        offset = binSize;
        averageExp = zeros([size(valsExp,1),numBins]);
        averageSim = zeros([size(valsExp,1),numBins]);
        centerBin = ceil(size(valsExp,2)/2);
        omegaResolution;
        offsets = ones(1,numBins)*90;
        for binNumb = -ceil(numBins/2)+1:floor(numBins/2)-1
            
            averageExp(:,binNumb+ceil(numBins/2)+1) = mean(valsExp(:,centerBin-ceil(binRes*(-binNumb+0.5)):centerBin-floor(binRes*(-binNumb-0.5))),2);
            averageSim(:,binNumb+ceil(numBins/2)+1) = mean(valsSim(:,centerBin-ceil(binRes*(-binNumb+0.5)):centerBin-floor(binRes*(-binNumb-0.5))),2);
            offsets(binNumb+ceil(numBins/2)+1) = 90+binNumb*offset;
        end
        
        % Do some normalisation
        %theta
        averageSim = averageSim-mean(averageSim(:),'omitnan'); % correct for std
        averageSim = averageSim/std(averageSim(:),'omitnan');
        %
        averageExp = averageExp-mean(averageExp(:),'omitnan'); % correct for std
        averageExp = averageExp/std(averageExp(:),'omitnan');
        
        averageSim = averageSim/3; % correct for std
        averageSim = averageSim*binSize; % correct for bin Size
        
        averageExp = averageExp/3; % correct for std
        averageExp = averageExp*binSize; % correct for bin Size
        %%
        if max(thetaPlot/degree)>=179.99
            xlims = [90 180];
        else
            xlims = [0 90];
        end
        if ~ doPlots == 0 
        figure
        
        
        f=figure;
        f.Position =[0 0 1900 1500];
        h = tiledlayout(6,1,"TileSpacing","none");
        nexttile(1)
        imagesc(thetaPlot/degree,omega/degree,valsExp');
        axis equal;
        ylim([90-BandWidth/degree, 90+BandWidth/degree])
        xlim(xlims) % only half of it
        colormap 'gray'
        clim([-4 4])
        
        
        nexttile(2);
        imagesc(thetaPlot/degree,omega/degree,valsSim');
        axis equal;
        ylim([90-BandWidth/degree, 90+BandWidth/degree])
        xlim(xlims) % only half of it
        colormap 'gray'
        clim([-4 4])
        
        ax3 = nexttile(3);
        imagesc(thetaPlot/degree,omega/degree,valsSim'-valsExp');
        axis equal;
        ylim([90-BandWidth/degree, 90+BandWidth/degree])
        xlim(xlims) % only half of it
        colormap(ax3,returnCmapDiv)
        clim([-4 4])
        
        
        nexttile([3 1])
        
        plot(thetaPlot/degree,averageExp+offsets,'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
        xlim(xlims)
        xlabel("azimuth angle")
        ylabel("band opening angle (°) ")
        hold on
        plot(thetaPlot/degree,averageSim+offsets,'Color',[0 1 1],'DisplayName','Simulation')
        % legend()
        exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparison_%s_%s.png',phase_name,char(m_plot(millerNumber)))))
        end
        %%
        if ~ doPlots == 0 
        f=figure;
        f.Position = [100 100 1500 1000];
        for i = 1:size(averageExp,2)
            if i == 1
                plot(thetaPlot/degree,averageExp(:,i)+offsets(i),'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
                hold on
                plot(thetaPlot/degree,averageSim(:,i)+offsets(i),'Color',[0 1 1],'DisplayName','Simulation')
            else
                plot(thetaPlot/degree,averageExp(:,i)+offsets(i),'Color',[0.5 0.5 0.5],'HandleVisibility','off')
                hold on
                plot(thetaPlot/degree,averageSim(:,i)+offsets(i),'Color',[0 1 1],'HandleVisibility','off')
            end
            
        end

        index = find(offsets==90);
        % phaseNumber
        % AllBands{phaseNumber} =  averageSim(:,index);
        % en
        xlim(xlims);
        ylim([90-BandWidth/degree, 90+BandWidth/degree]);
        
        legend()
        xlabel("azimuth angle")
        ylabel("band opening angle")
        exportgraphics(f,fullfile(InputUser.imDir,sprintf('BAndComparisonOnlyBands_%s_%s.png',phase_name,char(m_plot(millerNumber)))))
        end
        %% try to get the stuff 3d
    % end
        if ~ doPlots == 0 
        viewAngle = 20;
        viewAngle2 =90;
        zlimVals = 5;
        thetaLimMinus = xlims(1);
        thetaLimPlus = xlims(2);
        cmap = returnCmapDiv;
        cmap ='gray';
        [X,Y] = meshgrid(thetaPlot/degree,omega/degree);
        f=figure
        f.Position =[0 0 500 800];
        h = tiledlayout(1,3,"TileSpacing","none");
        nexttile(1)
        surf(X,Y,valsSim','EdgeColor','none')
        title('experimental')
        axis equal
        xlim([thetaLimMinus,thetaLimPlus])
        zlim([-zlimVals,zlimVals])
        clim([-zlimVals, zlimVals])
        view(viewAngle2,viewAngle)
        colormap(cmap)
        nexttile(2)
        surf(X,Y,valsExp','EdgeColor','none')
        axis equal
        xlim([thetaLimMinus,thetaLimPlus])
        zlim([-zlimVals,zlimVals])
        clim([-zlimVals, zlimVals])
        view(viewAngle2,viewAngle)
        colormap(cmap)
        title('simulated')
        
        nexttile(3)
        surf(X,Y,valsExp'-valsSim','EdgeColor','none')
        axis equal
        xlim([thetaLimMinus,thetaLimPlus])
        zlim([-zlimVals,zlimVals])
        clim([-zlimVals, zlimVals])
        view(viewAngle2,viewAngle)
        colormap(cmap)
        title('Difference')
        exportgraphics(f, fullfile(InputUser.imDir,...
            sprintf('BAndComparison3d_%s_%s.png', phase_name, char(m_plot(millerNumber)) )...
            ))
        
        end
        %%
        close all

        %% Data export


        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/bandValuesExp", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),averageExp,"single"); %
        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/thetaValues", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),thetaPlot); %
        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/bandValuesSim", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),averageSim,"single"); %
        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/openingAngles", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),omega,"single"); %
        if phaseNumber == 1 & InputUser.averagedOrisAfterRefinement ==1  %reduce amount of data 
            writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/reprojectedPatternExp", ...
                        InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),valsExp,"single"); %        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Bandprofiles/%s/reprojectedPatternExp", ...
        end
        writeDataToH5(InputUser.h5nameRepro,sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/BandprofilesAlong/%s/reprojectedPatternSim", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber))),valsExp,"single");

        %% prepare outputs
        [m_plot_sym,indexc,~] = unique(m_plot.symmetrise(),'noSymmetry','stable');
        [~,bandGroups]=ind2sub(size(m_plot.symmetrise()),indexc);
        
        bandValues  = zeros(omega_size,numel(m_plot_sym));
        bandNames = strings(numel(m_plot_sym),1);
        for bNum = 1: numel(m_plot_sym)
            bandNames(bNum) = string(char(m_plot_sym(bNum)));
        end
        bandAngle = omega/degree;
        % bandGroups = c;
        
        %%
        % for m_number = 1:numel(m_plot_sym)
        %     % tic
        %     m = m_plot_sym(m_number);
        %
        %     sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
        %     [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
        %     vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
        %     vals=reshape(vals,[],omega_size);
        %     sum_vals=sum(vals,1)/length(theta);
        %     bandValues(:,m_number)  = sum_vals;
        %     % toc
        % end
    end
    end
end


sprintf("/RTM/%s/RTMPCA/%i/Pat_%i/Bandprofiles/%s/thetaValues", ...
                    InputUser.Phase_Input{phaseNumber},tileNum,patNum,char(m_plot(millerNumber)))




% figure
% surf(Pat_sim_Nb6Ni7-Pat_sim_Nb7Ni6,edgecolor='None')
%
% colormap(cmap_div)
% clim([-1.5,1.5])
% view([20 80])
% title('difference_simulation')
% colormap(cmap_div)
%
% %%
% figure
% blurr=imgaussfilt(Pat_sim_Nb6Ni7-Pat_sim_Nb7Ni6,5);
% surf(blurr,edgecolor='None')
%
% colormap(cmap_div)
% clim([-1.5,1.5])
% view([20 80])
% title('difference_simulation')
% colormap(cmap_div)
%
% %%
% figure
%
% surf(Pat_sim_Nb6Ni7,edgecolor='None')
% colormap(cmap_div)
% clim([-3,3])
% view([20 80])
% title('Nb6Ni7')
% axis off
% %%
% figure
% surf(Pat_sim_Nb7Ni6,edgecolor='None')
% colormap(cmap_div)
% clim([-3,3])
% view([20 80])
% title('Nb7Ni6')
% axis off
% %%
% figure
% surf(Pat_exp,edgecolor='None')
% colormap(cmap_div)
% clim([-3,3])
% view([20 80])
% axis off
% %%
% figure
% blurr=imgaussfilt(Pat_sim_Nb7Ni6-Pat_exp,5);
% surf(blurr,edgecolor='None')
% surf(blurr,edgecolor='None')
% colormap(cmap_div)
% clim([-1,1])
% view([20 80])
% title('Nb7Ni6-Pat_exp')
% %%
% figure
% surf(Pat_Pat_exp,edgecolor='None')
% colormap(cmap_div)
% clim([-3,3])
% view([20 80])
% %%
% Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));
% figure
% Pat_exp=Pat_exp-mean(Pat_exp(:));
% Pat_exp=Pat_exp/std(Pat_exp(:));
% imagesc(Pat_exp)
% axis equal;
% axis xy;
% colormap gray;
% %% Save the image based on the same thresholds
% [ Pat_sim_pseudo ] = EBSP_gen(EBSD_geom_plot,RTM.Rz(60*degree)*g_refined,screen_int);
% file_name=fullfile(tt_path,'simulated_pattern.png')
% Pat_sim_export=Pat_sim;
% Pat_sim_export(Pat_sim_export>3.5)=3.5;
% Pat_sim_export(Pat_sim_export<-3.5)=-3.5;
% imwrite(flipud(rescale(Pat_sim_export)),file_name)
%
% file_name=fullfile(export_path,'experimental_pattern.png')
% Pat_exp_export=Pat_exp;
% Pat_exp_export(Pat_exp_export>3.5)=3.5;
% Pat_exp_export(Pat_exp_export<-3.5)=-3.5;
% imwrite(flipud(rescale(Pat_exp_export)),file_name)
%
% file_name=fullfile(export_path,'simulated_pattern_pseudo.png')
% Pat_sim_export=Pat_sim_pseudo;
% Pat_sim_export(Pat_sim_export>3.5)=3.5;
% Pat_sim_export(Pat_sim_export<-3.5)=-3.5;
% imwrite(flipud(rescale(Pat_sim_export)),file_name)
%
% %%
%
%
%
%
%
%
%
%
% figure
%
% imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_sim);
% clim([-3.5,3.5])
% axis equal;
% axis xy;
% % axis off
% colormap gray
%
%
% %%
% figure
%
% imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_exp);
% clim([-3.5,3.5])
% axis equal;
% axis xy;
% % axis off
% colormap gray




% % function label_bands_alt(pattern,EBSD_geometry,m,r,InputData)
% %     theta=(0:1:360)*degree;
% %     m_sym=m.symmetrise('unique','antipodal');
% %
% %     figure
% %     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern)
% %     colormap gray
% %     axis xy;
% %     axis off;
% %     axis equal
% %     hold on
% %     fix that we plot too many bands, and transparncy gets gone
% %     for i=1:size(m,2)
% %         %
% %         m_sym=m(i).symmetrise('unique','antipodal');
% %         %
% %         col_count=1;
% %         for j=1:size(m_sym,1)
% %             col_plot=InputData.cset(i,:);
% %             col_plot=[col_plot 1];
% %             if j==1
% %                 plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
% %             end
% %
% %             for k=1:3
% %                 % depending on case, plot band center or band edges
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
% %             m_sym(j)
% %             gn_coords=[plane.x./plane.z, plane.y./plane.z];
% %             only select coordinates that are in front
% %             gn_coords=gn_coords(plane.z>0.001,:);
% %             get rid of everything that is bigger than the pattern
% %             gn_coords=find_limits(gn_coords,EBSD_geometry);
% %             plot the bands
% %             if ~isempty(gn_coords)
% %                 col_plot=InputData.cset(i,:)*(1-(col_count/size(m_sym,1)));
% %                 lw=col_count*2
% %                 col_plot=[col_plot 0.5];
% %                 if k==1
% %                    plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
% %                 elseif k==2
% %                       l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m_sym(j)));
% %                 else
% %                     l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
% %                 end
% %
% %                 l=patchline()
% %                 l.lh=[0,0,0,0.5]
% %
% %             end
% %                hold on
% %         a
% %             end
% %             if ~isempty(gn_coords)
% %                 col_count=col_count+1;
% %             end
% %         disp([i j])
% %         end
% %     end
% %     legend()
% %         xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
% %         ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
% %         legend('FontSize',8)
% %         axis off
% % end
% %
% % function label_bands(pattern,EBSD_geometry,m,r,InputData)
% %     theta=(0:1:360)*degree;
% %     % m_sym=m.symmetrise('unique','antipodal');
% %
% %     % figure
% %     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern)
% %     clim([-3,3])
% %     colormap gray
% %     axis xy;
% %     axis off;
% %     axis equal
% %     hold on
% %     % fix that we plot too many bands, and transparncy gets gone
% %     for i=1:size(m,2)
% %         %%
%         m_sym=m(i).symmetrise('unique','antipodal');
%         %%
%         col_count=1;
%         % disp_name = char(m)
%         control =1;
%         for j=1:size(m_sym,1)
%
%             col_plot=InputData.cset(i,:);
%             col_plot=[col_plot 1];
%             % if j==1
%             %     plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
%             % end
%
%             for k=1:3
%                 %% depending on case, plot band center or band edges
%                 if k==1
%                 omega=pi/2;
%                 elseif k==2
%                     omega = pi/2-asin(InputData.lambda_p/(2*m_sym(j).dspacing));
%                 else
%                     omega = pi/2+asin(InputData.lambda_p/(2*m_sym(j).dspacing));
%                 end
%                 plane=circels(r*m_sym(j)./norm(r*m_sym(j)),theta,omega);
%
%             hold on
%             % m_sym(j)
%             gn_coords=[plane.x./plane.z, plane.y./plane.z];
%             % only select coordinates that are in front
%             gn_coords=gn_coords(plane.z>0.001,:);
%             % get rid of everything that is bigger than the pattern
%             gn_coords=find_limits(gn_coords,EBSD_geometry);
%             %plot the bands
%             if ~isempty(gn_coords)
%
%                 col_plot=InputData.cset(i,:);%*(1-(col_count/size(m_sym,1)));
%                 lw=2
%                 col_plot=[col_plot 0.5];
%                 if k==1
%                    plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
%                 elseif k==2
%                     if control == 1
%                         control =0;
%                         l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m(i)));
%                     else
%                         l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
%                     end
%                 else
%                     l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
%                 end
%
%                 % l=patchline()
%                 % l.lh=[0,0,0,0.5]
%
%             end
%                hold on
%         % a
%             end
%             if ~isempty(gn_coords)
%                 col_count=col_count+1;
%             end
%     %     disp([i j])
%         end
%     end
%     % legend()
%         xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
%         ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
%         legend('FontSize',8)
%         % axis off
% end

% function gn_coords=find_limits(gn_coords,EBSD_geometry)
% %%
%     y_width=EBSD_geometry.y_gn_max-EBSD_geometry.y_gn_min;
%     x_width=EBSD_geometry.x_gn_max-EBSD_geometry.x_gn_min;
%     gn_coords(abs(gn_coords)<1e-4)=0;
%     condition=gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.2) & gn_coords(:,2)<(EBSD_geometry.y_gn_max+y_width*0.2) & gn_coords(:,1)>(EBSD_geometry.x_gn_min-x_width*0.2) & gn_coords(:,1)<(EBSD_geometry.x_gn_max+x_width*0.2);
%     %%
%     gn_coords=gn_coords(condition,:);
% %     gn_coords=gn_coords(find(gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.1) & gn_coords(:,2)<(EBSD_geometry.y_gn_min+y_width*0.1)));
%    %%
% end
% function c= circels(normals,theta,omega)
% %     normals=normals(:);
% %     c_all=zeros([length(normals),length(theta),length(omega)])
% c_all=[];
%     for i = 1:length(normals)
% %%
%       n = normals.subSet(i);
%
%       % generate circles
%       c = axis2quat(n,theta) * axis2quat(orth(n),omega)*n;
% %       c=c';
% %       c_all=[c_all c];
%       % plot circles
% %       h = [h,line(c,varargin{:},'hold','doNotDraw')]; %#ok<AGROW>
%
%     end
%     %%
% end
% function [lambda_out] = calc_lambda(V_in)
% %CALC_LAMBDA calculate the wavelength of a particle, with relativistic correction
% %
% % V_in = voltage in eV
% %
% % lambda_out = output wavelength in m
% %
% % e.g.
% % [lambda_out] = calc_lambda(V_in)
% % where:
% % V_in=20E3; %voltage of the microscrope in eV
%
% %constants
% h=6.626E-34; %Planck constant in Js
% me=9.1E-31; %mass of electron in kg
% e_charge=1.602E-19; %charge of an electron in C
% c=2.998E8; %speed of light in m/s
%
% %equation with relativistic correction
% lambda_out = ( h / sqrt(2 * me * e_charge * V_in) ) / (sqrt(1 + e_charge * V_in /2/me/c^2)); % in m
%
% end
%
%
% %%

% function p = patchline(xs,ys,varargin)
% %https://de.mathworks.com/matlabcentral/fileexchange/36953-patchline
% % Plot lines as patches (efficiently)
% %
% % SYNTAX:
% %     patchline(xs,ys)
% %     patchline(xs,ys,zs,...)
% %     patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
% %     p = patchline(...)
% %
% % PROPERTIES:
% %     Accepts all parameter-values accepted by PATCH.
% %
% % DESCRIPTION:
% %     p = patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
% %         Takes a vector of x-values (xs) and a same-sized
% %         vector of y-values (ys). z-values (zs) are
% %         supported, but optional; if specified, zs must
% %         occupy the third input position. Takes all P-V
% %         pairs supported by PATCH. Returns in p the handle
% %         to the resulting patch object.
% %
% % NOTES:
% %     Note that we are drawing 0-thickness patches here,
% %     represented only by their edges. FACE PROPERTIES WILL
% %     NOT NOTICEABLY AFFECT THESE OBJECTS! (Modify the
% %     properties of the edges instead.)
% %
% %     LINUX (UNIX) USERS: One test-user found that this code
% %     worked well on his Windows machine, but crashed his
% %     Linux box. We traced the problem to an openGL issue;
% %     the problem can be fixed by calling 'opengl software'
% %     in your <http://www.mathworks.com/help/techdoc/ref/startup.html startup.m>.
% %     (That command is valid at startup, but not at runtime,
% %     on a unix machine.)
% %
% % EXAMPLES:
% %%% Example 1:
% %
% % n = 10;
% % xs = rand(n,1);
% % ys = rand(n,1);
% % zs = rand(n,1)*3;
% % plot3(xs,ys,zs,'r.')
% % xlabel('x');ylabel('y');zlabel('z');
% % p  = patchline(xs,ys,zs,'linestyle','--','edgecolor','g',...
% %     'linewidth',3,'edgealpha',0.2);
% %
% %%% Example 2: (Note "hold on" not necessary here!)
% %
% % t = 0:pi/64:4*pi;
% % p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
% % p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
% % l = legend('sine(t)','cosine(t)');
% % tmp = sort(findobj(l,'type','patch'));
% % for ii = 1:numel(tmp)
% %     set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
% % end
% %
% %%% Example 3 (requires Image Processing Toolbox):
% %%%   (NOTE that this is NOT the same as showing a transparent image on
% %%%         of the existing image. (That functionality is
% %%%         available using showMaskAsOverlay or imoverlay).
% %%%         Instead, patchline plots transparent lines over
% %%%         the image.)
% %
% % img = imread('rice.png');
% % imshow(img)
% % img = imtophat(img,strel('disk',15));
% % grains = im2bw(img,graythresh(img));
% % grains = bwareaopen(grains,10);
% % edges = edge(grains,'canny');
% % boundaries = bwboundaries(edges,'noholes');
% % cmap = jet(numel(boundaries));
% % ind = randperm(numel(boundaries));
% % for ii = 1:numel(boundaries)
% % patchline(boundaries{ii}(:,2),boundaries{ii}(:,1),...
% %     'edgealpha',0.2,'edgecolor',cmap(ind(ii),:),'linewidth',3);
% % end
% %
% % Written by Brett Shoelson, PhD
% % brett.shoelson@mathworks.com
% % 5/31/2012
% %
% % Revisions:
% % 6/26 Improved rice.png example, modified FEX image.
% % 11/30/30/2020 Fixed an issue on line 113 if varargin is empty.
% %
% % Copyright 2012 MathWorks, Inc.
% %
% % See also: patch, line, plot
% [zs,PVs] = parseInputs(varargin{:});
% if rem(numel(PVs),2) ~= 0
%     % Odd number of inputs!
%     error('patchline: Parameter-Values must be entered in valid pairs')
% end
% % Facecolor = 'k' is (essentially) ignored here, but syntactically necessary
% if isempty(zs)
%     p = patch([xs(:);NaN],[ys(:);NaN],'k');
% else
%     p = patch([xs(:);NaN],[ys(:);NaN],[zs(:);NaN],'k');
% end
% % Apply PV pairs
% for ii = 1:2:numel(PVs)
%     set(p,PVs{ii},PVs{ii+1})
% end
% if nargout == 0
%     clear p
% end
% end
% function [zs,PVs] = parseInputs(varargin)
% if ~isempty(varargin) && isnumeric(varargin{1})
%     zs = varargin{1};
%     PVs = varargin(2:end);
% else
%     PVs = varargin;
%     zs = [];
% end
% end

% function rot_mat=get_h_rot(h)