%% 17.10.2025 L. Berners
% In this skript we want to test the influence of potential misalignments
% on the comparison metrics between individual on-screen EBSD patterns.

% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;

InputUser.Phase_Folder=Phase_Folder
InputUser.Astro_loc=Astro_loc
warning('off','MATLAB:MKDIR:DirectoryExists')

% exportFolder = 'C:\Users\berners\sciebo\Paper6_OccupancyAndEffectsInNbNi\CompareMisalignment'
exportFolder = "C:\Users\lukas\sciebo\Paper6_OccupancyAndEffectsInNbNi\CompareMisalignment"
mkdir(exportFolder)
d = datetime()%t,'ConvertFrom','datenum')
d = d.Year + "-"+d.Month+"-"+d.Day+"_"+d.Hour+"-"+d.Minute +"-"+d.Second % datetime(d,'InputFormat','yyyy-MM-dd_HH-mm-ss-SSS')

%% Set up file names
% % InputUser.Phase_Input = {'Nb518Ni482_joubert_19_2025_08_01__11_26_28_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000'}
subFolder = 'NbNi_IL0.3_3aVar_ResTest'
% subFolder = "NbNi_IL_ResTest"

PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead) % -1 for debugging
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end
exportFolder = fullfile(exportFolder,subFolder)
mkdir(exportFolder)
diary( fullfile(exportFolder,sprintf("my_Diary_%s.md",d)))
% InputUser.Phase_Input  = {'NbNi_dt0.450_c26.7_a4.9_19.3kV_dmin0.020','NbNi_dt0.250_c26.7_a4.9_19.3kV_dmin0.020'}
% InputUser.Phase_Input = {'Ta4Al6Fe2_vesta','Ta4Al5Fe3_vesta','TaFeAl_20kV_hex'}
% InputUser.Phase_Input = {'YNi2'}

%% setup phase Information
for i = 1:numel(InputUser.Phase_Input)
    InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
     Phases = Phase_BuilderSimplified( InputUser.Phase_Input{i},InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    cs.mineral = Phases.xtalname
    % load ebsd map
    rtmConfig = defaultRTMConfig();
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end
%%
FHKLoptions.Voltage = double(Phases.ehistmin)*10^3; %TODO: read from phase file.
structConfig = defaultStructureFactorConfig(FHKLoptions);

rtmOptions.refine=1;
rtmOptions.clarity=1; % activate cross filter and quadrant equilib.
rtmConfigRefine = defaultRTMConfig(rtmOptions)
rtmConfigRefine.generalOptions.AlignEdax =1
% rtmConfigRefine.Settings_PCin.start = rtmConfig.Settings_PCin.start
rtmConfigRefine.gaussBlur=0


MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfigRefine    )
%%
save(fullfile(exportFolder,"MultPhaseStruct.mat"),"MultiPhaseStruct")
save(fullfile(exportFolder,"InputUser.mat"),"InputUser")
%%
rtmConfigRefine = defaultRTMConfig();
rtmConfigRefine.Settings_Cor.screensize = 514;
rtmConfigRefine.Settings_Cor.size = 514;
PC= [0.5 0.25 0.7]
rtmConfigRefine.PC = PC;
%%
% ori = orientation.rand(100);
ori = load(fullfile(exportFolder,"ori.mat"));
ori = ori.ori;
%%
save(fullfile(exportFolder,"ori.mat"),"ori");
%%
p1Name = InputUser.Phases{1};
MultiPhaseStruct(1).screen_int;
pNumCompare = 5;
metricName = 'NCC';
metrics = zeros(numel(ori),numel(InputUser.Phase_Input));
EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);  
gmatrix = ori(100).matrix;
ebsd_pat1 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNumCompare).screen_int);
figure
imagesc(ebsd_pat1)
axis equal xy
colormap gray
%% first get the "correlation without any deviations in PC and Ori
for oriNum = 1:numel(ori)
    disp(oriNum)
    gmatrix = ori(oriNum).matrix;
    ebsd_pat1 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNumCompare).screen_int);
    for pNum = 1:numel( InputUser.Phase_Input)
        

        ebsd_pat2 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNum).screen_int);
        %%
        % figure
        % imagesc(ebsd_pat1)
        % figure
        % imagesc(ebsd_pat2)
        %%
        metrics(oriNum,pNum) = getMatchMetric(ebsd_pat1,ebsd_pat2,metricName);
    end
end
%%
colorlist = [rwthPetrol(100);  rwthPetrol(75);  rwthPetrol(50); ...
            rwthTuerkis(50);  rwthTuerkis(75); rwthTuerkis(100);...
            rwthBlau(100); rwthBlau(75); rwthBlau(50);]
f=figure
f.Position=[100 100 400 400]
t = tiledlayout(3,3,"TileSpacing",'none')
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1)
binWidth = 0.002
binLimLow = 0.97
binSpacingFactor =9
for pNum = 1:numel( InputUser.Phase_Input)
    nexttile
    % [N, edges] = histcounts(metrics(:,pNum),[binLimLow:binWidth:1])
    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plName2 = extractILSpacingFromString(plName2)
    histogram(metrics(:,pNum),[binLimLow:binWidth:1],'DisplayName',sprintf('%s',plName2),'FaceAlpha',0.7,"FaceColor",colorlist(pNum,:))
    % bar((edges(1:end-1)+edges(2:end))/2+pNum*binWidth/binSpacingFactor, N, 1/binSpacingFactor,'DisplayName',sprintf('%s',plName2),'FaceAlpha',1)%binWidth/binSpacingFactor)
    hold on
    t=title(sprintf(plName2),'Position',[(binLimLow+1.005)/2,numel(ori)],'VerticalAlignment','top')
    set(t,'position',get(t,'position') -[0 0.025 0])
    % clear("t")
    xlim([binLimLow,1.005])
    ylim([0, numel(ori)])
    if pNum>6
        xlabel("NCC values")
        xticks(binLimLow:binWidth*5:1)
    else
        set(gca,'xtick',[])
    end
    if mod(pNum-1,3)==0
       ylabel("counts")
    else
       set(gca,'ytick',[])
    end
    hold on
end
% legend("Location","west")
sgtitle(sprintf("%s Nb on 3a vs.",plName1),"FontWeight","Bold")
% t.XLabel("NCC values")
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar_separated.png',subFolder)))
saveas(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar_separated.svg',subFolder)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar_separated.pdf',subFolder)))
% matrix.
%%
% colorlist = [rwthPetrol(100);  rwthPetrol(75);  rwthPetrol(50); ...
%             % rwthTuerkis(50);  rwthTuerkis(75); rwthTuerkis(100);...
%             % rwthBlau(100); rwthBlau(75); rwthBlau(50);]
figure
% t = tiledlayout("TileSpacing",'tight')
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1)
binWidth = 0.002
binLimLow = 0.97
binSpacingFactor =9
for pNum = 1:numel( InputUser.Phase_Input)
    % nexttile
    disp(oriNum)
    % [N, edges] = histcounts(metrics(:,pNum),[binLimLow:binWidth:1])
    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plName2 = extractILSpacingFromString(plName2)
    histogram(metrics(:,pNum),[binLimLow:binWidth:1],'DisplayName',sprintf('%s',plName2),'FaceAlpha',0.7)%,"FaceColor",colorlist(pNum,:))
    % bar((edges(1:end-1)+edges(2:end))/2+pNum*binWidth/binSpacingFactor, N, 1/binSpacingFactor,'DisplayName',sprintf('%s',plName2),'FaceAlpha',1)%binWidth/binSpacingFactor)
    hold on
    % title(sprintf(plName2))
    % xlim([binLimLow,1.005])
    % ylim([0, numel(ori)])
    % if pNum>6
    %     xlabel("NCC values")
    %     xticks(binLimLow:binWidth*10:1)
    % else
    %     set(gca,'xtick',[])
    % end
    % if mod(pNum-1,3)==0
    %    ylabel("counts")
    % else
    %    set(gca,'ytick',[])
    % end
    % hold on
end
legend("Location","west")
title(sprintf("%s Nb on 3a vs.",plName1))
% t.XLabel("NCC values")
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.png',subFolder)))
saveas(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.svg',subFolder)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.pdf',subFolder)))
% matrix.
%% Now try to vary the PC a bit get the "correlation without any deviations in PC and Ori
randRanges=0%[0.001 0.002 0.005]
% randRanges = [0.005
noises = [0 19.5 17 0]
% noises  =[0]

for randRange= randRanges 
for noise = noises
    diary on
    close all
    disp("###############################################")
    disp("       STarting LOOP")
    disp("random Range: ")

    disp(randRange)
    disp("")
    disp("")
    disp("Noise:")
    disp(noise)
    disp("")
    disp("")


    disp("###############################################")
    diary off
rtmConfigRefine.Settings_Cor.gaussBlur=0
offsets = zeros(3,numel(ori));
noiseRatios=zeros(1,numel(ori));
%%
% figure
% randNumbers = rand(1,3,100)
% histogram(randNumbers)
%%

% offsets
% noise =30
offsets = - randRange + rand(1,3,100)*2*randRange; %% needs to be done before loop otherwise it approaches gaussian distro
figure
% histogram(offsets)
%
% profile on
for oriNum = 1:numel(ori)
    disp(oriNum)
    EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);
    gmatrix = ori(oriNum).matrix;
    ebsd_pat1 = EBSP_gen(EBSD_geometry,gmatrix,MultiPhaseStruct(pNumCompare).screen_int);
    if noise == 0
        maxI = 2^16-1;
        meanSE = immse(ebsd_pat1,ebsd_pat1); %% Does not really make sense....
        noiseRatio = 10*log10((maxI)^2/meanSE);
    else
    [ebsd_pat1,noiseRatio] = addNoisePoisson(ebsd_pat1, noise);
    end
    ebsd_pat1 =EBSP_BGCor(ebsd_pat1,rtmConfigRefine.Settings_Cor);
    noiseRatios(oriNum) =noiseRatio;
    randOffset =  offsets(:,oriNum)'; %- randRange + rand(1,3)*2*randRange;
    % = randOffset;
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
%%
% figure
% imagesc(ebsd_pat2-ebsd_pat1)
% axis xy equal
% colormap( returnCmapDiv())
% % profile off
% % profile viewer

% saveas()
% exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.png',subFolder)))
imwrite(rescale(ebsd_pat1),fullfile(exportFolder,sprintf('ExampleNoisePattern_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
imwrite(rescale(ebsd_pat2),fullfile(exportFolder,sprintf('ExampleComparePattern_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
%%
writematrix(noiseRatios',fullfile(exportFolder,sprintf('Noise_ratios%s_pcVar%.4f_noise%.2f_comparePat%i.csv',subFolder,randRange,noise,pNumCompare)))
writematrix(squeeze(offsets)',fullfile(exportFolder,sprintf('randomOfssets%s_pcVar%.4f_noise%.2f_comparePat%i.csv',subFolder,randRange,noise,pNumCompare)))
writematrix(metrics,fullfile(exportFolder,sprintf('Metrics%s_pcVar%.4f_noise%.2f_comparePat%i.csv',subFolder,randRange,noise,pNumCompare)))
%%
% 


colorList = rwthViridis(numel(InputUser.Phase_Input));
% colorlist = colorlist(9)
% pNumCompare = 5
f=figure
f.Position = [100 100 400 400]
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1)
%
binLimLow = min(metrics(:))
binLimMax = max(metrics(:))
binLimMax =round(binLimMax,1)
if binLimLow< 0.95
    binLimLow = round(binLimLow,2)
elseif binLimLow <0.9
    binLimLow = round(binLimLow,1)

end
binWidth = (binLimMax-binLimLow)/20
%
for pNum = 1:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plName2 = extractILSpacingFromString(plName2)
    histogram(metrics(:,pNum),[binLimLow:binWidth:binLimMax],'DisplayName',sprintf('%s',plName2),'FaceAlpha',0.6,"FaceColor",colorList(pNum,:))
    hold on
end

title(sprintf(" IL %s vs.",plName1))
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
% exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_noVar.png',subFolder)))

exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_noLegend_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('Distribution_noLegend_%s_pcVar%.4f_noise%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_noLegend_%s_pcVar%.4f_noise%.2f.pdf',subFolder,randRange,noise)))
lgd = legend("Location","west")
title(lgd,'IL spacing')


exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f.pdf',subFolder,randRange,noise)))
 %%
% colorlist = [rwthPetrol(100);  rwthPetrol(75);  rwthPetrol(50); ...
%             rwthTuerkis(50);  rwthTuerkis(75); rwthTuerkis(100);...
%             rwthBlau(100); rwthBlau(75); rwthBlau(50);]

f=figure
f.Position=[100 100 400 400]
t = tiledlayout(3,3,"TileSpacing",'none')
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1)
% binWidth = 0.002
% binLimLow = 0.97
% binSpacingFactor =9
for pNum = 1:numel( InputUser.Phase_Input)
    nexttile
    % [N, edges] = histcounts(metrics(:,pNum),[binLimLow:binWidth:1])
    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plName2 = extractILSpacingFromString(plName2)
    histogram(metrics(:,pNum),[binLimLow:binWidth:1],'DisplayName',sprintf('%s',plName2),'FaceAlpha',1,"FaceColor",colorList(pNum,:))
    % bar((edges(1:end-1)+edges(2:end))/2+pNum*binWidth/binSpacingFactor, N, 1/binSpacingFactor,'DisplayName',sprintf('%s',plName2),'FaceAlpha',1)%binWidth/binSpacingFactor)
    hold on
    t=title(sprintf(plName2),'Position',[(binLimLow+binLimMax)/2,numel(ori)],'VerticalAlignment','top')
    % % % set(t,'position',get(t,'position') -[0 0.025 0])
    % clear("t")
    xlim([binLimLow,binLimMax+0.005])
    ylim([0, numel(ori)])
    if pNum>6
        xlabel("NCC values")
        % xticks(binLimLow:round((binLimMax-binLimLow)/2,3):binLimMax)
    else
        set(gca,'xtick',[])
    end
    if mod(pNum-1,3)==0
       ylabel("counts")
    else
       set(gca,'ytick',[])
    end
    hold on
end
% legend("Location","west")
sgtitle(sprintf("%s Nb on 3a vs.",plName1),"FontWeight","Bold")
% t.XLabel("NCC values")
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
% exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.png',subFolder,randRange,noise)))
% saveas(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.svg',subFolder,randRange,noise)))
% exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.pdf',subFolder,randRange,noise)))



exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('Distribution_%s_pcVar%.4f_noise%.2f_separated.pdf',subFolder,randRange,noise)))
%%

f=figure
f.Position=[100 100 400 400]
offsets = squeeze(offsets)
offsetsNormed = sqrt(offsets(1,:).^2+ offsets(2,:).^2+ offsets(3,:).^2)
[offsetsNormed_sorted,indexSorted] = sort(squeeze(offsetsNormed))
colorList = rwthViridis(numel(InputUser.Phase_Input))
for pNum = 2:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral)
    plName2 = extractILSpacingFromString(plName2)
    plot(offsetsNormed_sorted,metrics(indexSorted,pNum),'DisplayName',sprintf('%s',plName2),'LineWidth',2,"Color",colorList(pNum,:))
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
% f.Position=[100 100 400 400]
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1)
colorList = rwthViridis(numel(InputUser.Phase_Input))
f.Position = [100,100,400,400]
maxfp = 0;
for pNum = 1:numel( InputUser.Phase_Input)
    % pNum=4
    estimate = metrics(:,pNum) ;%- metrics(:,pNumCompare) ;
    if pNum ~= pNumCompare && noise ==0
        [~,fp,xfp] = kde(estimate);%,'NumPoints',100)%,NumPoints=100);
        if max(fp)>maxfp
            maxfp = max(fp(:));
        end
    end
% for pNum = 1:numel( InputUser.Phase_Input)
end


for pNum = 1:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral);
        plName2 = extractILSpacingFromString(plName2) 
    estimate = metrics(:,pNum)
    if pNum ==pNumCompare && noise ==0 && randRange==0
        xfp = ones(100,1) ;
        fp = linspace(0,maxfp,100);
        plot(xfp,fp,'DisplayName',sprintf('%s',plName2),'LineWidth',2,"Color",colorList(pNum,:))
    % end
    else
    % estimate = estimate(estimate>(estimate-0.3*std(estimate)))
    [~,fp,xfp] = kde(estimate',300);
    % histogram(metrics(:,pNum),[0.97:0.001:1],'DisplayName',sprintf('%s',plName2))
    plot(xfp,fp,'DisplayName',sprintf('%s',plName2),'LineWidth',2,"Color",colorList(pNum,:))
    end
    hold on
end
xlim([binLimLow, binLimMax]   )

title(sprintf("%s vs.",plName1))
xlabel('NCC')
ylabel('probability density')
f.Position=[100 100 400 400]
hold on
ax = gca
ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_noLegend_%s_pcVar%.4f_noise_%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_noLegend_%s_pcVar%.4f_noise_%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_noLegend_%s_pcVar%.4f_noise_%.2f.pdf',subFolder,randRange,noise)))

legend('Location','west')

exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_%s_pcVar%.4f_noise_%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_%s_pcVar%.4f_noise_%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_%s_pcVar%.4f_noise_%.2f.pdf',subFolder,randRange,noise)))
%%
f=figure
plName1 = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
plName1 = extractILSpacingFromString(plName1);
f.Position = [100,100,400,500]
% comparePhase=2
colorList = rwthViridis(numel(InputUser.Phase_Input))

maxfp = 0;
for pNum = 1:numel( InputUser.Phase_Input)
    % pNum=4
    estimate = metrics(:,pNum) - metrics(:,pNumCompare) ;
    if pNum ~= pNumCompare
        [~,fp,xfp] = kde(estimate);%,'NumPoints',100)%,NumPoints=100);
        if max(fp)>maxfp
            maxfp = max(fp(:));
        end
    end
% for pNum = 1:numel( InputUser.Phase_Input)
end
%
for pNum = 1:numel( InputUser.Phase_Input)

    plName2 = cleanNameLatex(InputUser.CSList{pNum}.mineral);
    plName2 = extractILSpacingFromString(plName2);
    estimate = metrics(:,pNum) - metrics(:,pNumCompare) ;
    if pNum == pNumCompare %%&& noise ==0
        xfp = zeros(100,1) ;
        fp = linspace(0,maxfp,100);
        plot(xfp,fp,'DisplayName',sprintf('%s',plName2),'LineWidth',2,"Color",colorList(pNum,:))
    else
        try
        % estimate = estimate(estimate>(estimate-0.3*std(estimate)))
        [~,fp,xfp] = kde(estimate,100);%,'NumPoints',100)%,NumPoints=100);
    
        % histogram(metrics(:,pNum),[0.97:0.001:1],'DisplayName',sprintf('%s',plName2))
        plot(xfp,fp,'DisplayName',sprintf('%s',plName2),'LineWidth',2,"Color",colorList(pNum,:))
        catch
            continue
        end
    end
    
    hold on
end
xlim([-0.03 0.002]   )

comparePhaseName = cleanNameLatex(InputUser.CSList{pNumCompare}.mineral)
title(sprintf("%s Nb on 3a vs. ",plName1))
xlabel('\Delta NCC')
ylabel('probability density')
hold on
ax = gca
f.Position = [100,100,400,400]
% ax.ColorOrder = defaultRWTHPalette();
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_noLegend_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_noLegend_%s_pcVar%.4f_noise%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_noLegend_%s_pcVar%.4f_noise%.2f.pdf',subFolder,randRange,noise)))
legend('Location','west')

exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_%s_pcVar%.4f_noise%.2f.png',subFolder,randRange,noise)))
saveas(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_%s_pcVar%.4f_noise%.2f.svg',subFolder,randRange,noise)))
exportgraphics(gcf,fullfile(exportFolder,sprintf('KDE_Estimation_NCCPointWiseDiff_%s_pcVar%.4f_noise%.2f.pdf',subFolder,randRange,noise)))
close all
%
end
end
%%
extractILSpacingFromString(InputUser.Phase_Input{1})
%%
function number_after_6c2= extractILSpacingFromString(str)
    pattern = '3aNb-?([0-9]*\.?[0-9]+)'; % Regular expression pattern
    matches = regexp(str, pattern, 'tokens'); % Extract matches
    
    if ~isempty(matches)
        number_after_6c2 = (matches{1}{1}); % Convert to double
        disp(number_after_6c2); % Display the result
    else
        disp('No match found.');
        number_after_6c2 = ""
    end
    % return number_after_6c2
end

