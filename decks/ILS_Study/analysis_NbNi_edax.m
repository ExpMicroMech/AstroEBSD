%% This function will load a h5oina file into MTEX and plot various plots that you might find useful
% it also has tools to load patterns
%
% written by Ben Britton 2023
% Use as is

%% File and folder infos
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;

%% Load Astro and MTEX if not loaded yet % I tried to set up a different script, so that one does not have to change the file name in every file
% InputUser.MTEX_loc='C:\Users\berners\mtex-5.7.0';
% InputUser.MTEX_loc='C:\Users\berners\MTEX_dev\mtex';
% try EBSD;
% % catch
% %     run(fullfile(MTEX_loc,"startup.m"));
% % end
%%
% try astro_loadcheck
% catch
%     run(fullfile(Astro_loc,"start_AstroEBSD.m"));
% end
% Folder that contains the h5oina file you want to load
% h5oina_folder='J:\\Users\\Berners_Lukas\\UBC\\NbNi\\pFib_ori\\EBSD\\pFib_NbNi\\h5oina';
osc_folder='I:\001Clara\FeTaAl_Laves_Christina\Ta4-Fe2-Al6_repolish\Area1'
osc_folder = "I:\MA_Dinah\Daten final\EBSD+EDX\ErNi2_S_D"
osc_folder = "J:\Users\Berners_Lukas\MasterThesisBang\EBSD\58Nb-42Ni_v1\Area1"
% osc_folder = 'J:\Users\Berners_Lukas\clara\LaNi2\LaNi2_SHL_2025_06_16\Area5_THS170'
% osc_folder = 'J:\Users\Berners_Lukas\Fib\58Nb_42Ni_v2'
% osc_folder = "E:\sciebo\IMM\BandProfiles\SiWafer\Area10-2"
%file that you want to load - omit the .h5oina
% osc_file='map20241016084212697.osc';
osc_file='NbNi_159_58Nb_42Ni_V1_nanoInd_repolish_Area 1_EBSD+EDS Map 1.ang'
% osc_file = "map20241016084212697.osc"
% osc_file = "Si_Si_with_scratch_thomas_calibrationStuff_Area 10_OIM Map 1.ang"

%location of MTEX 5.10.2
%if MTEX loaded this doesn't matter
%if MTEX is not loaded, then it will start MTEX up
% mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.10.2';
% mtex_location='C:\Users\rmb07\Documents\MATLAB\mtex-5.10.2';

%location of astroEBSD
% astro_location='C:\Users\rmb07\Documents\GitHub\AstroEBSD_v2';
% % astro_location='C:\Users\benja\OneDrive\Documents\Github\AstroEBSD';

%% Other settings
 % mb_length = 50; %micro bar length for plots - if you want to override, you can comment this out/clear this variable and the override will not happen

 
%% Load up AstroEBSD and MTEX

%start mtex if needed
% try EBSD;
% catch
%     run(fullfile(mtex_location,"startup.m"));
% end
% 
% %start astro if needed
% try astro_loadcheck
% catch
%     run(fullfile(astro_location,"start_AstroEBSD.m"));
% end

%% Load the H5 data
% [ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch] = load_h5oina_pm(h5oina_file,h5oina_folder);
% cs=loadCIF('Nb6Ni7')

% InputUser.Phase_Input={'C15'}%y,'Nb_em'}
InputUser.Phase_Input = {'Si_em'}
InputUser.Phase_Input = {'Nb6Ni7_em'}

[~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(1),Phase_Folder );
% [cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 
[cs,mineral,~]=loadCIF_atoms(RTM_info.cif_file);
%%
% cs = loadCIF(path_to_ciffile)
% cs=crystalSymmetry(cs.pointGroup,[norm(cs.axes(1)),norm(cs.axes(1)),norm(cs.axes(3))],[cs.alpha,cs.beta,cs.gamma],'x||a','Z||c*','mineral',cs.mineral)
cs = crystalSymmetry("-3m1",[4.9, 4.9, 26.7],'x||a','Z||c*','mineral','Nb-Ni mu phase')
%%
ebsd_original=EBSD.load(fullfile(osc_folder,osc_file),cs,'convertEuler2SpatialReferenceFrame','setting 2')
% ebsd_original = loadEBSD_angAstro(fullfile(osc_folder,osc_file),cs,'convertEuler2SpatialReferenceFrame','setting 2')
% ebsd_original = loadEBSd("J:\Users\Berners_Lukas\DetektorSEM_Comparison\Leo\51Nb-49Ni_v1\EBSD_51Nb-49Ni_v1")
% %if you want to plot the pattern matched data, you can change this here
%%
ebsd=ebsd_original;
ebsd("Nb6Ni7").CS=cs;
% cs = ebsd_original.CSList{2}
% ebsd=ebsd_patternmatched;
%% Set the plotting preferences - this has to be validated for your instrument

%this is set up for the pFIB and Oxford Insturments systems at UBC
setMTEXpref('xAxisDirection','east');           %aztec
setMTEXpref('zAxisDirection','IntoPlane');     %aztec
setMTEXpref('FontSize',12)
%% Now do some plots
ebsd = ebsd.gridify
%band contrast
figure;
%% 
plotx2east
plotzIntoPlane
try
plot(ebsd,ebsd.prop.iq); 
catch
    plot(ebsd,ebsd.prop.imagequality)
end
    colormap('gray');
if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_iqq.png',osc_file(1:end-4))),'Resolution',300)
%% Extract the phase to plot
plotx2east
plotzIntoPlane
figure;
try
plot(ebsd,ebsd.prop.ci); colormap('gray');
catch
    plot(ebsd,ebsd.prop.confidenceindex); colormap('gray');
end

if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ci.png',osc_file(1:end-4))),'Resolution',300)
%% Extract the phase to plot
% error('end')
%%
%select the mineral name
mineral_names = [];
for i=2:numel(ebsd.CSList)
    mineral_names= [mineral_names string(ebsd.CSList{i}.mineral)];
end
% mineral_name='Magnesium';

%% Plot the map of phases/indexed data
% figure;
% plot(ebsd');

plotx2east
plotzIntoPlane
%overlay with Band Contrast
figure;
try
plot(ebsd,ebsd.prop.iq); 
catch
    plot(ebsd,ebsd.prop.imagequality)
end
colormap('gray');
hold on;
plot(ebsd,'FaceAlpha','0.3');
hold on;
if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_iq.png',osc_file(1:end-4))),'Resolution',300)

%% IPF-x, y and z
for mineral_name = mineral_names
% compute the colors
try
    ipfKey = ipfHSVKey(ebsd(mineral_name));
    
    %plot the key
    figure;
    plot(ipfKey)
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipfColorKey_%s.png',osc_file(1:end-4),mineral_name)),'Resolution',300)
end
end
for mineral_name = mineral_names
    % if strcmp(mineral_name,"Nb6Ni7")
    try
    plotx2east
    plotzIntoPlane
    % now generate the IPF colour maps
    ipfKey.inversePoleFigureDirection = vector3d.X;
    colors_X = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors_Y = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors_Z = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    %
    %now plot the three IPF coloured maps
    figure;
    plotx2east
    plotzIntoPlane
    plot(ebsd(mineral_name),colors_X,'micronbar','off');
    title('IPF X')
    nextAxis
    plot(ebsd(mineral_name),colors_Y,'micronbar','off');
    title('IPF Y')
    nextAxis
    plot(ebsd(mineral_name),colors_Z,'micronbar','on');
    title('IPF Z')
    % end
    if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipf_%s.png',osc_file(1:end-4),mineral_name)),'Resolution',300)
    end
    end
for mineral_name = mineral_names
    try
    % now generate the IPF colour maps
    ipfKey.inversePoleFigureDirection = vector3d.X;
    colors_X = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors_Y = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors_Z = ipfKey.orientation2color(ebsd(mineral_name).orientations);
    figure;
    plotx2east
    plotzIntoPlane
    plot(ebsd(mineral_name),colors_X,'micronbar','off');
    title('IPF X')
    nextAxis
    plot(ebsd(mineral_name),colors_Y,'micronbar','off');
    title('IPF Y')
    nextAxis
    plot(ebsd(mineral_name),colors_Z,'micronbar','off');
    title('IPF Z')
    % end
    if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end %change the micronbar length to make it prettier/easier to read
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipf_%s_noMicronBar.png',osc_file(1:end-4),mineral_name)),'Resolution',300)
    end
end
%% Try to plot other properties
ebsdPropFNames = fieldnames(ebsd.prop)
for propNum = 1:numel(ebsdPropFNames)
    figure
    plot(ebsd,ebsd.prop.(ebsdPropFNames{propNum}))
    colormap gray
    colorbar
    title(ebsdPropFNames(propNum))
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_%s.png',osc_file(1:end-4),ebsdPropFNames{propNum})),'Resolution',300)
end
close all
%% Try to plot other properties - coordinates on
ebsdPropFNames = fieldnames(ebsd.prop)
for propNum = 1:numel(ebsdPropFNames)
    f=figure
    plot(ebsd,ebsd.prop.(ebsdPropFNames{propNum}),'coordinates','on')
    colormap gray
    colorbar
    title(ebsdPropFNames(propNum))
    f.Position=[100 100, 1200 1000]
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_%s_coords.png',osc_file(1:end-4),ebsdPropFNames{propNum})),'Resolution',300)%,'Padding','figure')
end
close all
% error('end')
%% Plot the PF - scatter

%plot planes as spots, and use a subset
% set
% figure
%%Miller(0,1,0,ebsd(mineral_name).CS)
for mineral_name = mineral_names
    plotPDF(ebsd(mineral_name).orientations,[Miller(1,1,0,ebsd(mineral_name).CS)  Miller(0,0,1,ebsd(mineral_name).CS)],'all','MarkerSize',1,'complete','upper')
    %you can plot all - but this takes a while..
    %plotPDF(ebsd(mineral_name).orientations,[Miller(1,1,0,ebsd(mineral_name).CS) Miller(0,1,0,ebsd(mineral_name).CS) Miller(0,0,1,ebsd(mineral_name).CS)],'Points','all')
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_texture_%s.png',osc_file(1:end-4),mineral_name)),'Resolution',300)
end

%%
Miller(1,0,-1,0,ebsd(mineral_name).CS).hkl
%% Plot the PF - scatter

%plot planes as spots, and use a subset
% set
% figure
%%Miller(0,1,0,ebsd(mineral_name).CS)

ebsd_original_2 = loadEBSD_angAstro(fullfile(osc_folder,osc_file),cs,'convertSpatial2EulerReferenceFrame','setting 2')
%%
figure
plot(ebsd_original_2,ebsd_original_2.prop.iq,'coordinates','on')
colormap gray
%%
figure
plot(ebsd_original_2,ebsd_original_2.orientations,'coordinates','on')
colormap gray
%%
figure
confirmation_index = 518995

for mineral_name = mineral_names
    plotPDF(ebsd_original_2(confirmation_index).orientations,[Miller(0,0,0,1,ebsd(mineral_name).CS)  Miller(1,0,-1,0,ebsd(mineral_name).CS),Miller(-2,1,1,0,ebsd(mineral_name).CS),Miller(-2,1,1,1,ebsd(mineral_name).CS)],'projection','eangle','all','MarkerSize',10,'complete','grid','on')%,'upper')
    %you can plot all - but this takes a while..
    % annotate(vector3d.Y,'MarkerFaceColor','green','label','y','backgroundColor','w','add2all')
    % annotate(vector3d.X,'MarkerFaceColor','green','label','x','backgroundColor','w','add2all')
    %plotPDF(ebsd(mineral_name).orientations,[Miller(1,1,0,ebsd(mineral_name).CS) Miller(0,1,0,ebsd(mineral_name).CS) Miller(0,0,1,ebsd(mineral_name).CS)],'Points','all')
    exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_texture_%s.png',osc_file(1:end-4),mineral_name)),'Resolution',300)
end

%% Plot the IPF - scatter

% figure
% plotIPDF(ebsd(mineral_name).orientations,[xvector yvector zvector])ebsd(confirmation_index).orientations
% 

% error('end')
%%
ebsd_grains=ebsd
[grains,ebsd_grains.grainId,ebsd_grains.mis2mean] = calcGrains(ebsd_grains('indexed'));

%% plot kam map
figure
plot(ebsd,ebsd.KAM('threshold',2.5*degree,'order',1) ./ degree,'micronbar','off')
setColorRange([0,2])
mtexColorbar
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',1.5)
hold off
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_KAM.png',osc_file(1:end-4))),'Resolution',600)
% end
figure
plot(ebsd,ebsd.KAM('threshold',1*degree,'order',1) ./ degree,'micronbar','off')
setColorRange([0,2])
mtexColorbar
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',1.5)
hold off
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_KAM_smallerRange.png',osc_file(1:end-4))),'Resolution',600)

figure
plot(ebsd,ebsd.KAM('threshold',2.5*degree,'order',2) ./ degree,'micronbar','off')
setColorRange([0,2])
mtexColorbar
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',1.5)
hold off
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_KAM_order2.png',osc_file(1:end-4))),'Resolution',600)
% end
figure
plot(ebsd,ebsd.KAM('threshold',1*degree,'order',2) ./ degree,'micronbar','off')
setColorRange([0,2])
mtexColorbar
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',1.5)
hold off
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_KAM_smallerRange_order2.png',osc_file(1:end-4))),'Resolution',600)
%%
%%
ebsd_grains=ebsd
ebsd_grains = ebsd_grains(ebsd_grains.prop.ci>0.1)
[grains,ebsd_grains.grainId,ebsd_grains.mis2mean] = calcGrains(ebsd_grains('indexed'));

%%
ebsd_grains(grains(grains.numPixel<=5)) = [];
[grains,ebsd_grains.grainId] = calcGrains(ebsd_grains('indexed'));
%%
figure
plot(grains)
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_grains.png',osc_file(1:end-4))),'Resolution',300)
%%
% ori=orientation.random()


%% Plotting schmid factors
% SSListNormal = [{""}]
% mu phase
SSListNormal = {{0, 0, 0, 1}'  , {1,0,-1,0}'    ,  {1,1,-2,0}'   ,  {1,-1,0,5}',    {1, 0, -1, 1}',    {1, 1, -2, 2}', {1,1,-2, 2}'}
SSListBurgers = {{1, 1, -2, 0}', {1, -2, 1,0 }' , {1, -1, 0, 0}' , {1, 1, -2, 0}', {1, -2, 1,0 }',  {1, 1, -2, -3 }',  {1, -1, 0, 0}'}

% Mg
% Basal,                            prsimatic     pyramidal a slip   pyramidal c+a  % 
% SSListNormal = {{0, 0, 0, 1}'  , {1,0,-1,0}'    ,  {1,0,-1,1}'   ,  {1,0,-1,1}',    {2,-1,-1,2}',    {1,0,-1,2}'}%, {1,1,-2, 2}'}
% SSListBurgers = {{1, 1, -2, 0}', {1, -2, 1,0 }' , {-1,2,-1,0}' , {-2,1,1,3}', {-2,1,1,3 }',  {1,0,-1,1}'}%#,  {1, -1, 0, 0}'}

% cs =ebsd.CS
for slipSysNum = 1: numel(SSListNormal)

    b = Miller(SSListBurgers{slipSysNum},cs,'UVTW')
    n = Miller(SSListNormal{slipSysNum},cs,'HKIL')
    % sS = slipSystem
    sS = slipSystem(b,n)
    
    sSSym = sS.symmetrise('antipodal')
    % r = vector3d.Z
    
    % sSSym.SchmidFactor(r)
    %
% end
    
% end 
    sSLocal  = grains.meanOrientation*sSSym
    %%
    sigma = stressTensor.uniaxial(vector3d.Z)
    SF = sSLocal.SchmidFactor(sigma);
    [SFMax,active] = max(SF,[],2);
    
    figure
    % MtexFigure('figSize','large')
    newMtexFigure('figSize','large')
    plot(grains,SFMax,'micronbar','off','linewidth',0.5)
    mtexColorbar
    clim([0,0.5])
    exportgraphics(gcf,fullfile(osc_folder,...
         sprintf('%s_SF__%s_%s.png',osc_file(1:end-4),mineral_name,char(sS))),...
         'Resolution',300)
end

%%

 % error('end')








%%









% error('end')
gB = grains("Nb6Ni7").boundary
f=figure
f.Position = [100 100 300 250]
% histogram(gB.misorientation.angle./degree,90)
% xlabel('misorientation angle (degree)')
hold on
plotAngleDistribution(gB("Nb6Ni7","Nb6Ni7").misorientation,'Displayname','Misorientation between grains','numBins',45)
%
hold on
plotAngleDistribution(cs,cs,'Displayname','uniform Misorientation distributionf')
% ylim(0,10)
f.Position = [100 100 300 250]
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_bB_misori.png',osc_file(1:end-4))),'Resolution',300)
legend('Misorientation between grains','uniform Misorientation distributionf')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_bB_misori_wLegend.png',osc_file(1:end-4))),'Resolution',300)
% f.Position = [100 100 300 300]
error('end')
%%
gB = grains.boundary
figure
% histogram(gB.misorientation.angle./degree,90)
% xlabel('misorientation angle (degree)')
hold on
plotAngleDistribution(gB.misorientation,'Displayname','Misorientation between grains','numBins',90)
hold on
plotAngleDistribution(cs,cs,'Displayname','uniform Misorientation distributionf')
ylim([0,10])
legend('Misorientation between grains','uniform Misorientation distributionf')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_bB_misori_boundMOd.png',osc_file(1:end-4))),'Resolution',300)
%%
figure
axes = gB.misorientation.axis

% define a color key
colorKey = HSVDirectionKey(axes);

% compute colors
color = colorKey.direction2color(axes);

hold on
plot(gB,'lineColor','black','linewidth',0.5) % some black background for contrast
plot(gB,color,'linewidth',0.5)
hold off
mtexColorbar('visible','off')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipf_misori_gB_map.png',osc_file(1:end-4))),'Resolution',600)
%%
figure
plot(colorKey)
hold on
plot(axes,'MarkerSize',0.1,'MarkerColor','black','all')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipf_misori_axes.png',osc_file(1:end-4))),'Resolution',300)
%%
ind = gB.misorientation.angle>59*degree & gB.misorientation.angle<61*degree;
mori = gB.misorientation(ind);
mori_mean = mean(mori,'robust')
miller1=round2Miller(mori_mean)

ind2 = gB.misorientation.angle>67*degree & gB.misorientation.angle<69*degree;
mori2 = gB.misorientation(ind2);
mori_mean2 = mean(mori2,'robust')
miller2=round2Miller(mori_mean2)
%%
% mb_length=100
gB_select=gB(ind2)
f=figure
% set(gcf, 'Position', [100 100 500 2000]);
ax=plot(gB,gB.misorientation.angle./degree,'linewidth',0.5)
clim(gca,[65,75])
mtexColorbar
if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gb60deg.png',osc_file(1:end-4))),'Resolution',600)
%%
gB_select=gB(ind2)
f=figure
f.Position= [100 100 500 2000]
hold on
plot(gB,'lineWidth',0.25)
hold on
ax=plot(gB_select,gB_select.misorientation.angle./degree,'linewidth',0.5)
if exist('mb_length','var'); mp = MTEX_mb_fix(mb_length); end
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gball.png',osc_file(1:end-4))),'Resolution',600)
% clim(gca,[65,75])
mtexColorbar
%%
figure
axes = gB(ind2).misorientation.axis

% define a color key
colorKey = HSVDirectionKey(axes);

% compute colors
color = colorKey.direction2color(axes);

hold on
plot(gB,'lineColor','black','linewidth',1) % some black background for contrast
hold on
plot(gB(ind2),color,'linewidth',2)
hold off
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gtwin2_highlighted.png',osc_file(1:end-4))),'Resolution',600)
% mtexColorbar('off')
%%
CS=ebsd.CS
figure
plot(colorKey)
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipfkey.png',osc_file(1:end-4))),'Resolution',300)

%%
figure
plot(colorKey)
hold on

% plot(axes,'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.3,'MarkerColor','black','MarkerSize',5)
% m_plot=Miller(-1,2-1,0,CS,'UVTW')
% m_plot=Miller(1,0,-1,2,'uvtw',CS)
m_plot_plane=Miller(1,0,-1,2,'hkil',CS)
hold on
m_plot_plane=Miller(0,1,-1,2,'hkil',CS)
m_plot_dir = Miller(0,2,-2,1,'hkil',CS)
m_plot=round(Miller(m_plot_plane.UVTW(1),m_plot_plane.UVTW(2),m_plot_plane.UVTW(3),m_plot_plane.UVTW(4),cs,'UVTW'))
plot(m_plot.symmetrise(),'MarkerColor','red','MarkerSize',50,'MarkerFaceAlpha',0.1,'label',char(m_plot))
plot(m_plot_plane.symmetrise(),'antipodal','DisplayName',char(m_plot),'plane')
plot(m_plot_dir.symmetrise(),'antipodal','label',char(m_plot_dir))
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_ipf_with_mori_annotated.png',osc_file(1:end-4))),'Resolution',300)
% exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gball.png',osc_file(1:end-4))),'Resolution',600)
%%

twinning=miller1;
% twinning = orientation.map(Miller(0,1,-1,-1,CS),Miller(-1,1,0,-1,CS),...
  % Miller(1,0,-1,1,CS,'uvw'),Miller(1,0,-1,-1,CS,'uvw'))
twinning2=orientation.byAxisAngle(Miller(10,0,-10,1,CS),68.5*degree)
% twinning2=miller2
%%
round(twinning.axis)
twinning.angle / degree 
%%

round(twinning2.axis)
twinning2.angle / degree
% the rotational angle

%%
figure

isTwinning = angle(gB.misorientation,twinning) < 2*degree;
twinBoundary = gB(isTwinning)

isTwinning2 = angle(gB.misorientation,twinning2) < 2*degree;
twinBoundary2 = gB(ind2)

% plot the twinning boundaries
plot(ebsd,ebsd.prop.iq)
colormap gray;
hold on
plot(grains,grains.meanOrientation,'FaceAlpha',0.3)
%plot(ebsd('indexed'),ebsd('indexed').orientations)
hold on
%plot(gB_MgMg,angle(gB_MgMg.misorientation,twinning),'linewidth',4)
% plot(twinBoundary,'linecolor','w','linewidth',2,'displayName','twin boundary, 60°, 0001')
plot(twinBoundary2,'linecolor',[1 0 1],'linewidth',1,'displayName','twin boundary, 68.5°, 10-12')
hold off

exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_tw_AxAngle.png',osc_file(1:end-4))),'Resolution',300)


%%
grainids_unique=unique(twinBoundary2.grainId,'rows');
%%
moris=[]
ori_inv=inv(grains(grainids_unique(:,1)).meanOrientation)
mori=ori_inv.*grains(grainids_unique(:,2)).meanOrientation
%%
hkl=Miller({0,0,0,1},{1,1,-2,0},CS)
figure
plotPDF(mori,hkl,'complete','all','upper')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_cluster.png',osc_file(1:end-4))),'Resolution',300)

% end
%% Plot the angle of certain Miller indices with the z plane
figure;
plotx2east
plotzIntoPlane
hkl = Miller({1,-2,0},{2,-1,0},{1,0,-1,0},CS)
compareVec = vector3d.Z
angs1 = angle(grains.meanOrientation*hkl(1).symmetrise(),compareVec)
[~,angs1_ind] = min(abs(90-angs1) ,[],2)
angle1 = angs1(angs1_ind)

angs2 = angle(grains.meanOrientation*hkl(2).symmetrise(),compareVec)
[~,angs2_ind] = min(abs(90-angs2) ,[],2)
angle2 = angs2(angs2_ind)

angs3 = angle(grains.meanOrientation*hkl(3).symmetrise(),compareVec)
[~,angs3_ind] = min(abs(90-angs3) ,[],2)
angle3 = angs3(angs3_ind)

clims = [75 105]
%%
plot(grains,angle1/degree,'micronbar','off');
colorbar()
clim(clims)
colormap(returnCmapDiv)
title(sprintf("Angle %s and %s", char(hkl(1)), char(compareVec) ))
nextAxis
plot(grains,angle2/degree,'micronbar','off');
title(sprintf("Angle %s and %s", char(hkl(2)), char(compareVec) ))
clim(clims)

colormap(returnCmapDiv)
colorbar()
nextAxis
plot(grains,angle3/degree,'micronbar','on');
colorbar()
title(sprintf("Angle %s and %s", char(hkl(3)), char(compareVec) ))
clim(clims)
colormap(returnCmapDiv)             
exportgraphics(gcf,fullfile(osc_folder,sprintf('%sAnlgeToZ.png',osc_file(1:end-4))),'Resolution',300)
%%
mori*Miller(1,0,0,cs)
mori = normalize(mori)
[cluster,center]=calcCluster(mori*normalize(Miller(1,0,0,cs)))
%%
cS = crystalShape.hex(ebsd.CS)
ori_0=orientation.byEuler(0,0,0,cs)
mori_0=inv(ori_0)*ori_0
figure
hkl=Miller({0,0,0,1},{1,1,-2,0},CS)
plotPDF(mori(cluster==1),hkl,'complete','all','upper')
hold on
plot(m_plot_plane.symmetrise,'plane','add2all')
% plot(mori_0,mori_0*cS,'add2all')
hold on
% plot(center)%,0.02*cS*mori(cluster==1),'add2all')
%%
figure
plot(ebsd,ebsd.prop.iq)
colormap gray;
hold on
plot(grains,grains.meanOrientation,'FaceAlpha',0.3)
hold on
% plot(grains(grainids_unique(cluster==1,1)))%,grains(grainids_unique(cluster==1,1)).meanOrientation)
hold on
% plot(grains(grainids_unique(cluster==1,2)))%,grains(grainids_unique(cluster==1,2)).meanOrientation)
boundary=twinBoundary2(ismember(twinBoundary2.grainId(:,1),grainids_unique(cluster==1,1)))
% plot
hold on
plot(boundary,'linecolor',[0 1 1],'linewidth',1)
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gB_clustered.png',osc_file(1:end-4))),'Resolution',300)
%%
% for
mmm=Miller(1,0,-1,2,cs)
cS = crystalShape.hex(ebsd.CS)
cS_copy=cS
% cS=crystalShape(mmm,1,1)
N = Miller(1,0,-1,2,cs)
plot(cS,N,'faceColor','blue','faceAlpha',0.5,'edgeColor','k')

sS = [slipSystem.pyramidal2CA(ebsd.CS), ...
  slipSystem.pyramidalA(ebsd.CS)]
%%

 %%
new_grain_ids=grainids_unique(cluster==1,1);
new_grain_ids2=grainids_unique(cluster==1,2);
i=1
figure
plot(grains(new_grain_ids(i)).meanOrientation*cS,'faceAlpha',0.2)
hold on
% plot(cS,sS(2),'faceColor','blue')

plotInnerFace(grains(new_grain_ids(i)).meanOrientation*cS,mmm,'faceColor','red')
drawNow(gcm,'final')
%%
 temp_boundary=grains.boundary(grains(new_grain_ids2(1)),grains(new_grain_ids(1)))
 angle(mean(temp_boundary.direction),grains(new_grain_ids2(i)).meanOrientation*mmm.symmetrise())/degree
 angle(mean(temp_boundary.direction),grains(new_grain_ids(i)).meanOrientation*mmm.symmetrise())/degree
%%
i=1
figure
plot(mean(temp_boundary.direction))
hold on 
plot(grains(new_grain_ids(i)).meanOrientation*mmm.symmetrise)
%%
% error('end')
mmm=Miller(1,0,-1,2,cs)
mm_sym=mmm.symmetrise()
cNum = 2
% selected_grains=grains(grainids_unique(cluster==1,1));
new_grain_ids=grainids_unique(:,1)%(cluster==cNum,1);
new_grain_ids2=grainids_unique(:,2)%(cluster==cNum,2);
% sel_grains2=selected_grains(selected_grains.grainSize>200)
for i=1:numel(new_grain_ids)%(cluster==cNum,1))
    if grains(new_grain_ids(i)).grainSize>2000 && grains(new_grain_ids2(i)).grainSize>2000
        figure
        plot(grains(new_grain_ids(i)),grains(new_grain_ids(i)).meanOrientation)
        hold on
        plot(grains(new_grain_ids2(i)),grains(new_grain_ids2(i)).meanOrientation)
        hold on 
        
        plot(grains(new_grain_ids(i)),grains(new_grain_ids(i)).meanOrientation*cS,'FaceAlpha',0.3)
        
        an_sel=angle(mean(temp_boundary.direction),grains(new_grain_ids(i)).meanOrientation*mmm.symmetrise())/degree
        sym_m_plot=mm_sym(find(min(an_sel)))
        sym_m_plot=Miller(sym_m_plot,cs)
        plot(grains(new_grain_ids(i)),grains(new_grain_ids(i)).meanOrientation*sym_m_plot,'faceColor','red')
        
        % plot(grains(new_grain_ids2(i)),grains(new_grain_ids2(i)).meanOrientation*cS,'FaceAlpha',0.3)
        an_sel=angle(mean(temp_boundary.direction),grains(new_grain_ids2(i)).meanOrientation*mmm.symmetrise())/degree
        sym_m_plot=mm_sym(find(min(an_sel)))
        sym_m_plot=Miller(sym_m_plot,cs)
        plot(grains(new_grain_ids2(i)),cS,grains(new_grain_ids2(i)).meanOrientation*sym_m_plot,'faceColor','red','FaceAlpha',0.3)
        
        
        
        
        % plotInnerFace(grains(new_grain_ids2(i)).meanOrientation*cS,sym_m_plot,'faceColor','red','parent',g2)
        exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gB_pair%04i_g%04i_g%04i.png',osc_file(1:end-4),i,new_grain_ids(i),new_grain_ids2(i))),'Resolution',300)
        close all
        disp(i)
        figure
        plotPDF(grains(new_grain_ids2(i)).meanOrientation,Miller({1,1,-2,0},{1,-2,1,0},{1,-1,0,2},{1,0,-1,2},{1,1,0,2},{0,2,-2,1},cs),'grid','on')
        hold on
        plotPDF(grains(new_grain_ids(i)).meanOrientation,Miller({1,1,-2,0},{1,-2,1,0},{1,-1,0,2},{1,1,0,2},{0,2,-2,1},cs))
        exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_gB_pair_PF%04i_g%04i_g%04i.png',osc_file(1:end-4),i,new_grain_ids(i),new_grain_ids2(i))),'Resolution',300)

        % error('end')

    end
end
    % grains(grainids_unique(cluster==1,1)).grainSize>200
%%
figure
cS = crystalShape.hex(ebsd.CS)
plot(grains,grains.meanOrientation)
hold on
plot(grains(twinBoundary2.grainId(:)),cS)%'FaceColor','none','linewidth',2)
plot(twinBoundary2,'linecolor','w','linewidth',6,'displayName','twin boundary, 68.5°, 10-12')
exportgraphics(gcf,fullfile(osc_folder,sprintf('%s_csOnIPF%i.png',osc_file(1:end-4),i)),'Resolution',300)