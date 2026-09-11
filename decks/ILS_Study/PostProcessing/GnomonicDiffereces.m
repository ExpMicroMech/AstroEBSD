clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).

%% --------- Specify Phases ----------------------------------------
InputUser.Phase_Folder=Phase_Folder

% InputUser.Phase_Input={'Nb_em','Nb_em_10kV','Nb_em_5kV','Ni_em_10kV','Ni','Ni_em_5kV','NbNi3_em','NbNi3_em_5kV','NbNi3_em_10kV','Nb6Ni7_em','Nb6Ni7_em_5kV','Nb6Ni7_em_10kV'};
% InputUser.Phase_Input={'Cu','Ni_em20','10Cu90Ni','20Cu80Ni','30Cu70Ni','40Cu60Ni','50Cu50Ni','60Cu40Ni','70Cu30Ni','80Cu20Ni','90Cu10Ni'}
% InputUser.Phase_Input={'Nb7Ni6_em','Nb6Ni7_em'}%,Ni_em20','10Cu90Ni','20Cu80Ni','30Cu70Ni','40Cu60Ni','50Cu50Ni','60Cu40Ni','70Cu30Ni','80Cu20Ni','90Cu10Ni'}
% InputUser.Phase_Input={'Au3Ni','75Au25Ni'}
% InputUser.Phase_Input={'Nb7Ni6_em_lattice_fixed','Nb6Ni7_em_lattice_fixed'}
% InputUser.Phase_Input={'Si_HR_lowd'}
% InputUser.Phase_Input = {'NbNi_dt0.400_c26.7_a4.9_3aNb_c1-8_c2-20'}%,'NbNi_dt0.400_c26.7_a4.9_3aNb_c2-12','NbNi_dt0.400_c26.7_a4.9_3aNb_c2-20'}
% InputUser.Phase_Input = {'La7Ni16_30kV','La15Ni32_30kV','LaNi2_30kV'}
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.5_a4.9_20kV_2025_07_22__10_59_34_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% InputUser.Phase_Input = {'NbNi_dt0.400_c26.9_a4.9_20kV_2025_07_21__16_01_10_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.250_c26.9_a4.9_20kV_2025_07_21__16_01_05_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.250_c26.5_a4.9_20kV_2025_07_21__16_01_00_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.400_c26.5_a4.9_20kV_2025_07_21__16_00_52_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% InputUser.Phase_Input = {'Nb518Ni482_joubert_20.0kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.000_npx-1000',...
%                          'Nb518Ni482_joubert_20.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.800_npx-1000',...
%                         'Nb518Ni482_joubert_20.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.300_npx-1000',...
%                         'Nb518Ni482_joubert_19.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.800_npx-1000',...
%                         'Nb518Ni482_joubert_19.6kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.600_npx-1000',...
%                         'Nb518Ni482_joubert_19.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...}
%                         'Nb518Ni482_joubert_18.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.800_npx-1000',...}
%                         'Nb518Ni482_joubert_18.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.300_npx-1000',...
%                         'Nb518Ni482_joubert_17.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-17.800_npx-1000'
% InputUser.Phase_Input ={'Nb518Ni482_joubert_20_2025_07_30__21_40_49_dmin-0.020_c1-8.000_c2-30.000_ehistmin-19.800_npx-1000',...
%    'Nb518Ni482_joubert_20_2025_07_30__19_17_43_dmin-0.020_c1-6.000_c2-30.000_ehistmin-19.800_npx-1000',...
%    'Nb518Ni482_joubert_20.0kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.000_npx-1000'}
% InputUser.Phase_Input={'Nb569Ni431_joubert_20kV_2025_07_24__08_56_35_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}%,'Nb569Ni431_joubert',}
% InputUser.Phase_Input = {'Ta4Al5Fe3_vesta','Ta4Al6Fe2_vesta'}%,'Nb518Ni482_joubert_20_2025_07_31__16_58_58_dmin-0.020_c1-6.000_c2-20.000_ehistmin-19.800_npx-1000'}
subFolder = 'YNi2_SLO'
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))
% temp = PhasesToRead{2}

InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead)
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end

temp = InputUser.Phase_Input{2}
% InputUser.Phase_Input{2} =InputUser.Phase_Input{1}
% InputUser.Phase_Input{1} = temp
% InputUser.Phase_Input={'Nb569Ni431_joubert_20kV_2025_07_24__08_56_35_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}%,'Nb569Ni431_joubert',}
% PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/Nb518Ni482_joubert_*_2025_08_01_*.pha'))'
% PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/Y*Ni*_dmin*.pha'))'

% InputUser.Phase_Input = cell(size(PhasesToRead));
% for i=1:numel(PhasesToRead)
%     InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
% end
% InputUser.Phase_Input= {'75Au25Ni','Au3Ni'}
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.7_a4.9_19_2025_08_04__10_23_33_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.275_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.325_c26.7_a4.9_19_2025_08_04__10_23_35_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.300_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.350_c26.7_a4.9_19_2025_08_04__10_23_36_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.375_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.400_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     }
% cat(PhasesToRead.name)
% InputUser.Phase_Input={''}
% InputUser.Phase_Input={'Au','Ni_em20','AuNi3','10Au90Ni','20Au80Ni','30Au70Ni','40Au60Ni','50Au50Ni','60Au40Ni','70Au30Ni','80Au20Ni','90Au10Ni'}
root_path=fullfile("E:\sciebo\IMM\EBSD\YNi2\compareGnomonics_Gfilt",subFolder)
root_path = fullfile("I:\sciebo\EBSD\YNi2\compareGnomonics_Gfilt\YNi2_dbwCompare")
% root_path = 'I:\002Simulation_EBSD\YNiCompare'
% root_path = 
mkdir(root_path)
% root_path='E:\RWTH_Aachen\Promo\CIF\ecp'

%% save some stuff to structured data
InputUser.h5name = fullfile(root_path,sprintf('GnomonicDifferences_%s.hdf5',subFolder))

try
    makeNewH5(InputUser.h5name)
catch
    warning("File could not be created, mak sure, you have setup git properly on your machine")
end



%% ---------------------------- Set up Zone Axes to put in view-------------
% h=Miller({0,0,1},{1,1,1},{0,1,1},cs);
% h=Miller('UVTW',cs)
% hList = {{0 0 0 1}',{-1,0,1,0}',{-1,2,-1,0}',{1, 1, -2, 1}',{3, 0, -3, 1}',}
hList = {{0 0 1},{1 0 0},{0 1 0},{0 0 -1},{-1 0 0},{0 -1 0}}
hList = {{0 0 1},{1 1 0},{1 1 2},{0 0 -1},{1 1 1}}%,{0 -1 0}}
hType = 'UVW';%fullfil
hnames = strings(numel(hList),numel(InputUser.Phase_Input))
emptyHs = cell(numel(InputUser.Phase_Input),1)
% h is for the zone axis along which we want to observe the ecpori=[]


numOris = numel(hList)








%% Set up plotting resolution and parameters ------------------------------
% [InputData.lambda] = calc_lambda(InputData.V_in);

%scale the wavelength to the same units as the dhkl
% InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale; 
%lambda normalised by the scale factor of the unit cell
patterns=zeros(514,514,numOris); 
%foruth number for size of millers to plot
rtmConfigRefine = defaultRTMConfig();
rtmConfigRefine.Settings_Cor.screensize = 514;
rtmConfigRefine.Settings_Cor.size = 514;
rtmConfigRefine.Settings_Cor.gaussBlur = 1;



%% radon transform settings
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak hunt
Settings_Rad.max_peaks=12; %max number of peaks to return
Settings_Rad.num_peak=20; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=6; %search size in theta steps
Settings_Rad.rho_search_per=0.1; %radon search in fractions
Settings_Rad.min_peak_width=0.002;


% rtmConfigRefineRadon


rtmConfigRefineRadon = rtmConfigRefine;
rtmConfigRefineRadon.Settings_Cor.gaussianRange = 2;
rtmConfigRefineRadon.Settings_Rad = Settings_Rad;


PC= [0.5 0.5 0.65]
PC_offset = [0.000 0.000 0.000]
rrmtConfigRefine.PC = PC;



try
    writeStructToH5(rtmConfigRefine, InputUser.h5name,"/UserInfo")
catch
 warning("Could not write Input Info to h5")
end

try
    writeStructToH5(InputUser, InputUser.h5name,"/PlottingConfig")
catch
 warning("Could not write Input Info to h5")
end




%% ------------- Plot all patterns (gnomonic and sphere) -----------
phaseSize = numel(InputUser.Phase_Input);

%% 
% f
for phaseNumber=1:phaseSize
    %% ------------------ set up Phase information ----------------------
    % try
    phase_name=InputUser.Phase_Input{phaseNumber};
% end
    %% 
    if iscell(phase_name)
        phase_name = phase_name{1};
    end
    phase_namePlot = cleanNameLatex(phase_name);
    try
    [ ~,~,~,~,~, RTM_info ] = ...
            Phase_Builder_RTM( {phase_name},InputUser.Phase_Folder );
    catch
        warning('Phase %s could not be determined', phase_name)
        continue;
    end


% end
    % -------------- get FHKL params-----------------------------------
    RTM_info.Voltage = 19.5E3;
    structConfig = defaultStructureFactorConfig(RTM_info);
    structConfig.Phase_Folder = InputUser.Phase_Folder;
    structConfig.maxInd = 10;
    %%
    % R
    [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
    %%
    [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
    [m_plot,I] =cleanFHKLDuplicates(m_plot,I);
    structConfig.IPlot = I;
% end
    [cs,~]=loadCIF(RTM_info.cif_file);

    [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
    %
    try
        RootExport_path=fullfile(root_path,phase_name,sprintf("%.2fkV",RTM_info.Voltage));
    catch
        RootExport_path=fullfile(root_path,phase_name,"defalutkV");
    end
    mkdir(RootExport_path);
    % 
% end
    %% ------------ Plot a sphere -------------------
    export_path = fullfile(RootExport_path,'sphere');
    mkdir(export_path);
    if phaseNumber ==2
        PC = PC +PC_offset
    end
    % disp('until heree')
    EBSD_geometry=EBSP_Gnom(rtmConfigRefine.Settings_Cor,PC);  
    % disp('until heree')
    [x,y,z]=sphere(2000);
    i_data=Cube_Sample(x,y,z,screen_int,screen_int.isHex);
    i_data = i_data-mean(i_data(:));
    i_data = i_data/std(i_data(:));
    % disp('until heree')
% end
    % f = figure;
    % surf(x,y,z,i_data,'EdgeColor','none')
    % colormap gray;
    % axis equal;
    % title(phase_namePlot)
    % clim([-4 4])
    % exportgraphics(f,fullfile(export_path,sprintf('%s_sphere.png', ...
    %          phase_name)))
    % print(f, ['-d' 'png'], '-r920', ...
    %         fullfile(export_path,sprintf('%s_sphere.png', ...
    %          phase_name)));
    close all
    % disp('until heree')
% end
% end
    %% ---------------- Build rotation matrices -------------------------
    clear ori;
    ori=[];
    h = Miller(hList{:},hType,cs);
    hnames = [];
    for m_num=1:numOris
        rot_mat=get_h_rot(h(m_num));
        % rot_mat=get_h_rot
        ori=[ori rot_mat];
        
        hnames = [hnames string(char(h(m_num)))];
    end
    emptyHs{phaseNumber} =h;
% en
    
    %% Export images
    export_path = fullfile(RootExport_path,'gnomonic')
    mkdir(export_path);
    export_pathBands = fullfile(RootExport_path,'gnomonicWBands')
    mkdir(export_pathBands);
% end
    % hnames = []
    for oriNum=1:numOris
        
        r=inv(ori(oriNum).matrix);
        % euler=ori(oriNum).Euler();
        pattern=EBSP_gen(EBSD_geometry,r,screen_int,0);
        pattern = EBSP_BGCor(pattern,rtmConfigRefine.Settings_Cor);
        % pattern = pattern-mean(pattern(:));
        % pattern = pattern./std(pattern(:));
        % patt
        patterns(:,:,oriNum)=pattern;

        fname=fullfile(export_path, sprintf( '%s_gnomonic_%s_normalised.png' ,phase_name, char(h(oriNum))));
        imwrite(normaliseForImwrite((flipud(pattern))),fname)
        fname=fullfile(export_path, sprintf( '%s_gnomonic_%s.png' ,phase_name, char(h(oriNum))));
        imwrite(rescale((flipud(pattern))),fname)

        figure
        pPattern(pattern,EBSD_geometry)
        clim([-3.5,3.5])
        xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        ylim([EBSD_geometry.y_gn_min,EBSD_geometry.y_gn_max])
        % 
        % label_bandsStructureFactor3(pattern,EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),structConfig)
        % title(sprintf("%s \n %s",phase_namePlot, char(h(oriNum))) )
        % fname2 = fullfile(export_pathBands, sprintf( '%s_gnomonic_%s_differnce.png' ,phase_name, char(h(oriNum))));
        % exportgraphics(gcf,fname2)
        phase_name = InputUser.Phase_Input{phaseNumber}
 close all     
        %%
    end
    %% --------- save Patterns to structured h5 ----------------
    try 
        voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
    catch
        voltageString = "defaultVoltage"
    end
    %%
    h5PhaseInfo=    Phase_BuilderSimplified(phase_name,Phase_Folder)
    h5PhaseInfo.phaName =string(phase_name);
    writeStructToH5(h5PhaseInfo,InputUser.h5name,sprintf("/%s/%s/SimulationInfo",voltageString,h5PhaseInfo.phaName))

    writeDataToH5(InputUser.h5name,sprintf("/%s/%s/Patterns",voltageString,h5PhaseInfo.phaName),...
        squeeze(patterns));
    writeDataToH5(InputUser.h5name,sprintf("/%s/%s/ZoneAxes",voltageString,h5PhaseInfo.phaName),hnames,"string")
    % writeDataToH5(InputUser.h5name,sprintf("/%s/Patterns",phase_name),...
        % squeeze(patterns(:,:,:)),"single");
    % hnames(hnum,phaseNumber) = string(char(h(hnum)));
end

 








%% --------------- compare differences --------------------------
comparisonPhase = 4; % aVoid magic numers :D





% set up some more colormap and plotting variables
% limits
differenceLimits = [-2,2];
differenceLimitsNeg = [-2,0];
differenceLimitsPos = [0,2]
ssimLimits = [-1,1];


%colormaps
differenceColormap = returnCmapDiv;
imageColormap = 'gray';

[cmap_blues]=cbrewer('seq','Blues',30);
cmap_blues = flipud(cmap_blues);
cmap_blues(cmap_blues<0)=0;

[cmap_reds]=cbrewer('seq','Reds',30);
% cmap_reds=flipud(cmap_reds)
cmap_reds(cmap_reds<0)=0;
[cmap_spectral]=cbrewer('div','Spectral',30);
% cmap_reds=flipud(cmap_reds)
cmap_spectral(cmap_spectral<0)=0;


ssim_colormap=cbrewer('div','PuOr',30);
% cmap_reds=flipud(cmap_reds)
ssim_colormap(ssim_colormap<0)=0;
ssim_colormap(ssim_colormap>1)=1;

NCCValues = zeros(numel(hList),1);
SSIMValues = zeros(numel(hList),1);
SSIMLuminanceValues = zeros(numel(hList),1);
SSIMStructureValues = zeros(numel(hList),1);
SSIMContrastValues  = zeros(numel(hList),1);
phase_nameFirst = InputUser.Phase_Input{comparisonPhase};
phase_nameFirstPlot = cleanNameLatex(phase_nameFirst);

   [ ~,~,~,~,~, RTM_info ] = ...
            Phase_Builder_RTM( {phase_nameFirst},InputUser.Phase_Folder );
 try 
        voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
    catch
        voltageString = "defaultVoltage"
 end
   % voltageString = "19.30kV"
 h5PhaseInfo=    Phase_BuilderSimplified(phase_nameFirst,Phase_Folder)
%%
comparisonPatterns = double(h5read(InputUser.h5name,sprintf("/%s/%s/Patterns",voltageString,phase_nameFirst)));
%%
for phaseNumber=1:numel(InputUser.Phase_Input)
   
    
    
    phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name);
    mkdir(RootExport_path);
    %% load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
  % h5PhaseInfo
     % voltageString = "19.30kV"
    h5PhaseInfo =  Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)
    patterns = double(h5read(InputUser.h5name,sprintf("/%s/%s/Patterns",voltageString,string(phase_name))));
    for oriNum= 1 : numOris



        

   
        export_path = fullfile(RootExport_path,'Difference');
        mkdir(export_path);
        
        ncc = getMatchMetric(comparisonPatterns(:,:,oriNum),...
                    patterns(:,:,oriNum),"NCC");
        
        NCCValues(:,oriNum) = max(ncc(:));
        %% ------------  difference-----------------
        h = emptyHs{phaseNumber}
        diff_pat = comparisonPatterns(:,:,oriNum) - ...
            patterns(:,:,oriNum);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference.png'...
                        ,phase_name, char(h(oriNum))));
        plotTitle = [char(h(oriNum)),sprintf('\n%s \n- %s',...
                        phase_nameFirstPlot,phase_namePlot)];
        plotComparison(diff_pat,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
        
        % export gray value image
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_gray.png'...
                        ,phase_name, char(h(oriNum))));
        imwrite(normaliseForImwrite(diff_pat),fname2);
        
     %% ------------   Radon Transform-----------------
        theta = -90:0.2:89; % bad practise but I am lazy.
        % rtmConfigRefineRadon
        % EBSP_BGCor()
        % BW_comp = edge(EBSP_BGCor(comparisonPatterns(:,:,oriNum));
        % BW_pat = edge(EBSP_BGCor(patterns(:,:,oriNum))
        [R_comp,xp] = radon(EBSP_BGCor(comparisonPatterns(:,:,oriNum),...
                    rtmConfigRefineRadon.Settings_Cor),theta);
        [R_pat,xp] = radon(EBSP_BGCor(patterns(:,:,oriNum),...
                    rtmConfigRefineRadon.Settings_Cor),theta);
        %% 
        figure
        R_comp = R_comp-mean(R_comp(:));
        R_comp = R_comp./std(R_comp(:));
        R_pat = R_pat-mean(R_pat(:));
        R_pat = R_pat./std(R_pat(:));
        % [R,xp] = radon(Pat_sim,-90:0.2:89);
        imagesc(theta,xp,R_comp)
        xlabel('\theta (degrees)')
        ylabel('x''')
        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s.png'...
                        ,phase_nameFirst, char(h(oriNum))));        
        plotTitle = [char(h(oriNum)),sprintf('\n%s',...
                        phase_nameFirstPlot)];
        title(plotTitle);
        xlabel('\theta (degrees)')
        ylabel('x''')
        colormap gray
        exportgraphics(gcf,fname2)

        % h = emptyHs{phaseNumber}
        
        % diff_pat = comparisonPatterns(:,:,oriNum) - ...
        %     patterns(:,:,oriNum);

        % plotComparison(R_comp,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
        
        imagesc(theta,xp,R_pat)
        xlabel('\theta (degrees)')
        ylabel('x''')
        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s.png'...
                        ,phase_name, char(h(oriNum))));
        plotTitle = [char(h(oriNum)),sprintf('\n%s',...
                        phase_namePlot)];
        title(plotTitle);
        xlabel('\theta (degrees)')
        ylabel('x''')
        colormap gray
        exportgraphics(gcf,fname2)
        % plotComparison(R_pat,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
        

        % diff_pat_r = R_comp - R_pat;
        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s_difference.png'...
                        ,phase_name, char(h(oriNum))));
        plotTitle = [char(h(oriNum)),sprintf('\n%s \n- %s',...
                       phase_nameFirstPlot,phase_namePlot)];
        % plotComparison(diff_pat,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
 

        imagesc(theta,xp,R_comp-R_pat)
        xlabel('\theta (degrees)')
        ylabel('x''')
        % fname2=fullfile(export_path, sprintf( '%s_radonTransform%s.png'...
        %                 ,phase_nameFirst, char(h(oriNum))));        
        plotTitle = [char(h(oriNum)),sprintf('\n%s',...
                        phase_nameFirstPlot)];
        colormap(differenceColormap)
        title(plotTitle);
        xlabel('\theta (degrees)')
        ylabel('x''')
        % colormap gray
        

        exportgraphics(gcf,fname2)

        colorbar
         fname2=fullfile(export_path, sprintf( '%s_radonTransform%s_difference_cbar.png'...
                        ,phase_name, char(h(oriNum))));
        exportgraphics(gcf,fname2)
        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s_difference_gray.png'...
                        ,phase_name, char(h(oriNum))));

        imwrite(normaliseForImwrite(diff_pat),fname2);
        % % % % %% ------------  difference Radon Transform-----------------
        % % % % h = emptyHs{phaseNumber}
        % % % % diff_pat = comparisonPatterns(:,:,oriNum) - ...
        % % % %     patterns(:,:,oriNum);
        % % % % fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference.png'...
        % % % %                 ,phase_name, char(h(oriNum))));
        % % % % plotTitle = [char(h(oriNum)),sprintf('\n%s \n- %s',...
        % % % %                 phase_nameFirstPlot,phase_namePlot)];
        % % % % plotComparison(diff_pat,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
        % % % % 
        % % % % % export gray value image
        % % % % fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_gray.png'...
        % % % %                 ,phase_name, char(h(oriNum))));
        % % % % imwrite(normaliseForImwrite(diff_pat),fname2);
        

        %% label peaks from radon diff on the pat
        %%
        pat_in = EBSP_BGCor(comparisonPatterns(:,:,oriNum),...
                    rtmConfigRefineRadon.Settings_Cor)-...
                    EBSP_BGCor(patterns(:,:,oriNum),rtmConfigRefineRadon.Settings_Cor);

        [ Peak_Centre,Peak_Quality,Peak_Set,R_EBSP,R_Edge,R_RHO,R_theta ] = EBSP_RadHunt( pat_in,rtmConfigRefineRadon.Settings_Rad );
        % r_ratio = size(R_EBSP,2)/size(pat_in,1);

        figure
        imagesc(R_theta(:),R_RHO(:),R_EBSP);
        colormap gray
        hold on 
        scatter(Peak_Centre(:,1),Peak_Centre(:,2),100,"s",'MarkerEdgeColor',[0 1  1],'LineWidth',2)
        for i=1:size(Peak_Centre)
            text(Peak_Centre(i,1)-7,Peak_Centre(i,2)-10,sprintf("%i",i),'Color',[0 1 1])
        end
        % exportgraphics('radon_Peaks')

        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s_peaks_labelled.png'...
                        ,phase_name, char(h(oriNum))));
        exportgraphics(gcf,fname2)
        %
        [x,y] = pol2cart(deg2rad(Peak_Centre(:,1)-90),500);
      % [xp,yp] = pol2cart(Peak_Centre(:,1),50);
        % [x_text,y_text] = pol2cart(deg2rad(Peak_Centre(:,1)-90),500);
        
        figure
            
        imagesc(flipud(pat_in))
        colormap(differenceColormap)
        % axis xy
        hold on 
        centerX = ceil(size(pat_in,2)/2);
        centerY = ceil(size(pat_in,1)/2);
    
        offset_x  = Peak_Centre(:,2).* sin(deg2rad(Peak_Centre(:,1)-90))
        offset_y = Peak_Centre(:,2).* cos(deg2rad(Peak_Centre(:,1)-90))
    
        for i=1:size(Peak_Centre,1)
            plot([centerX-x(i)+offset_x(i) centerX+x(i)+offset_x(i)], ...
            [centerY+y(i)+offset_y(i) centerY-y(i)+offset_y(i)],... %+y(i)+offset_y(i) centerY-y(i)-offset_y(i)], ...
            "k",LineWidth=2)
            % scatter([centerX-x(i)+offset_x(i) centerX+x(i) + offset_x(i)], ...
            % [centerY-y(i)+offset_y(i) centerY+y(i)+offset_y(i)])
            text(centerX+offset_x(i)+3,centerY+offset_y(i)+3,sprintf("%i",i),'Color','k','fontSize',18)
        end
        % for i=1:size(Peak_Centre)
        %         text(Peak_Centre(i,1)-7,Peak_Centre(i,2)-10,sprintf("%i",i),'Color',[0 1 1])
        % end
        fname2=fullfile(export_path, sprintf( '%s_radonTransform%s_pattern_bands_labelled.png'...
                        ,phase_name, char(h(oriNum))));
        exportgraphics(gcf,fname2)
        %%
        % figure
        % imagesc(R_EBSP)
        % colormap gray
        %% -------------- difference pos and neg separate ----------------



        plotTitle =[char(h(oriNum)),sprintf('%s \n- %s',phase_nameFirstPlot,phase_namePlot)];

        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_minus.png'...
                    ,phase_name, char(h(oriNum))));

        plotComparison(diff_pat,cmap_reds,differenceLimitsNeg,plotTitle,fname2,rtmConfigRefine.Settings_Cor);

        


        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_plus.png' ...
                    ,phase_name, char(h(oriNum))));     

        plotTitle = [char(h(oriNum)),sprintf('\n%s \n- %s',phase_nameFirstPlot,phase_name)];

        plotComparison(diff_pat,cmap_blues,differenceLimitsPos,plotTitle,fname2,rtmConfigRefine.Settings_Cor);

        %% ---------- calculate similarity metrics

        [ssimAll,ssimAllMap] = ssim(comparisonPatterns(:,:,oriNum), ...
                    patterns(:,:,oriNum));

        [luminance,luminanceMap] = ssim(comparisonPatterns(:,:,oriNum),...
                    patterns(:,:,oriNum),'exponents',[1 0 0]);

        [contrast, contrastMap] = ssim(comparisonPatterns(:,:,oriNum),...
                    patterns(:,:,oriNum),'exponents',[0 1 0] ) ;

        [structureVal, ssimstructureMap] = ssim(comparisonPatterns(:,:,oriNum),...
                    patterns(:,:,oriNum),'exponents', [0 0 1] ) ;

        export_path = fullfile(RootExport_path,'SSIM');
        mkdir(export_path);

        SSIMValues(:,oriNum) = ssimAll;
        SSIMLuminanceValues(:,oriNum) = luminance;
        SSIMStructureValues(:,oriNum) = structureVal;
        SSIMContrastValues(:,oriNum) = contrast  ;
        %% --------- Plot SSIM mappings -------------------------


        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimAll.png' ,...
                phase_name, char(h(oriNum))));  

        plotTitle = [char(h(oriNum)),sprintf('\nSSIM : %.4f \n %s \n %s',...
                    ssimAll,phase_nameFirstPlot,phase_name)];

        plotComparison(ssimAllMap,ssim_colormap,ssimLimits,plotTitle,...
                    fname2,rtmConfigRefine.Settings_Cor);
        

        %luminance

        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimLuminance.png' ...
                    ,phase_name, char(h(oriNum))));

        plotTitle = [char(h(oriNum)),sprintf('\nSSIM Luminance : %.4f \n %s \n %s',...
                    luminance,phase_nameFirstPlot,phase_name)];

        plotComparison(luminanceMap,ssim_colormap,ssimLimits,plotTitle,...
                    fname2,rtmConfigRefine.Settings_Cor);
        


        % contrast
      

        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimContrast.png' ...
                    ,phase_name, char(h(oriNum))));

        plotTitle = [char(h(oriNum)),sprintf('\nSSIM Contrast : %.4f \n %s \n %s', ...
                    contrast,phase_nameFirstPlot,phase_name)];

        plotComparison(contrastMap,ssim_colormap,ssimLimits,plotTitle,fname2,...
                    rtmConfigRefine.Settings_Cor);




        % structure
        

        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimStructure.png' ...
                    ,phase_name, char(h(oriNum))));

        plotTitle = [char(h(oriNum)),sprintf('\nSSIM Structure : %.4f \n %s \n %s', ...
                    structureVal,phase_nameFirstPlot,phase_name)];

        plotComparison(ssimstructureMap,ssim_colormap,ssimLimits,plotTitle,fname2,...
                    rtmConfigRefine.Settings_Cor);
    writeDataToH5(InputUser.h5name,sprintf("/%s/NCCValues",phase_name),NCCValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMContrastValues",phase_name),SSIMContrastValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMLuminanceValue",phase_name),SSIMLuminanceValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMStructureValue",phase_name),SSIMStructureValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMValues",phase_name),SSIMValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/DataOrder",phase_name),"The Values should be 1d order. each row corresponding to a unique orientation. The Values are references against phase 1")
    writeDataToH5(InputUser.h5name,sprintf("/%s/referencePhase",phase_name),phase_nameFirstPlot)

    writeDataToH5(InputUser.h5name,sprintf("/%s/comparePhase",phase_nameFirst),phase_nameFirst,"String")

          close all      
    end

end

% 
% for pNumber = 1:numel(InputUser.Phase_Input)
%     %% --------- write more data to h5 struct ------------------
%     phase_name = InputUser.Phase_Input{pNumber};
% 
% 
%     close all
% end
function plotComparison(plotImage,cmap_in,colorlims,plotTitle,fname,SettingsCor)
        f=figure;
        f.Position=[100 100 1000 1000];
        imagesc(plotImage);      
        title(plotTitle);
        colormap(cmap_in)
        axis('equal','xy')
        xlim([0,SettingsCor.size])
        ylim([0,SettingsCor.size])
        clim(colorlims)
        colorbar()
        exportgraphics(gcf,fname)
end

%%
%     figure
%     f=figure
%     h_ecp = axes('Units','normalized','Position',[0 0 1 1],'Parent',f);
%     % h_ecp=imagesc(pattern)
%     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern,'Parent',h_ecp)
%     axis(h_ecp,'equal','xy','tight');
%     hold(h_ecp,'on');
%     colormap gray;
%      cset=InputData.cset;
%       sf=sf_calc(EBSD_geometry);
%
%      % cset=InputData.cset;
% % if h_ecp_bands.Value == 1
% hold on
%     if ~isempty(Family_List)
%         for n=1:size(Family_List,1)
%             HKLs=Family_List{n}(:,2:4);
%
%             for p=1:size(HKLs,1)
%                 %calculate the bands
%                 [bands] = Cone_Build(HKLs(p,:),Crystal_UCell{1}.Astar,phaseNumbernputData.lambda_p,r,cone_params);
%                 %plot them on the EBSP/ECP
%                 %%
%                 hold on
%                 if p==1
%                     pBand(bands,cset(n,:),h_ecp,sf,'display',HKLs(p,:))
%                 else
%                     pBand(bands,cset(n,:),h_ecp,sf)
%                 end
%             end
%         end
%     end
% % % end
% h_ecp.XLim=([EBSD_geometry.x_screen(1),EBSD_geometry.x_screen(end)]);
% h_ecp.YLim=([EBSD_geometry.y_screen(1),EBSD_geometry.y_screen(end)]);
% h_ecp.FontSize=8;
% legend()