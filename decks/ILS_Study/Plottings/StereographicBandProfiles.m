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
subFolder = 'YNi2_dbwCompare'
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))
InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead)
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end
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
root_path=fullfile("I:\sciebo\EBSD\YNi2\Bandprofiles",subFolder)
% root_path = fullfile('I:\002Simulation_EBSD\NbNiBandProfiles',subFolder)
% root_path='E:\RWTH_Aachen\Promo\CIF\ecp'
InputUser.ExportPath = root_path
%% save some stuff to structured data
InputUser.h5name = fullfile(root_path,sprintf('StereogrpahicBandProfiles_manBandNames_2.hdf5'))

try
    makeNewH5(InputUser.h5name)
catch
    warning("File could not be created, mak sure, you have setup git properly on your machine")
end



%% ---------------------------- Set up Zone Axes to put in view-------------
% h=Miller({0,0,1},{1,1,1},{0,1,1},cs);
% h=Miller('UVTW',cs)
hList = {{0 0 0 1}',{-1,0,1,0}',{-1,2,-1,0}',{1, 1, -2, 1}',{3, 0, -3, 1}',}
hList = {{0 0 4}',{4 0 4}', {4 4 4}', {2 2 6}', {2 2 4}',{6 6 10}'}


hType = 'UVW';
hnames = strings(numel(hList),numel(InputUser.Phase_Input))
InputUser.mPlotNames = hList

%%
emptyHs = cell(numel(InputUser.Phase_Input),1)
% h is for the zone axis along which we want to observe the ecpori=[]


numOris = numel(hList)


for i = 1:numel(InputUser.Phase_Input)
     InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % load ebsd map
    rtmConfig = defaultRTMConfig()
    
    cs = crystalAlignEdax(cs);
    % for the emsoft master patterns we do not want the crystal Align Edax
    InputUser.CSList{i}=cs;
end

%% -------------- try to setup miller indices externally---------------------------
RTM_info.Voltage = 19.5E3;
structConfig = defaultStructureFactorConfig(RTM_info);
% structConfig.Phase_Folder = InputUser.Phase_Folder;
structConfig.maxInd = 10;
% structConfig = defaultStructureFactorConfig(RTM_info);
structConfig.Phase_Folder = InputUser.Phase_Folder;
structConfig.Phase_Folder=Phase_Folder   
structConfig.lambda = calc_lambda(structConfig.V_in);
structConfig.lambda_p=structConfig.lambda/structConfig.lattice_param_scale;
[m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{1}}, structConfig.lambda_p ,structConfig);
[m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
if ~isfield(InputUser,'numHKLPlot')
    InputUser.numHKLPLot=10;
end
[m_plot,I] =cleanFHKLDuplicates(m_plot,I,InputUser);
%%
structConfig.IPlot = I;
InputUser.mPlotNames = {}
for i = 1:numel(m_plot)
    InputUser.mPlotNames{i} = {m_plot(i).h m_plot(i).k m_plot(i).i m_plot(i).l}';
end
InputUser.mPlotNames = hList
%%
InputUser.htype = 'HKL';
% Miller(InputUser.mPlotNames{:},InputUser.htype,m_plot.CS);
 % {{0 0 0 1}',{-1,0,1,0}',{-1,2,-1,0}',{1, 1, -2, 1}',{3, 0, -3, 1}',}

%% Set up plotting resolution and parameters ------------------------------
% [InputData.lambda] = calc_lambda(InputData.V_in);

%scale the wavelength to the same units as the dhkl
% InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale; 
%lambda normalised by the scale factor of the unit cell
patterns=zeros(2048,2048,numOris); 
% foruth number for size of millers to plot
rtmConfigRefine = defaultRTMConfig();
rtmConfigRefine.Settings_Cor.screensize = 2000; %hardcoded for plotting
rtmConfigRefine.Settings_Cor.size = 2000;
% PC= [0.5 0.25 0.7]
% rrmtConfigRefine.PC = PC;



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

    structConfig = defaultStructureFactorConfig(RTM_info);
    structConfig.maxInd = 10;
    structConfig.Phase_Folder = InputUser.Phase_Folder;
    [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
    [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
    [m_plot,I] =cleanFHKLDuplicates(m_plot,I);
    structConfig.IPlot = I;
% end
    [cs,~]=loadCIF(RTM_info.cif_file);
%%
    [stereoPattern] = Cube_GenerateStereo(RTM_info.bin_file,RTM_info.isHex);
    %%
    %
    try
        RootExport_path=fullfile(root_path,phase_name,sprintf("%.2fkV",RTM_info.Voltage));
    catch
        RootExport_path=fullfile(root_path,phase_name,"defalutkV");
    end
    mkdir(RootExport_path);
    % 

    % SampleBandProfilesFromReprojectionsStereo(InputUser,structConfig)

    
   %% --------- save patterns to structured h5 ----------------
    try 
        voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
    catch
        voltageString = "defaultVoltage"
    end
    %%
    h5PhaseInfo=    Phase_BuilderSimplified(phase_name,Phase_Folder)
    h5PhaseInfo.phaName =string(phase_name);
    writeStructToH5(h5PhaseInfo,InputUser.h5name,sprintf("/%s/%s/SimulationInfo",voltageString,phase_name))
    % writeStructToH5(h5PhaseInfo,stereoPattern,sprintf("/%s/%s/StereoPatternRepro",voltageString,h5PhaseInfo.Name))
    writeDataToH5(InputUser.h5name,sprintf("/%s/%s/StereoPatternRepro",voltageString,phase_name),...
        squeeze(stereoPattern));
    % writeDataToH5(InputUser.h5name,sprintf("/%s/%s/ZoneAxes",voltageString,h5PhaseInfo.Name),hnames,"string")
    % % writeDataToH5(InputUser.h5name,sprintf("/%s/Patterns",phase_name),...
    %     % squeeze(patterns(:,:,:)),"single");
    % hnames(hnum,phaseNumber) = string(char(h(hnum)));

     structConfig.voltageString = voltageString
     structConfig.h5PhaseInfo = h5PhaseInfo 
     %%

%%
 % h5PhaseInfo=    Phase_BuilderSimplified(phase_name,Phase_Folder)
InputUser.h5nameRepro = InputUser.h5name
% structConfig.voltageString = "defaultVoltage"
 % structConfig.h5PhaseInfo = h5PhaseInfo 
%% ----------------------- Extract Band Profiles ---------------------------
structConfig.Phase_Folder = InputUser.Phase_Folder;
structConfig.theta=(-180:0.1:180)*degree;
structConfig.BandWidth = 5*degree;
InputUser.BandWidth = 5;
structConfig.omegaResolution = 0.05; %in degree
structConfig.binSize = 0.5 %in degree;
InputUser.AlignEMsoft = 1;
structConfig.h5PhaseInfo.Name = phase_name
% InputUser.mPlotNames = hList
SampleBandProfilesFromReprojectionsStereo(InputUser,structConfig)

end
%% ------------------------ EXtract perpendicular Band Profiles ----------------
% FHKLoptions.Voltage = 20E3;
% structConfig = defaultStructureFactorConfig(FHKLoptions);     
% phaseNumber = 1
structConfig.Phase_Folder = InputUser.Phase_Folder;
structConfig.theta=(-180:0.1:180)*degree;
structConfig.BandWidth = 5*degree;
structConfig.omegaResolution = 0.05; %in degree
structConfig.binSize = 0.5 %in degree;
InputUser.AlignEMsoft = 1;
structConfig.h5PhaseInfo.Name = phase_name
%%
SampleBandProfilesFromReprojectionsStereoPerp(InputUser,structConfig)
% a



%% --------------- compare differences --------------------------
comparisonPhase = 1; % aVoid magic numers :D





% set up some more colormap and plotting variables
% limits
differenceLimits = [-0.5,0.5];
differenceLimitsNeg = [-0.5,0];
differenceLimitsPos = [0,0.5]
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

NCCValues = zeros(numel(hList));
SSIMValues = zeros(numel(hList));
SSIMLuminanceValues = zeros(numel(hList));
SSIMStructureValues = zeros(numel(hList));
SSIMContrastValues  = zeros(numel(hList));
phase_nameFirst = InputUser.Phase_Input{comparisonPhase};
phase_nameFirstPlot = cleanNameLatex(phase_nameFirst);


   [ ~,~,~,~,~, RTM_info ] = ...
            Phase_Builder_RTM( {phase_nameFirst},InputUser.Phase_Folder );
 try 
        voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
    catch
        voltageString = "defaultVoltage"
 end
 h5PhaseInfo=    Phase_BuilderSimplified(phase_nameFirst,Phase_Folder)
% phase_name = InputUser.Phase_Input{phaseNumber}
comparisonPatterns =double(h5read(InputUser.h5name,sprintf("/%s/%s/StereoPatternRepro", ...
                    voltageString,phase_nameFirst)));
%% Plot the differences 
for phaseNumber=2:numel(InputUser.Phase_Input)
   
    
    
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
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)

    patterns = double(h5read(InputUser.h5name,sprintf("/%s/%s/StereoPatternRepro", ...
                    voltageString,phase_name)));
    % for oriNum= 1 : numOris



        

   
        export_path = fullfile(RootExport_path,'DifferenceStereographic');
        mkdir(export_path);
        
        ncc = normxcorr2(comparisonPatterns,...
                    patterns());

        NCCValues = max(ncc(:));
        %% ------------  difference-----------------
        h = emptyHs{phaseNumber}
        diff_pat = comparisonPatterns - ...
            patterns;
        fname2=fullfile(export_path, sprintf( '%s_stereographic_difference.png'...
                        ,phase_name ));
        plotTitle = [sprintf('\n%s \n- %s',...
                        phase_nameFirstPlot,phase_namePlot)];
        plotComparison(diff_pat,differenceColormap,differenceLimits,plotTitle,fname2,rtmConfigRefine.Settings_Cor);
        
        % export gray value image
        fname2=fullfile(export_path, sprintf( '%s_stereographic_difference_gray.png'...
                        ,phase_name));
        imwrite(normaliseForImwrite(diff_pat),fname2);


        %% -------------- difference pos and neg separate ----------------



        plotTitle =[sprintf('%s \n- %s',phase_nameFirstPlot,phase_namePlot)];

        fname2=fullfile(export_path, sprintf( '%s_stereographic_difference_minus.png'...
                    ,phase_name));

        plotComparison(diff_pat,cmap_reds,differenceLimitsNeg,plotTitle,fname2,rtmConfigRefine.Settings_Cor);

        


        fname2=fullfile(export_path, sprintf( '%s_stereographic_difference_plus.png' ...
                    ,phase_name));     

        plotTitle = [sprintf('\n%s \n- %s',phase_nameFirstPlot,phase_name)];

        plotComparison(diff_pat,cmap_blues,differenceLimitsPos,plotTitle,fname2,rtmConfigRefine.Settings_Cor);

        %% ---------- calculate similarity metrics

        [ssimAll,ssimAllMap] = ssim(comparisonPatterns, ...
                    patterns);

        [luminance,luminanceMap] = ssim(comparisonPatterns,...
                    patterns,'exponents',[1 0 0]);

        [contrast, contrastMap] = ssim(comparisonPatterns,...
                    patterns,'exponents',[0 1 0] ) ;

        [structureVal, ssimstructureMap] = ssim(comparisonPatterns,...
                    patterns,'exponents', [0 0 1] ) ;

        export_path = fullfile(RootExport_path,'SSIM');
        mkdir(export_path);

        SSIMValues= ssimAll;
        SSIMLuminanceValues = luminance;
        SSIMStructureValues= structureVal;
        SSIMContrastValues = contrast  ;
        %% --------- Plot SSIM mappings -------------------------


        fname2=fullfile(export_path, sprintf( '%s_stereographic_ssimAll.png' ,...
                phase_name));  

        plotTitle = [sprintf('\nSSIM : %.4f \n %s \n %s',...
                    ssimAll,phase_nameFirstPlot,phase_name)];

        plotComparison(ssimAllMap,ssim_colormap,ssimLimits,plotTitle,...
                    fname2,rtmConfigRefine.Settings_Cor);
        

        %luminance

        fname2=fullfile(export_path, sprintf( '%s_stereographic_ssimLuminance.png' ...
                    ,phase_name));

        plotTitle = [sprintf('\nSSIM Luminance : %.4f \n %s \n %s',...
                    luminance,phase_nameFirstPlot,phase_name)];

        plotComparison(luminanceMap,ssim_colormap,ssimLimits,plotTitle,...
                    fname2,rtmConfigRefine.Settings_Cor);
        


        % contrast
      

        fname2=fullfile(export_path, sprintf( '%s_stereographic_ssimContrast.png' ...
                    ,phase_name));

        plotTitle = [sprintf('\nSSIM Contrast : %.4f \n %s \n %s', ...
                    contrast,phase_nameFirstPlot,phase_name)];

        plotComparison(contrastMap,ssim_colormap,ssimLimits,plotTitle,fname2,...
                    rtmConfigRefine.Settings_Cor);




        % structure
        

        fname2=fullfile(export_path, sprintf( '%s_stereographic_ssimStructure.png' ...
                    ,phase_name));

        plotTitle = [sprintf('\nSSIM Structure : %.4f \n %s \n %s', ...
                    structureVal,phase_nameFirstPlot,phase_name)];

        plotComparison(ssimstructureMap,ssim_colormap,ssimLimits,plotTitle,fname2,...
                    rtmConfigRefine.Settings_Cor);
    writeDataToH5(InputUser.h5name,sprintf("/%s/NCCValues",phase_name),NCCValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMContrastValues",phase_name),SSIMContrastValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMLuminanceValue",phase_name),SSIMLuminanceValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMStructureValue",phase_name),SSIMStructureValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/SSIMValues",phase_name),SSIMValues)
    writeDataToH5(InputUser.h5name,sprintf("/%s/comparePhase",phase_nameFirst),phase_nameFirst,"String")

          close all      
    % end

end
%% Plot the bands flattened
% Plot the differences  
% H5InfoStruc = readS
mSize = 5;
for mNumber = 1:mSize
    %%
    h = tiledlayout(1,numel(InputUser.Phase_Input),'TileSpacing','none')
        f=figure;
    f.Position =[0 0 1900 240*numel(InputUser.Phase_Input)];

    for phaseNumber=1:numel(InputUser.Phase_Input)
   
    
    
    phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name,'ProfilesAlongBand');
    mkdir(RootExport_path);
    %% load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)

    bandImage = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandImage", ...
                    voltageString,phase_name,mNumber)));
    thetaPlot = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandParallelAngle", ...
                    voltageString,phase_name,mNumber)));
    omega = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/openingAngles", ...
                    voltageString,phase_name,mNumber)));
                    % voltageString,h5PhaseInfo.Name,mNumber)));
    bandName = h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandName", ...
                    voltageString,phase_name,mNumber));
    % for oriNum= 1 : numOris
    nexttile(phaseNumber)
    
%%
    % figure
    if max(thetaPlot/degree)>=179.99
        xlims = [90 180]
    else
        xlims = [0 90]
    end
    % figure
    
    
    ylims = [min(omega) max(omega)]
    % h = tiledlayout(6,1,"TileSpacing","none")
    % nexttile(1)
    imagesc(thetaPlot/degree,omega/degree,bandImage');
    axis equal;
    ylim(ylims/degree)
    xlim(xlims) % only half of it
    colormap 'gray'
    clim([-4 4])
    phaseNamePlotLoc = char(cleanNameLatex(InputUser.Phase_Input{phaseNumber})) %make it indexable

    maxLen=20
    for titleNum=1:ceil(size(phaseNamePlotLoc,2)/maxLen)
        try
            text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:(titleNum)*maxLen),'Units','normalized');
        catch
            text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:end),'Units','normalized');
        end
        % try
        %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
        % catch
        %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
        % end
    end
    % text(-0.14,1,phaseNamePlotLoc(20:30),'Units','normalized');
    text(-0.12,0.4,bandName,'Units','normalized');
    if ~(phaseNumber == numel(InputUser.Phase_Input))
        axis off
    end
    
    % nexttile(2);
    % imagesc(thetaPlot/degree,omega/degree,valsSim');
    % axis equal;
    % ylim([90-BandWidth/degree, 90+BandWidth/degree])
    % xlim(xlims) % only half of it
    % colormap 'gray'
    % clim([-4 4])
    % 
    % ax3 = nexttile(3);
    % imagesc(thetaPlot/degree,omega/degree,valsSim'-valsExp');
    % axis equal;
    % ylim([90-BandWidth/degree, 90+BandWidth/degree])
    % xlim(xlims) % only half of it
    % colormap(ax3,returnCmapDiv)
    % clim([-4 4])
    % 
    % 
    % nexttile([3 1])
    % 
    % plot(thetaPlot/degree,averageExp+offsets,'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
    % xlim(xlims)
    % xlabel("azimuth angle")
    % ylabel("band opening angle (°) ")
    % hold on
    % plot(thetaPlot/degree,averageSim+offsets,'Color',[0 1 1],'DisplayName','Simulation')
    % % legend()
    

        
    end
    exportgraphics(f,fullfile(RootExport_path,sprintf('BAndComparison_%s.png',bandName) ))
    close all
%%
end
%% Plot Gradients
mSize = 5;
for mNumber = 1:mSize
    %%
    h = tiledlayout(1,numel(InputUser.Phase_Input),'TileSpacing','none')
        f=figure;
    f.Position =[0 0 1900 240*numel(InputUser.Phase_Input)];

    for phaseNumber=1:numel(InputUser.Phase_Input)
   
    
    
    phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name,'ProfilesAlongBand');
    mkdir(RootExport_path);
    %% load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)

    bandImage = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandImage", ...
                    voltageString,phase_name,mNumber)));
    [bandImage,~] = gradient(bandImage);
    thetaPlot = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandParallelAngle", ...
                    voltageString,phase_name,mNumber)));
    omega = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/openingAngles", ...
                    voltageString,phase_name,mNumber)));
                    % voltageString,h5PhaseInfo.Name,mNumber)));
    bandName = h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandName", ...
                    voltageString,phase_name,mNumber));
    % for oriNum= 1 : numOris
    nexttile(phaseNumber)
    
%%
    % figure
    if max(thetaPlot/degree)>=179.99
        xlims = [90 180]
    else
        xlims = [0 90]
    end
    % figure
    
    
    ylims = [min(omega) max(omega)]
    % h = tiledlayout(6,1,"TileSpacing","none")
    % nexttile(1)
    imagesc(thetaPlot/degree,omega/degree,bandImage');
    axis equal;
    ylim(ylims/degree)
    xlim(xlims) % only half of it
    % colormap(returnCmapDiv)
    colormap 'gray'
    clim([-1 1])
    phaseNamePlotLoc = char(cleanNameLatex(h5PhaseInfo.Name)) %make it indexable

    maxLen=20
    for titleNum=1:ceil(size(phaseNamePlotLoc,2)/maxLen)
        try
            text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:(titleNum)*maxLen),'Units','normalized');
        catch
            text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:end),'Units','normalized');
        end
        % try
        %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
        % catch
        %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
        % end
    end
    % text(-0.14,1,phaseNamePlotLoc(20:30),'Units','normalized');
    text(-0.12,0.4,bandName,'Units','normalized');
    if ~(phaseNumber == numel(InputUser.Phase_Input))
        axis off
    end
    
    % nexttile(2);
    % imagesc(thetaPlot/degree,omega/degree,valsSim');
    % axis equal;
    % ylim([90-BandWidth/degree, 90+BandWidth/degree])
    % xlim(xlims) % only half of it
    % colormap 'gray'
    % clim([-4 4])
    % 
    % ax3 = nexttile(3);
    % imagesc(thetaPlot/degree,omega/degree,valsSim'-valsExp');
    % axis equal;
    % ylim([90-BandWidth/degree, 90+BandWidth/degree])
    % xlim(xlims) % only half of it
    % colormap(ax3,returnCmapDiv)
    % clim([-4 4])
    % 
    % 
    % nexttile([3 1])
    % 
    % plot(thetaPlot/degree,averageExp+offsets,'Color',[0.5 0.5 0.5],'DisplayName','Experiment')
    % xlim(xlims)
    % xlabel("azimuth angle")
    % ylabel("band opening angle (°) ")
    % hold on
    % plot(thetaPlot/degree,averageSim+offsets,'Color',[0 1 1],'DisplayName','Simulation')
    % % legend()
    

        
    end
    exportgraphics(f,fullfile(RootExport_path,sprintf('BAndComparison_Gradients%s.png',bandName) ))
    close all
%%
end



%% plot the profiles
%% 
colormapProfiles=cbrewer('seq','YlGn',numel(InputUser.Phase_Input)+2);
colormapProfiles =colormapProfiles(3:end,:)
colormapProfiles = rescale(colormapProfiles,0,0.8)
% cmap_reds=flipud(cmap_reds)
% colormapProfiles(colormapProfiles<0)=0;
% colormapProfiles(colormapProfiles>1)=1; 
mSize = 5;
phaseOffset = 0.03;
%
for mNumber = 1:mSize
    %%
    % h = tiledlayout(1,numel(InputUser.Phase_Input),'TileSpacing','none')
        f=figure;
    f.Position =[0 0 1900 240*numel(InputUser.Phase_Input)];
    hold on
    for phaseNumber=1:numel(InputUser.Phase_Input)
   
    
    
    phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name,'ProfilesAlongBand');
    mkdir(RootExport_path);
    % load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)

    bandValues = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandValues", ...
                    voltageString,phase_name,mNumber)));
    offsets = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/offsets", ...
                    voltageString,phase_name,mNumber)));
    thetaPlot = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandParallelAngle", ...
                    voltageString,phase_name,mNumber)));
    omega = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/openingAngles", ...
                    voltageString,phase_name,mNumber)));
                    % voltageString,h5PhaseInfo.Name,mNumber)));
    bandName = h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandName", ...
                    voltageString,phase_name,mNumber));
    % for oriNum= 1 : numOris
    % nexttile(phaseNumber)
    
%
    % figure
    if max(thetaPlot/degree)>=179.99
        xlims = [90 180]
    else
        xlims = [0 90]
    end
    % figure
   phaseNamePlotLoc = char(cleanNameLatex(h5PhaseInfo.Name)) %make it indexable 
    
    ylims = [min(omega) max(omega)]
    % h = tiledlayout(6,1,"TileSpacing","none")
    % nexttile(1)
    % imagesc(thetaPlot/degree,omega/degree,bandImage');
    plot(-1,0,'color',colormapProfiles(phaseNumber,:),'DisplayName',phase_namePlot);
    plot(thetaPlot/degree,bandValues+offsets-90+phaseOffset*phaseNumber,'color',colormapProfiles(phaseNumber,:),'HandleVisibility','off');
    % axis equal;
    % ylim(ylims/degree)
    xlim(xlims) % only half of it
    colormap 'gray'
    clim([-4 4])


    % maxLen=20
    % for titleNum=1:ceil(size(phaseNamePlotLoc,2)/maxLen)
    %     try
    %         text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:(titleNum)*maxLen),'Units','normalized');
    %     catch
    %         text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:end),'Units','normalized');
    %     end
    %     % try
    %     %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
    %     % catch
    %     %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
    %     % end
    % end
    % % text(-0.14,1,phaseNamePlotLoc(20:30),'Units','normalized');
    % text(-0.12,0.4,bandName,'Units','normalized');
    % if ~(phaseNumber == numel(InputUser.Phase_Input))
    %     axis off
    % end

        % title(bandName,phaseNamePlotLoc)
    end
    legend()
    exportgraphics(f,fullfile(RootExport_path,sprintf('BAndComparison_Profiles%s.png',bandName) ))
    close all
%
end





%% Plot the differences
for mNumber = 1:mSize
    %%
    phaseNumber=1
    % h = tiledlayout(1,numel(InputUser.Phase_Input),'TileSpacing','none')
        f=figure;
    f.Position =[0 0 1900 240*numel(InputUser.Phase_Input)];
    hold on
        phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name,'ProfilesAlongBand');
    mkdir(RootExport_path);
    % load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)




    bandValuesCompare = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandValues", ...
                    voltageString,h5PhaseInfo.Name,mNumber)));
    for phaseNumber=1:numel(InputUser.Phase_Input)
   
    
    
    phase_name = InputUser.Phase_Input{phaseNumber};
    phase_namePlot = cleanNameLatex(phase_name);
     % phase_name = InputUser.pahse_
    RootExport_path=fullfile(root_path,phase_name,'ProfilesAlongBand');
    mkdir(RootExport_path);
    % load Patterns
    [ ~,~,~,~,~, RTM_info ] = ...
    Phase_Builder_RTM( InputUser.Phase_Input(phaseNumber),InputUser.Phase_Folder );
     try 
            voltageString = sprintf("%.2fkV",RTM_info.Voltage/1e3)
        catch
            voltageString = "defaultVoltage"
     end
     h5PhaseInfo=    Phase_BuilderSimplified(InputUser.Phase_Input{phaseNumber},Phase_Folder)

    bandValues = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandValues", ...
                    voltageString,h5PhaseInfo.Name,mNumber)));
    offsets = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/offsets", ...
                    voltageString,h5PhaseInfo.Name,mNumber)));
    thetaPlot = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandParallelAngle", ...
                    voltageString,h5PhaseInfo.Name,mNumber)));
    omega = double(h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/openingAngles", ...
                    voltageString,h5PhaseInfo.Name,mNumber)));
                    % voltageString,h5PhaseInfo.Name,mNumber)));
    bandName = h5read(InputUser.h5name,sprintf("/%s/%s/BandprofilesPerp/%i/bandName", ...
                    voltageString,h5PhaseInfo.Name,mNumber));
    % for oriNum= 1 : numOris
    % nexttile(phaseNumber)
    
%
    % figure
    if max(thetaPlot/degree)>=179.99
        xlims = [90 180]
    else
        xlims = [0 90]
    end
    % figure
    
    
    ylims = [min(omega) max(omega)]
    % h = tiledlayout(6,1,"TileSpacing","none")
    % nexttile(1)
    % imagesc(thetaPlot/degree,omega/degree,bandImage');
    try
    plot(thetaPlot/degree,bandValuesCompare-bandValues+(offsets-90)/10,'color',colormapProfiles(phaseNumber,:),'HandleVisibility','off');
    plot(-1,0,'color',colormapProfiles(phaseNumber,:),'DisplayName',phaseNamePlotLoc);
    % plot(thetaPlot/degree,bandValues+offsets-90+phaseOffset*phaseNumber,'color',colormapProfiles(phaseNumber,:),'HandleVisibility','off');
    % axis equal;
    % ylim(ylims/degree)
    xlim(xlims) % only half of it
    colormap 'gray'
    clim([-4 4])
    phaseNamePlotLoc = char(cleanNameLatex(h5PhaseInfo.Name)) %make it indexable
    end
    % maxLen=20
    % for titleNum=1:ceil(size(phaseNamePlotLoc,2)/maxLen)
    %     try
    %         text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:(titleNum)*maxLen),'Units','normalized');
    %     catch
    %         text(-0.12,1-(titleNum-1)*0.1,phaseNamePlotLoc((titleNum-1)*maxLen+1:end),'Units','normalized');
    %     end
    %     % try
    %     %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
    %     % catch
    %     %     text(-0.14,1,phaseNamePlotLoc(10:20),'Units','normalized');
    %     % end
    % end
    % % text(-0.14,1,phaseNamePlotLoc(20:30),'Units','normalized');
    % text(-0.12,0.4,bandName,'Units','normalized');
    % if ~(phaseNumber == numel(InputUser.Phase_Input))
    %     axis off
    % end

        % title(bandName,phaseNamePlotLoc)
    end
    legend()
    exportgraphics(f,fullfile(RootExport_path,sprintf('BAndComparison_Differnece%s.png',bandName) ))
    close all
%
end


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

