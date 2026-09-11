%% L. Berners (Beginning of 2025), This deck uses a combination of the fourier based 
% refined template matching and the real space based optimisation routine.
% You need to specify a phase and a rough estimate of the pattern centre,
% and the deck will try to refine the PC and OR at the same time.
% In this deck, we will utilise hdf5 files for export of the pattern
% matched Data.
%
%Furhter, we will try to outsurce the standard config.
% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;

InputUser.Phase_Folder=Phase_Folder
InputUser.Astro_loc=Astro_loc
%%
% load('I:\001Clara\PCA\51Nb-49Ni_V2_septemberbatch\_EBSD\vt0,2\extracted_patterns.mat')
% set up color coding
 [cmap_div]=cbrewer('div','RdBu',30);
cmap_div(cmap_div<0)=0;
%%

InputData.cset=[1,0,0; % Colors of the Bands in rgb (1 =255)
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0];
InputData.cset=[InputData.cset;InputData.cset*0.5];

InputData.V_in=20E3; %voltage in eV
InputData.lattice_param_scale=1E-10; %lattice paramters are provided in this unit, often Angstrom (1E-10);
i=1
%% Deal with the wavelegnth

[InputData.lambda] = calc_lambda(InputData.V_in);

%scale the wavelength to the same units as the dhkl
InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale;
% InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path
%%
rtmConfig = defaultRTMConfig();
RTM = rtmConfig.RTM
%%
%location of the binary files used for RTM
% InputUser.PatternLoc='J:\Users\Berners_Lukas\clara\Energies_Si_wafer\Si_Wafer_2024_12_09\Scan\corrected_THL_240.tif'
% InputUser.PatternLoc = 'J:\Users\Berners_Lukas\clara\Energies_Nb42NiV2_september\20kV10nA\9mm\100ms\maxingOut\THL400_0.3s_spot\corrected.tif'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 6__EBSD+EDS Map 2\_EBSD\vt0,2\Tile1\RotPCA_n29_w1,00e+00'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 2__EBSD+EDS Map 1\_EBSD\vt0,1\Tile1\RotPCA_n30_w1,00e+00'%
InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__58Nb-Ni_V2_september__Area 1__EBSD+EDS Map 2\_EBSD\vt0,2\Tile2\RotPCA_n21_w1,00e+00'%
InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\YNi2MaDinah\PCA\MA_Dinah_YNi2_S_D_repolish_Area 3_EBSD+EDS Map 3\_EBSD\vt0,1\Tile2\RotPCA_n10_w1,00e+00'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\MA_Dinah\Daten final\EBSD+EDX\YNi2_S_D\PCA\vt0,05\Tile1\RotPCA_n4_w1,00e+00'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\LaNi2\PCA\LaNi2_LaNi2 vibropol_Area 5_EBSD+EDS Map 1\_EBSD\vt0,1\Tile1\RotPCA_n18_w1,00e+00'
% InputUser.PatternFolder = 'J:\Users\Berners_Lukas\clara\LaNi2\PCA\LaNi2_LaNi2 vibropol_Area 5_EBSD+EDS Map 1\_EBSD\vt0,1\Tile1\RotPCA_n23_w1,00e+00'
% InputUser.PatternList=["_Pattern_1.png","_Pattern_9.png"]%,"_Pattern_3.png","_Pattern_4.png","_Pattern_5.png","_Pattern_9.png","_Pattern_14.png"]
InputUser.PatternList = ["_Pattern_1.png","_Pattern_2.png"]
% InputUser.PatternList = ["_Pattern_14.png"]%,"_Pattern_7.png","_Pattern_16.png"]
% InputUser.PatternList=["Pattern1.tif"]%;,"Pattern2.tif","Pattern3.tif","Pattern4.tif","Pattern5.tif","Pattern6.tif"]
% export_path='I:\001Clara\EnergyFiltering\NbNi\RTM'
% export_path = 'J:\Users\Berners_Lukas\clara\LaNi2\PCA\LaNi2_LaNi2 vibropol_Area 5_EBSD+EDS Map 1'
export_path = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__Nb-49Ni_V2_september__Area 2__EBSD+EDS Map 1\'%
export_path = 'J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\PCA\NbNi__58Nb-Ni_V2_september__Area 1__EBSD+EDS Map 2'
export_path = 'J:\Users\Berners_Lukas\clara\YNi2MaDinah\PCA\MA_Dinah_YNi2_S_D_repolish_Area 3_EBSD+EDS Map 3'





% Settings_PCin.start=[0.5 0.28 0.60];
% Settings_PCin.start = [0.52 0.4 0.6400];
Settings_PCin.start = [0.529145, 1-0.679621, 0.693581]
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );
%%
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%%
% [ Data_InputMap ]=PC_square( MapInfo.EBSPData, MapInfo.Data_InputMap,Settings_CorX );

%%

% InputUser.Phase_Input  = {'LaNi2_30kV','La7Ni16_30kV','La15Ni32_30kV'};
% InputUser.Phase_Input = {'Nb518Ni482_joubert_20.0kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.000_npx-1000',...
%                          'Nb518Ni482_joubert_20.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.800_npx-1000',...
%                         'Nb518Ni482_joubert_20.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.300_npx-1000',...
%                         'Nb518Ni482_joubert_19.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.800_npx-1000',...
%                         'Nb518Ni482_joubert_19.6kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.600_npx-1000',...
%                         'Nb518Ni482_joubert_19.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...}
%                         'Nb518Ni482_joubert_18.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.800_npx-1000',...}
%                         'Nb518Ni482_joubert_18.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.300_npx-1000',...
%                         'Nb518Ni482_joubert_17.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-17.800_npx-1000'}
% InputUser.Phase_Input = {'Nb569Ni431_joubert_20kV_2025_07_24__08_56_35_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.7_a4.9_3aNb_20kV_2025_07_22__17_06_04_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000','Nb518Ni482_joubert_20_2025_07_31__21_42_11_dmin-0.020_c1-4.000_c2-8.000_ehistmin-19.800_npx-1000','Nb518Ni482_joubert_20_2025_07_31__21_40_54_dmin-0.050_c1-6.000_c2-30.000_c3-100.000','Nb518Ni482_joubert_20_2025_07_31__17_55_50_dmin-0.020_c1-6.000_c2-20.000_ehistmin-19.800_npx-1000','Nb518Ni482_joubert_20_2025_07_31__21_39_41_dmin-0.050_c1-6.000_c2-20.000_ehistmin-19.800_npx-1000'}
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.7_a4.9_3aNb_20kV_2025_07_22__17_06_04_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.250_c26.7_a4.9_3aNi_20kV_2025_07_22__17_06_03_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.400_c26.7_a4.9','NbNi_dt0.300_c26.7_a4.9_3a1.00Ni_20kV_2025_07_22__17_06_03_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000','NbNi_dt0.300_c26.7_a4.9_3a0.00Ni_20kV_2025_07_22__17_05_59_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% 'Nb518Ni482_joubert_20_2025_07_31__11_19_18_dmin-0.010_c1-4.000_c2-12.000_ehistmin-19.800_npx-1000
% InputUser.Phase_Input = {'Nb518Ni482_joubert_20_2025_07_30__19_17_43_dmin-0.020_c1-6.000_c2-30.000_ehistmin-19.800_npx-1000'}
% InputUser.Phase_Input = {'YNi2'}
InputUser.maxInd = 30;
PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/Y*Ni*20kV*.pha'))'

InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead)
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.7_a4.9_19_2025_08_04__10_23_33_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.275_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.325_c26.7_a4.9_19_2025_08_04__10_23_35_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.300_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.350_c26.7_a4.9_19_2025_08_04__10_23_36_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.375_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.400_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     }
% f_j * exp( - 2*pis*i * (h*x(j) + k*y(j)+ l*z(j) ))
% InputUser
InputUser.cutoff = 0.05
InputUser.minD = 0.2
InputUser.AlignEdax =1 

 cutoff = 0.2       

extractions  = struct();
%%

for phaseNumber = 1:size(InputUser.Phase_Input,2)
    % RTM.screensize = 128; %size of the library patterns and the resize of the raw EBSPs
    % Settings_CorX.size=RTM.screensize; 
    % InputUser.Phase_Input = {'Nb6Ni7_em'}
    [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{phaseNumber}},RTM.Phase_Folder);
    cs_phase=loadCIF_atoms(RTM_info.cif_file);
    cs_phase =  crystalAlignEdax(cs_phase);
    [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
    [EBSD_geometry] = EBSP_Gnom( PatternInfo,Settings_PCin.start);
    %populate SO(3)
    
    % InputUser.Phase_Input  = {'Nb6Ni7_em'}; 
    % InputUser.Phase_Input  = {'Ta4Al5Fe3_vesta'}; 
    % [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{phsaeNumber}},RTM.Phase_Folder);
    % cs_phase=loadCIF(RTM_info.cif_file);
    % [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
    % [ EBSD_geometry] = EBSP_Gnom( PatternInfo,Settings_PCin.start);
    % %%%
    
    
    %% test the interpolant
    % EBSD_pat=
    
    EBSD_geometry=EBSP_Gnom(PatternInfo,Settings_PCin.start);
    rot=rotation.byEuler(0,0,0);
    r=rot.matrix;
    euler=rot.Euler();
    pattern_simulated=EBSP_gen(EBSD_geometry,r,screen_int,0);
    figure
    imagesc(pattern_simulated)
    axis xy;
    axis equal;
    colormap gray;
    
    %%
    [ library_G ] = SO3_rotmat_gen(cs_phase,RTM.Sampling_Freq);
    
    %%
    [ template_library ] = Library_Gen(EBSD_geometry,screen_int,RTM_info.isHex,library_G,2,SettingsXCF);
    %%
    
    radial_adjustment=fPHFac(EBSD_geometry,Settings_CorX,screen_int,RTM,SettingsXCF);

    for patternNumber = 1: size(InputUser.PatternList,2)
            %%
        RTM.screensize = 64; %size of the library patterns and the resize of the raw EBSPs
        Settings_CorX.size=RTM.screensize;
        InputUser.PatternLoc = fullfile(InputUser.PatternFolder,InputUser.PatternList(patternNumber))
        %     figure
        %     imshow(template_library(:,:,5),'Border','tight')
        %%
        % 
        InputUser.SampleTilt=0
        % InputUser.CameraTilt=5.7
        % MicroscopeData.SampleTilt=InputUser.SampleTilt
        % MicroscopeData.CameraTilt=InputUser.CameraTilt
        % MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180
        MicroscopeData.TotalTilt=0
        Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
        %%
        P_test=uint8(1)
        [~,XCF_data_fill] = refine_prep(zeros(Settings_CorX.size,Settings_CorX.size),SettingsXCF,RTM);
        % Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));
        Pat_exp =flipud(imread(InputUser.PatternLoc));
        if size(Pat_exp,3)>1
            Pat_exp = double(Pat_exp(:,:,1));
        end
        Pat_exp = double(Pat_exp)
        %%
        [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,Settings_CorX);
        Pat_Ref=Pat_exp_cor;
        % [Pat_Ref] = Pat_exp_cor;
        %                 %%
        [Pat_Ref_BG] = EBSP_BGCor( Pat_Ref,Settings_CorX );
        figure
        imagesc(Pat_Ref_BG)
        colormap gray
        axis xy
        axis equal
        pause(3)
        % error('end')
        % Full_PC(P_test,:)=[PCX(P_r,P_c),PCY(P_r,P_c),DD(P_r,P_c)];
        %%
        [Pat_Ref_r] = refine_prep(Pat_Ref_BG,SettingsXCF,RTM);
        %%
        % figure
        % imshow(Pat_Ref_r.FFT,'Border','tight')
        %%
        [G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r(P_test).FFT,SettingsXCF,XCF_data_fill,0);
        
        %%
        % [G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r(P_test).FFT,SettingsXCF,XCF_data_fill,0);
        %%
        % Pat_exp=flipud(double(imread(InputUser.PatternLoc)));
        Pat_exp =flipud(imread(InputUser.PatternLoc));
        Pat_exp = double(Pat_exp);
        if size(Pat_exp,3)>1
            Pat_exp = double(Pat_exp(:,:,1));
        end
        %%
        [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,Settings_CorX);
        
        %%
        num_lib=1
                        %generate the geometry for this location
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,Settings_PCin.start);
        %refine
        [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,Settings_PCin.start,G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);
        
        %store the orientation
        G_Out(:,:,P_test,num_lib)=G_Refined;
        %store the PH
        PH_Out(P_test,num_lib)=regout_R(4)./radial_adjustment;
        % send(D,P_test)
        %%
        % [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,pc_best);
        [ Pat_sim ] = EBSP_gen(EBSD_geom,G_Refined,screen_int);
        %%
        RTM.isHex=RTM_info.isHex;
        Refine.num_pts=6; %number of points to use for the refinement step
        Refine.ss=0.2; %initial probe volume
        Refine.p=2; %order of polynomial to fit to tetrahedron
        Refine.n_its=10; %number of interations
        Refine.debug=1; 
        
        RTM_info.screensize=PatternInfo.ScreenHeight
        RTM_setup=RTM %% at some point i got lost which info needs to be in which one
        % RTM_setup.InputMode=InputUser.Mode;
        % RTM_setup.issorted=InputUser.issorted;
        RTM_setup.PatternFlip=0;
        RTM_setup.apply_gauss=0;
        %%
        % [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );
        [G_Refined,PC_best] = fPCRefine(Settings_PCin.start,Pat_Ref_BG,G_Refined,PatternInfo,Refine,SettingsXCF,screen_int,RTM,RTM_setup)
        
        %%
        % G_Refined
        % figure
        % regout_R
        % [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,PC_best,G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);
        
        %%
        PatternInfo_new.size=[514, 514]
        [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,PC_best);
        [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_Refined,screen_int);
        
        figure
        imagesc(Pat_sim);
        axis equal;
        axis xy;
        %%
        eangs=conv_G_to_EA(G_Refined);
        RTM.screensize = 256; %size of the library patterns and the resize of the raw EBSPs
        Settings_CorX.size=RTM.screensize;
        Pat_exp=double(flipud(imread(InputUser.PatternLoc)));
        if size(Pat_exp,3)>1
            Pat_exp = double(Pat_exp(:,:,1));
        end
        [ Pat_Ref ] = EBSP_BGCor( Pat_exp,Settings_CorX);
        
        %%
        [pc_best,eangs_out,phbest,LPR]=optimize_pc_ori(PC_best,eangs,Pat_Ref,Settings_CorX,screen_int,'Mode','fminsearch','refinePC')
        %%
        [pc_best,eangs_out,phbest,LPR]=optimize_pc_ori(pc_best,eangs_out,Pat_Ref,Settings_CorX,screen_int,'Mode','fminsearch')
        g_refined=conv_EA_to_G(eangs_out)
        %%
        
        extractions.phase_index(patternNumber,phaseNumber)=phaseNumber;
        extractions.rotation_data_refined(:,:,patternNumber,phaseNumber)=g_refined;
        extractions.eangs_out(:,patternNumber,phaseNumber)=eangs_out;
        extractions.refined_pc(:,patternNumber,phaseNumber)=pc_best;
        extractions.peak_hight(:,patternNumber,phaseNumber)=phbest;
        extractions.rotation_data_refined_Fourier(:,:,patternNumber,phaseNumber)=G_Refined;
        extractions.refined_pc_fourrier(:,patternNumber,phaseNumber)=PC_best;
        % extractions.LPR(:,patternNumber,phaseNumber)=LPR;

        [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,pc_best);
        [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,g_refined,screen_int);
        %%
            % cross_size=Settings_CorX.cross;
            % center=uint16(size(EBSP2)/2);
            % % slice=[center(1)-cross_size:center(1)+cross_size,center(2)-cross_size,center(2)+cross_size]
            % % EBSP2(slice)=0;
            % % EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");
            % mask=ones(size(EBSP2));
            % mask(center(1)-cross_size+1:center(1)+cross_size,:)=0;
            % mask(:,center(2)-cross_size+1:center(2)+cross_size)=0;
            % % EBSP2(slice)=0;
            % EBSP2(~mask==0)=(EBSP2(~mask==0)-mean(EBSP2(~mask==0),"omitnan"))/std(EBSP2(~mask==0),"omitnan");
        %%
        % figure
        % imagesc(mask-flipud(mask))
        
        %%
        Settings_Cor_plot=Settings_CorX
        %
        Settings_Cor_plot.size=size(Pat_exp)
        Settings_Cor_plot.cross=2
        Settings_Cor_plot.quadrants=0
        [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,pc_best);
        % [ EBSD_geom_plot ] = EBSP_Gnom( PatternInfo_new,PC_best);
        [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,g_refined,screen_int);
        % [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_Refined,screen_int);
        Pat_sim=EBSP_BGCor(Pat_sim,Settings_Cor_plot);
        % [ Pat_sim_pseudo ] = EBSP_gen(EBSD_geom_plot,RTM.Rz(60*degree)*g_refined,screen_int);
        % Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));
        Pat_exp=double(flipud((imread(InputUser.PatternLoc))));
        if size(Pat_exp,3)>1
            Pat_exp = double(Pat_exp(:,:,1));
        end
        % figure
        Pat_exp_cor=EBSP_BGCor(Pat_exp,Settings_Cor_plot);
        % Pat_exp=Pat_exp-mean(Pat_exp(:));
        % Pat_exp=Pat_exp/std(Pat_exp(:));
        %%
        % imwrite(flipud(rescale(Pat_sim)),fullfile(export_path,['refined_pattern.tif']))
        
        [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, InputData.lambda_p ,InputUser);

%%

        [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,cutoff)
        m_plot = m_plot(1:10)
        if m_plot.CS.numSym>=24 %cubic :D
            I = I(1:5)
            m_plot = m_plot(1:5)
        else
            I = I(1:10)
        end
        InputData.IPlot = I;
        %%Pat_exp
        InputUser.Pattern = Pat_exp;
        InputUser.geometry = EBSD_geom_plot
        InputUser.gmatrix_refined = conv_EA_to_G(eangs_out) 

%%
        % [gridded_interpolant, InputUser] = PatternToInterpolant(InputUser,'ApplyMask',InputUser.PatternMask,'bgCor');
        InputUser.pattern_info = struct();
        InputUser.pattern_info.size = size( InputUser.Pattern );
        InputUser.PatternMask = ones(size(InputUser.Pattern));
        InputUser.PC_in = pc_best;
        colorbar()
        %%
        try
        InputUser= rmfield(InputUser,'PatternMask')
        end
        InputUser.cs = cs_phase;
        InputUserSim = InputUser;
        InputUserSim.Pattern = Pat_sim;
        [gridded_interpolant, InputUser] = PatternToInterpolant(InputUser,'gaussian',1.5,'bgCor');
        InputUserSim.Pattern = Pat_sim;

%
        [gridded_interpolantSim, InputUserSim] = PatternToInterpolant(InputUserSim,'gaussian',1.5,'bgCor');
      %
        figure
        imagesc(InputUser.Pattern_masked)
        axis xy
        axis equal
        clim([-3.3,3.3])
        colormap gray
%% 
        figure
        imagesc(InputUserSim.Pattern_masked)
        colorbar
        clim([-3.3,3.3])
           axis xy
        axis equal
        colormap gray
       %%
        img_out=reprojections_discrete(InputUser,gridded_interpolant,'img_size',2000,'fundamental');
        imwrite(img_out,fullfile(export_path,[sprintf('resampledPattern_exp%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
        tic
        img_outSim=reprojections_discrete(InputUserSim,gridded_interpolantSim,'img_size',2000,'fundamental');
        imwrite(img_outSim,fullfile(export_path,[sprintf('resampledPattern_Sim%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
toc
        %%
        img_outSim = img_outSim -mean(img_outSim(:));
        img_outSim = img_outSim/std(img_outSim(:));
        img_out = img_out -mean(img_out(:));
        img_out = img_out/std(img_out(:));
        figure
        imagesc(img_out)
        axis('equal','xy')
        colormap gray
        clim([-3 3])
        % %%
        % img_out = ( img_out-mean(img_out(:)) );
        % img_out = img_out / std( img_out(:));
        %%
        f = figure
        f.Position = [100 100 1100 1100]
        f.Units = 'pixels'
        imagesc(img_outSim-img_out)
        axis equal
        axis xy
        axis off
        colormap(cmap_div)
        %
        % figure
        % imagesc(diffImg)
        exportgraphics(gcf,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',600)
        % diffImg = img_outSim-img_out;%,cmap_div
        % max(diffImg(:))
        % min(diffImg(:))
        %  figure
        % imagesc(diffImg)
        % % exportgraphics(gcf,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',600)
        % rescale(diffImg)
        % %%
        % diffImg =rescale(diffImg)
        % %%
        % imwrite(diffImg,fullfile(export_path,[sprintf('resampledPattern_Difference%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]))
        %%
        % figure
        % imshow(diffImg)
        %% Plot the refinement results
     f=   figure
        % f.Position = [100 100 1500 600]
        f.Units = 'pixels'
        % figure
        subplot(2,3,1)
        imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim);
        axis equal;
        axis xy;
        axis off
        title('Simulated Pattern')
        clim([-3.5,3.5])
        ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        % figure
        % ax5=subplot(2,3,4)
        text(0,-0.07, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_out(1)/degree,eangs_out(2)/degree,eangs_out(3)/degree),'units','normalized','fontsize',14)
        % text(0, -0.14, sprintf('LPR= %.4f , %.4f, %.4f', LPR(1),LPR(2),LPR(3)),'units','normalized','fontsize',14)
        text(0, -0.21, sprintf('PC= %.4f , %.4f, %.4f', pc_best(1),pc_best(2),pc_best(3)),'units','normalized','fontsize',14)
        text(0, -0.28, sprintf('NCC= %.4f', phbest),'units','normalized','fontsize',14)
        % imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pat_sim_pseudo)
        
        % axis equal;
        % axis xy;
        % colormap(ax2,'gray');
        % ylim([EBSD_geometry.y_gn_min,EBSD_geometry.y_gn_max])
        % % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        % axis off
        % title('Pseudosymmetry')
        % clim([-3.5,3.5])
        % m = Miller({1,1,-2,0},cs_phase)
        % m=Miller({0,0,0,1},{2,0,-2,0},{2,2,-4,0},cs_phase)
        % for j=1:size(ori,2)
                % r=ori(j).matrix;
                % euler=ori(j).Euler();
                % r = 
        pattern =Pat_exp_cor;
        % fname=fullfile(export_path, sprintf( '%s_gnomonic_%.1f_%.1f_%.1f_%i.png' ,InputUser.Phase_Input{i}, euler(1),euler(2),euler(3),j));
        % imwrite(rescale(flipud(pattern)),fname)
        % fname = 'I:\001Clara\EnergyFiltering\Si\RTM\bandsOnPattern.png'
        % title(sprintf('Rz = %.0f °, Rx =%.0f°, Rz2 = %.0f°',angles(j),angles(i),0))
        % label_bands(pattern,EBSD_geom_plot,m,inv(rotation.byMatrix(g_refined)),InputData)
        % fname2=fullfile(export_path, sprintf( '%s_gnomonic_%.1f_%.1f_%.1f_%i_with_bands.png' ,InputUser.Phase_Input{i}, euler(1),euler(2),euler(3),j));
        % exportgraphics(gcf,fname2,'resolution',300)
        
                ax2=subplot(2,3,4)
        
        imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim)
        label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
        % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
        axis equal;
        axis xy;
        colormap(ax2,'gray');
        ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        axis off
        title('Simulated. Pattern','fontsize',14)
        clim([-3.5,3.5])
        % axis(app.UIAxes, 'tight')
        % set(gcf,'Resize','off')
        
        
        
        
        ax2=subplot(2,3,2)
        
        imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
        % label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
        axis equal;
        axis xy;
        colormap(ax2,'gray');
        ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        axis off
        title('Exp. Pattern','fontsize',14)
        clim([-3.5,3.5])
        % axis(app.UIAxes, 'tight')
        % set(gcf,'Resize','off')
        

        ax2=subplot(2,3,5)
        
        imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
        label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
        % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
        axis equal;
        axis xy;
        colormap(ax2,'gray');
        ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        axis off
        title('Exp. Pattern','fontsize',14)
        clim([-3.5,3.5])
        % axis(app.UIAxes, 'tight')
        % set(gcf,'Resize','off')
        
        
        ax=subplot(2,3,3)%subplot(1,3,3)
        % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
        imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,InputUser.Pattern_masked-InputUserSim.Pattern_masked)
        axis equal;
        axis xy;
        colormap gray;
        ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        clim([-3,3])
        % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        axis off
        title('Difference')
        colormap(ax,cmap_div)
        colorbar('position',[0.91,0.63,0.01,0.25])        
        
        % set(gcf, 'Position' , [100 000 1500 900])
        %%
        % ax=subplot(2,3,6)%subplot(1,3,3)
        % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,normalize(Pat_exp_cor)./normalize(Pat_sim))
        % axis equal;
        % axis xy;
        % colormap gray;
        % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
        % clim([0.01,10])
        % % caxis(log([0,10]));
        % % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
        % axis off
        % title('Division')
        % colormap(ax,cmap_div)
        % colorbar('position',[0.91,0.13,0.01,0.25])        
        % % set(ax6,'colorscale','log')
        set(gcf, 'Position' , [100 000 1500 900])
        f.Units = 'pixels'
        exportgraphics(gcf, fullfile(export_path,[sprintf('refined_pattern_Pattern%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
        % error('end')
        close all
        %%
        numBins = 3
binSize = 100/numBins
bins = prctile(Pat_sim(:),[0:binSize:100])
binsExp = prctile(Pat_exp_cor(:),[0:binSize:100])
% f=figure
f= figure
f.Position=[20,20,1900,1200]
h = tiledlayout(3,numBins, 'TileSpacing', 'none', 'Padding', 'none');
for i = 1:numBins
% subplot(5,3,i*3-2)
greyAx = nexttile(i)
% if bins(i)<0
    % setValue = 
plotImage = (imgaussfilt(Pat_sim,1.5));
plotImage(plotImage<bins(i))=0;
plotImage(plotImage>bins(i+1))=0;
imagesc(plotImage)
clim([bins(i),bins(i+1)])
axis equal
axis xy
axis off
colormap(greyAx,'gray')

% subplot(5,3,i*3-1)
greyAx= nexttile(i+numBins)
plotImageExp = (imgaussfilt(Pat_exp_cor,1.5));
plotImageExp = imgaussfilt(plotImageExp);
plotImageExp(plotImageExp<binsExp(i))=0;
plotImageExp(plotImageExp>binsExp(i+1))=0;
imagesc(plotImageExp)
clim([binsExp(i),binsExp(i+1)])
axis equal
axis off
axis xy
colormap(greyAx,'gray')
ax = nexttile(i+numBins*2)
% plotImageExp = (imgaussfilt(Pat_exp_cor,1.5));
% plotImageExp = imgaussfilt(plotImageExp);
% plotImageExp(plotImageExp<bins(i))=-1000;
% plotImageExp(plotImageExp>bins(i+1))=-1000;
% fusedImage = imfuse(plotImageExp,plotImage);
imagesc(plotImageExp-plotImage)
% imagesc(fusedImage)
clim([-3 3])
colormap(ax,cmap_div)
colorbar
axis equal
axis off
axis xy
end
f.Position=[20,20,1900,1200]
exportgraphics(gcf, fullfile(export_path,[sprintf('refined_pattern_Pattern_Thresholded%i_%s.tif',patternNumber,InputUser.Phase_Input{phaseNumber})]),'Resolution',300)
close all
    end
end
%%
%%% end
%%


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
% file_name=fullfile(export_path,'simulated_pattern.png')
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