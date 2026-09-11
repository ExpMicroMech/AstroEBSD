function [PhaseParameters,FFTParameters,EBSD_geometry] = prepareRefinementMultiPhase(InputUser,RTM)

    [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );    
    %%
    [ EBSD_geometry] = EBSP_Gnom( InputUser.PatternInfo,InputUser.PC);
    if isfield(InputUser,'PatternInfoRTM')
        [ EBSD_geometry_RTM] = EBSP_Gnom( InputUser.PatternInfoRTM,InputUser.PC);
    else
        EBSD_geometry_RTM = EBSD_geometry;
    end
    for pNum = 1: size(InputUser.Phase_Input,2)
           %%
            % InputUser.Phase_Input = {'Nb6Ni7_em'}
            [ ~,~,~,~,~, RTM_info(pNum) ] = Phase_Builder_RTM({InputUser.Phase_Input{pNum}}, RTM.Phase_Folder);
            phases(pNum).csPhase=loadCIF_atoms(RTM_info(pNum).cif_file);
            [screen_int(pNum)] = Cube_Generate(RTM_info(pNum).bin_file,RTM_info(pNum).isHex);
            % 
            %populate SO(3)
            
            
            [ libraries(pNum).library_G ] = SO3_rotmat_gen(phases(pNum).csPhase,RTM.Sampling_Freq);
            %%
            [ template(pNum).template_library ] = Library_Gen(EBSD_geometry_RTM,screen_int(pNum),RTM_info(pNum).isHex, libraries(pNum).library_G ,2,SettingsXCF);
            radialAjdustments(pNum)=fPHFac(EBSD_geometry_RTM,InputUser.SettingsCorRTM,screen_int(pNum),RTM,SettingsXCF);
            %%
            [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
            cs = crystalAlignEdax(cs);
            % InputUser.CS = {'notIndexed',cs};
            InputUser.AlignEdax = 1;
            CSList(pNum) = cs;
    end

    [~,XCF_data_fill] = refine_prep(zeros(RTM.screensize,RTM.screensize),SettingsXCF,RTM);

    FFTParameters.radialAdjustments =radialAjdustments;
    FFTParameters.libraries = libraries;
    FFTParameters.template = template;
    FFTParameters.SettingsXCF = SettingsXCF;
    FFTParameters.SettingsXCF2 = SettingsXCF2;
    FFTParameters.XCF_data_fill = XCF_data_fill;
    FFTParameters.EBSD_geometry = EBSD_geometry_RTM;
    PhaseParameters.phase = phases;
    PhaseParameters.screen_int = screen_int;
    PhaseParameters.CSList = CSList;
    PhaseParameters.RTM_info = RTM_info;
end
