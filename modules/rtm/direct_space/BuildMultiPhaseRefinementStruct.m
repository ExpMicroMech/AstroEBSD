function MultiPhaseStruct = BuildMultiPhaseRefinementStruct(InputUser,rtmConfig)
    Phase_Input = InputUser.Phase_Input;
   MultiPhaseStruct = struct();
    for phaseNumber = 1:size(Phase_Input,2)
        % try:

        [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {Phase_Input{phaseNumber}},rtmConfig.RTM.Phase_Folder);


        cs_phase = loadCIF_atoms(RTM_info.cif_file);


        if rtmConfig.generalOptions.AlignEdax == 1


            cs_phase =  crystalAlignEdax(cs_phase);


        end

   % ) disp(RTM_info
    disp("Debug:")
    %...
        disp(phaseNumber)
    
%%
        % disp(RTM_info.bin_file)
        [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
        %%
        % InputUser.Phase_Input
        % rtmConfig.Settings_Cor

        MultiPhaseStruct(phaseNumber).screen_int = screen_int;
        MultiPhaseStruct(phaseNumber).cs_phase = cs_phase;
        try
        [EBSD_geometry] = EBSP_Gnom(rtmConfig.Settings_Cor,rtmConfig.Settings_PCin.start);
        MultiPhaseStruct(phaseNumber).EBSD_geometry = EBSD_geometry;
        catch
            warning('EBSD geometry not build')
        end
        % %% test the interpolant      
        % EBSD_geometry=EBSP_Gnom(PatternInfo,Settings_PCin.start);
        % rot=rotation.byEuler(0,0,0);
        % r=rot.matrix;
        % euler=rot.Euler();
        % pattern_simulated=EBSP_gen(EBSD_geometry,r,screen_int,0);
        % figure
        % imagesc(pattern_simulated)
        % axis xy;
        % axis equal;
        % colormap gray;
        % catch
            % warning(1)
    end
    try
        writeStructToH5(MultiPhaseStruct, InputUser.h5name,"/MultiPhaseStruct")
    catch
     warning("Could not write MultiPhaseStruct to h5")
    end
end