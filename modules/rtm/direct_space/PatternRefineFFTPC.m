function PatternRefineFFT(InputUser,rtmConfigRefine,MultiPhaseStruct) 
    % PatternRefineFFT Refines patterns using Fast Fourier Transform (FFT) methods.
    % (see Foden, RTM Paper)
    %
    % This function processes input data to refine patterns based on specified
    % configurations and phase structures. It reads experimental patterns,
    % applies background corrections, refines the orientation of the patterns and the projection center,
    % and saves the refined results into an HDF5 file.
    %
    % Parameters:
    %   InputUser (struct): A structure containing user input data including 
    %                       details about the experimental setup and a
    %                       pattern list.
    %       - h5name: Name of the HDF5 file where results will be saved.
    %       - PatternList: List of patterns to process.
    %       - Phase_Input: Cell array of phase identifiers for processing.
    %
    %   rtmConfigRefine (struct): Configuration settings for refinement,
    %                              including correction settings and RTM parameters.
    %                               Can be generated from default
    %
    %   MultiPhaseStruct (struct): Structure containing information about multiple 
    %                              phases, including screen intensity and hexagonal 
    %                              configuration settings for each phase.
    %
    % Returns:
    %   None: The function saves refined data directly into an HDF5 file specified 
    %         in InputUser.h5name without returning any values.
    %
    % Example:
    %   PatternRefineFFT(userInputData, configSettings, multiPhaseInfo);
    %


[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( rtmConfigRefine.Settings_Cor.size, rtmConfigRefine.RTM.LPTsize );
Detector_tilt = rtmConfigRefine.RTM.Rx(InputUser.MicroscopeData.TotalTilt);
numPatterns = numel(InputUser.PatternList)

for phaseNumber = 1:size(InputUser.Phase_Input,2)
    %%
    gNew = zeros(3,3,numPatterns);
    Library_PHList = zeros(1,numPatterns);
    PC_best_List =  zeros(3, numPatterns);
    G_Out = zeros(size(gNew));

    G_SO3 = h5read(InputUser.h5name,sprintf('/RTM/%s/RTMCoarse/G_SO3',InputUser.Phase_Input{phaseNumber}));

    for patternNumber = 1: numPatterns
        
        Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
        %%
        [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,rtmConfigRefine.Settings_Cor);
        [Pat_Ref_BG] = EBSP_BGCor( Pat_exp_cor,rtmConfigRefine.Settings_Cor );
        [Pat_Ref_r] = refine_prep(Pat_Ref_BG,SettingsXCF,rtmConfigRefine.RTM);
        [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                        rtmConfigRefine.Settings_PCin.start );
        %refine
        [G_Refined,regout_R] = refine5(Pat_Ref_r, EBSD_geom,...
                        rtmConfigRefine.Settings_PCin.start, G_SO3(:,:,patternNumber),...
                        SettingsXCF, MultiPhaseStruct(phaseNumber).screen_int,...
                        MultiPhaseStruct(phaseNumber).screen_int.isHex,...
                        rtmConfigRefine.RTM);
        %% setup PC refinement
        rtmConfigRefine.RTM.isHex=MultiPhaseStruct(phaseNumber).screen_int.isHex;
             
        RTM_setup=rtmConfigRefine.RTM; %% at some point i got lost which info needs to be in which one
        RTM_setup.PatternFlip=0;
        RTM_setup.apply_gauss=0;
        %% fit pattern centre
        tic
        [G_Refined,PC_best] = fPCRefine( rtmConfigRefine.Settings_PCin.start,...
                Pat_Ref_BG, G_Refined, rtmConfigRefine.Settings_Cor,...
                rtmConfigRefine.Refine, SettingsXCF,...
                MultiPhaseStruct(phaseNumber).screen_int, rtmConfigRefine.RTM, RTM_setup);
        toc
        %% refine again orientaiton, to get optimal peak hight
         [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                        PC_best);
        [G_Refined,regout_R] = refine5(Pat_Ref_r, EBSD_geom,...
                        PC_best, G_Refined,...
                        SettingsXCF, MultiPhaseStruct(phaseNumber).screen_int,...
                        MultiPhaseStruct(phaseNumber).screen_int.isHex,...
                        rtmConfigRefine.RTM);
        gTemp = G_Refined*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
        gNew(:,:,patternNumber) = gTemp';

        G_Out(:,:,patternNumber)=G_Refined;
        Library_PHList(:,patternNumber) = regout_R(4);
        PC_best_List(:,patternNumber) = PC_best;
        %%
    end
    %%
    % astro to mtex eangs
    eangs_out = rotation.byMatrix(gNew);
    eangs_out_export = zeros(3,size(eangs_out,1));
    for i = 1:size(eangs_out,2)
       eangs_out_export(i,1) =     eangs_out(i).phi1;
       eangs_out_export(i,2) =     eangs_out(i).Phi;
       eangs_out_export(i,3) =     eangs_out(i).phi2;
    end

    %%
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/Refined_PH", ...
                InputUser.Phase_Input{phaseNumber}), Library_PHList);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/eangs_out", ...
                    InputUser.Phase_Input{phaseNumber}), eangs_out_export);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/G_Out", ...
                        InputUser.Phase_Input{phaseNumber}), G_Out);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/PC", ...
                        InputUser.Phase_Input{phaseNumber}), PC_best_List);
end
end