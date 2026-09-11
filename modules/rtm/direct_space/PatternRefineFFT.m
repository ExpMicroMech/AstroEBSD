function PatternRefineFFT(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum) 
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

   if nargin <4
        tileNum = 1
    end

    %% options 
    if ~isfield(rtmConfigRefine,'options')
        rtmConfigRefine.options = struct();
    end
    if ~isfield(rtmConfigRefine.options,'TestPseudoSymmetry')
        rtmConfigRefine.options.TestPseudoSymmetry =0;
    end
    if ~isfield(rtmConfigRefine.options,'refinePCFFT')
        rtmConfigRefine.options.refinePCFFT =0;
    end
    if rtmConfigRefine.options.TestPseudoSymmetry==1
        try
            numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
            rotPseudo = rtmConfigRefine.options.rotPseudo;
        catch
            rtmConfigRefine.options.TestPseudoSymmetry=0
            warning('Option Pseudosymmetry specified, but no matrix given, continuing without Pseudosymmetry')
        end
    end

    % fftS
    [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( rtmConfigRefine.Settings_Cor.size(1), rtmConfigRefine.RTM.LPTsize );
    Detector_tilt = rtmConfigRefine.RTM.Rx(InputUser.MicroscopeData.TotalTilt);

    if ~isfield( InputUser, "inputEulerMode")
        InputUser.inputEulerMode = 0;
    end

    if strcmp(InputUser.inputEulerMode, "commercial")
        G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercialDetector",tileNum));
    % else
        % G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/gmatrixCommercialDetector",tileNum));
    end

  %%

    PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
                    tileNum));
    PCBest = zeros(size(PC));
    % phaseId = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/Phase", ...
                    % tileNum)); %
    % % numPatterns = size(Patterns,3);
    numPatterns = size(PC,2);
    
    if ~(InputUser.PCA_Mode == "Full_Map")
        Patterns =h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
        numPatterns = size(Patterns,3);
    else
        Patterns = zeros(1,1,numPatterns); %Variable needs to exist. Otherwise parfor complains.
    end



        %% loop refinement
    %%
    for phaseNumber = 1:numel(InputUser.Phase_Input)%,2)
        %%
         if ~strcmp(InputUser.inputEulerMode, "commercial")
            G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/G_SO3",InputUser.Phase_Input{phaseNumber},tileNum));
        end
        gNew = zeros(3,3,numPatterns);
        Library_PHList = zeros(1,numPatterns);
        PC_best_List =  zeros(3, numPatterns);
        G_Out = zeros(size(gNew));
         if rtmConfigRefine.options.TestPseudoSymmetry==1
            gPseudo = zeros(3,3,numPatterns,numPseudo);
            phPseudo =  zeros(1, numPatterns,numPseudo) ;
        end
      %%
        %TODO: group InputUser.Phases, in such a way, that we can assign
        %different candidates to each phase number.
    % writeDataToH5(InputUser.h5name,, gmatrix_in);
    rtmConfigRefine.RTM.isHex = MultiPhaseStruct(phaseNumber).screen_int.isHex;
    for patternNumber =1:numPatterns   
        % tic
        disp(sprintf('running pattern Number %i',patternNumber))
             if strcmp(InputUser.PCA_Mode, "Full_Map")
                Pat_exp = bReadEBSP(InputUser.MapInfo.EBSPData,patternNumber);
            else
                Pat_exp = Patterns(:,:,patternNumber);
            end

            % phase = phaseId()
            % Pat_exp = Patterns(:)
            %%
            % % for now check if phaseId
            [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,rtmConfigRefine.Settings_Cor);
            % [Pat_Ref_BG] = EBSP_BGCor( Pat_exp_cor,rtmConfigRefine.Settings_Cor );
            [Pat_Ref_r] = refine_prep(Pat_exp_cor,SettingsXCF,rtmConfigRefine.RTM);
            [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                            PC(:,patternNumber));
            % figure
            % imagesc(Pat_exp_cor)
            % %refine
            %%
            if rtmConfigRefine.options.TestPseudoSymmetry ==1
                [G_RefinedPseudo, regout_RPseudo] = testPseudoSymmetry(rtmConfigRefine,MultiPhaseStruct(phaseNumber),...
                Pat_Ref_r, PC(:,patternNumber), ...
                G_SO3(:,:,patternNumber), SettingsXCF);
                PHPseudo = regout_RPseudo(4,:);
                [maxPH,maxPos]=max(PHPseudo(:));
                G_Refined = G_RefinedPseudo(:, : ,maxPos);
                
                regout_R = regout_RPseudo(:, maxPos);
                gPseudo(:,:,patternNumber,:) = G_RefinedPseudo;
                phPseudo(:,patternNumber,:) = PHPseudo;
            else
            
            [G_Refined,regout_R] = refine5(Pat_Ref_r, EBSD_geom,...
                            PC(:,patternNumber), G_SO3(:,:,patternNumber),...
                            SettingsXCF, MultiPhaseStruct(phaseNumber).screen_int,...
                            MultiPhaseStruct(phaseNumber).screen_int.isHex,...
                            rtmConfigRefine.RTM);
            if rtmConfigRefine.options.refinePCFFT ==1

                   [G_Refined, PCBest(:,patternNumber)] = fPCRefine(PC(:,patternNumber),Pat_exp_cor, ...
                        G_SO3(:,:,patternNumber), rtmConfigRefine.Settings_Cor, ...
                        rtmConfigRefine.Refine, SettingsXCF, MultiPhaseStruct(phaseNumber).screen_int,...
                        rtmConfigRefine.RTM,rtmConfigRefine.RTM);
                    [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                            PCBest(:,patternNumber));
                    [G_Refined,regout_R] = refine5(Pat_Ref_r, EBSD_geom,...
                            PCBest(:,patternNumber), G_Refined,...
                            SettingsXCF, MultiPhaseStruct(phaseNumber).screen_int,...
                            MultiPhaseStruct(phaseNumber).screen_int.isHex,...
                            rtmConfigRefine.RTM);
            end
            end
            %%
            EBSD_geom = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                            PC(:,1));
            % [ Pat_sim_in ] = EBSP_gen(EBSD_geom, G_SO3(:,:,4),MultiPhaseStruct(phaseNumber).screen_int);
            [ Pat_sim ] = EBSP_gen(EBSD_geom,G_Refined,MultiPhaseStruct(phaseNumber).screen_int);
%             figure
%             imagesc(EBSP_BGCor(Pat_sim,rtmConfigRefine.Settings_Cor))
%             figure
%             imagesc(EBSP_BGCor(Patterns(:,:,1),rtmConfigRefine.Settings_Cor) )
            %%
            % &&
            % tic
            gTemp = G_SO3(:,:,patternNumber)*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
            gNew(:,:,patternNumber) = gTemp';
            % toc
            G_Out(:,:,patternNumber)=G_Refined;
            Library_PHList(:,patternNumber) = regout_R(4);
            % PC_best_List(:,patternNumber) = PC_best;
            % toc
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
        disp(sprintf("-----------Finished Tile %i      ----------------",tileNum))
        disp(sprintf("--------- mean peak hight : %.4f ----------------",mean(Library_PHList)))
        %%
        exportString = InputUser.Phase_Input{phaseNumber};
        if isfield(InputUser,"paramString")
            exportString =  sprintf("%s/%s", exportString,InputUser.paramString);
        end
        if rtmConfigRefine.options.refinePCFFT ==1
            writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/PCBest", ...
                    exportString,tileNum), PCBest);
        end
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/Refined_PH", ...
                    exportString,tileNum), Library_PHList);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/eangs_out", ...
                        exportString,tileNum), eangs_out_export);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/G_Out", ...
                            exportString,tileNum), G_Out);
        
        if rtmConfigRefine.options.TestPseudoSymmetry ==1
                  writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/G_OutPseudo", ...
                            exportString,tileNum), gPseudo);
                  writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/PH_Pseudo", ...
                            exportString,tileNum), phPseudo);
        end
        end
end