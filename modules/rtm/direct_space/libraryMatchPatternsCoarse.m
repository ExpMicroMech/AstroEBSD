function libraryMatchPatternsCoarse(InputUser,rtmConfig,MultiPhaseStruct,tileNum)
    

    % switch nargin
        % % case 4
            % pass
        % % otherwise
            tileNum=1
    % end
    if nargin <4
        tileNum = 1
    end
    [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( ...
                    rtmConfig.RTM.screensize, rtmConfig.RTM.LPTsize );
    for fn = fieldnames(rtmConfig.Settings_Cor)'
       SettingsXCF.(fn{1}) =rtmConfig.Settings_Cor.(fn{1});
    end
    [~,XCF_data_fill] = refine_prep(zeros(rtmConfig.Settings_Cor.size,rtmConfig.Settings_Cor.size), ...
                    SettingsXCF, rtmConfig.RTM);
       % Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
    if strcmp(InputUser.PCA_Mode, "Full_Map")
        warning('DEBUG / TESTING mode. RUNNIGN WITh only one pattern')
        Patterns = zeros(1,1,1)
    else

        Patterns = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
    end
        try
        PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
            tileNum));
    catch
        PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCStatic", ...
            tileNum));
    end
    %     PC = v
    % end
    numPatterns = size(Patterns,3);
    for phaseNumber = 1:numel(InputUser.Phase_Input)
    
        %% Build the library
        tic
        [ library_G ] = SO3_rotmat_gen(MultiPhaseStruct(phaseNumber).cs_phase,...
                    rtmConfig.RTM.Sampling_Freq);
        [ template_library ] = Library_Gen(MultiPhaseStruct(phaseNumber).EBSD_geometry ,...
                    MultiPhaseStruct(phaseNumber).screen_int, ...
                    MultiPhaseStruct(phaseNumber).screen_int.isHex, library_G, 2, SettingsXCF);
        radial_adjustment=fPHFac(  MultiPhaseStruct(phaseNumber).EBSD_geometry,... 
                     rtmConfig.Settings_Cor, MultiPhaseStruct(phaseNumber).screen_int,...
                     rtmConfig.RTM,SettingsXCF);
        toc
        %%
        % numPatterns = numel(InputUser.PatternList)
        G_SO3_List = zeros([3,3,numPatterns]);
        Library_PHList = zeros([1, numPatterns]);
        for patternNumber = 1: numPatterns
                %%
            if strcmp(InputUser.PCA_Mode, "Full_Map")
                 Pat_exp= bReadEBSP(InputUser.MapInfo.EBSPData,1);
            else
                Pat_exp = Patterns(:,:,patternNumber);
            end

            %%
            [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,rtmConfig.Settings_Cor);
            [Pat_Ref_BG] = EBSP_BGCor( Pat_exp_cor,rtmConfig.Settings_Cor );
            [Pat_Ref_r] = refine_prep(Pat_Ref_BG,SettingsXCF,rtmConfig.RTM);
            %%
            tic
            [G_SO3_List(:,:,patternNumber), Library_PHList(patternNumber)] = ...,
                fLibraryTest(template_library, library_G, Pat_Ref_r.FFT,...
                SettingsXCF, XCF_data_fill, 1);
            toc
        end
        %%
            writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/G_SO3",...
                            InputUser.Phase_Input{phaseNumber},tileNum), [G_SO3_List]);
                
            writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/Library_PH", ...
                        InputUser.Phase_Input{phaseNumber},tileNum), Library_PHList);
    
             % writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/Library_PH", ...
             %            InputUser.Phase_Input{phaseNumber},tileNum), InputUser.PatternList);
    
            % writeDataToH5(h5name,sprintf("%s/RTM/NCC", ...
            %             InputUser.Phase_Input{phaseNumber}),phbest);
    end
    toc
end