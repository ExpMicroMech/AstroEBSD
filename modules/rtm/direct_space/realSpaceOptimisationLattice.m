function realSpaceOptimisationLattice(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions)
% numPatterns = numel(InputUser.PatternList)
    if nargin <=4
         optimisationOptions.Display = 'iter';
    end


    try
        PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
        tileNum));
    catch
        PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCStatic", ...
        tileNum));
    end

    if ~strcmp(InputUser.inputEulerMode, "commercial")
        InputDataset = "RTMRefined";
    elseif (InputUser.fitPC==0)
        InputDataset = "RTMRefined";
    else
        InputDataset = "RTMPCA";
    end

    OutputDataset = 'RTMRefinedLattice';

    Patterns = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));


numPatterns = size(Patterns,3);
for phaseNumber = 1:numel(InputUser.Phase_Input)
    %%   
    if strcmp(InputDataset,'RTMRefined')
        PC = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/PC", ...
        InputUser.Phase_Input{phaseNumber},tileNum));
    end
    Detector_tilt = rtmConfigRefine.RTM.Rx(InputUser.MicroscopeData.TotalTilt);
        gNew = zeros(3,3,numPatterns);
        PCBestList = zeros(3,numPatterns);
        NCCBestList = zeros(1,numPatterns);
        
        G_RefinedList = h5read(InputUser.h5name,sprintf('/RTM/%s/%s/%i/G_Out',InputUser.Phase_Input{phaseNumber},InputDataset,tileNum));
        PC_bestList = PC;
        parfor patternNumber = 1: numPatterns
        %%
        % G_Refined
        % figure
        % regout_R
        % [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,PC_best,G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);
        G_Refined = G_RefinedList(:,:,patternNumber);
        PC_best = PC_bestList(:,patternNumber);
        %%
        [ EBSD_geom_plot ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,PC_best);
        [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_Refined,MultiPhaseStruct(phaseNumber).screen_int);

        eangs = conv_G_to_EA(G_Refined);
        % Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
        Pat_exp = Patterns(:,:,patternNumber);
        [ Pat_Ref ] = EBSP_BGCor( Pat_exp,rtmConfigRefine.Settings_Cor);
        
        %% pc refinement
        if ~any(isnan(eangs))
        [pc_best, eangs_out, NCCBest,LPRRatio]=optimize_pc_ori(PC_best,...
                    eangs,Pat_Ref,...
                    rtmConfigRefine.Settings_Cor,...
                    MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,'Mode','fmincon','LPRatio');
        %% Ori refinement
        [pc_best, eangs_out, NCCBest,~]=optimize_pc_ori(pc_best, eangs_out, Pat_Ref,...
                        rtmConfigRefine.Settings_Cor, MultiPhaseStruct(phaseNumber).screen_int,...
                        optimisationOptions,'Mode','fmincon');
        else
            pc_best = [0 0 0];
            eangs_out = [NaN NaN NaN];
            NCCBest = 0;
            LPRRatio = [0 0 0]
        end
%%
        gmatrix_out = conv_EA_to_G(eangs_out);
        gTemp = gmatrix_out*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
        gNew(:,:,patternNumber) = gTemp';
        NCCBestList(:,patternNumber) = NCCBest;
        PCBestList(:,patternNumber) = pc_best
        G_Out(:,:,patternNumber)=gmatrix_out;
        LatticeRatios(:,patternNumber) = LPRRatio;
        end
        %%
        % eangs to mtex
        eangs_out = rotation.byMatrix(gNew);
        % the quaternions where causing problems witht he parfor...
        eangs_out_export = zeros(3,numel(eangs_out));
        for i = 1:numPatterns
           eangs_out_export(1,i) =     eangs_out(i).phi1;
           eangs_out_export(2,i) =     eangs_out(i).Phi;
           eangs_out_export(3,i) =     eangs_out(i).phi2;
        end
        %%

        disp("------------------------------------------------------------")
        disp("------------------------------------------------------------")
        disp("############################################################")
        fprintf("Finished refinement of Latticeparameters for %s Tile %i\n",...
                    InputUser.Phase_Input{phaseNumber},tileNum)
        fprintf('mean NCC: %.3f\n',mean(NCCBestList(:)));
        fprintf('mean variation Lattice (x,y,z): %.4f, %.4f, %.4f\n',...
                mean(LatticeRatios(1,:)),mean(LatticeRatios(2,:)),...
                mean(LatticeRatios(3,:)));
        fprintf('std deviation Lattice Params (x,y,z): %.4f, %.4f, %.4f\n',...
                std(LatticeRatios(1,:)),std(LatticeRatios(2,:)),...
                std(LatticeRatios(3,:)));
        

        disp("############################################################")
        disp("------------------------------------------------------------")
        disp("------------------------------------------------------------")
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/RefinedNCC", ...
                        InputUser.Phase_Input{phaseNumber},OutputDataset,tileNum), NCCBestList);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/eangs_out", ...
                        InputUser.Phase_Input{phaseNumber},OutputDataset,tileNum), eangs_out_export);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/G_Out", ...
                         InputUser.Phase_Input{phaseNumber},OutputDataset,tileNum), G_Out);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/PC", ...
                          InputUser.Phase_Input{phaseNumber},OutputDataset,tileNum), PCBestList);
        writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/LatticeRatio", ...
                          InputUser.Phase_Input{phaseNumber},OutputDataset,tileNum), LatticeRatios);


end
