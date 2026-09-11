function realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions)
% numPatterns = numel(InputUser.PatternList)
%  switch nargin
%     case 4
%         tileNum =1
%      otherwise
%         optimisationOptions.Display = 'iter';
% end
if nargin<=4
    optimisationOptions.Display = 'iter';
end
% if nargin<=4

[varargs,varargNum,optimisationOptions] =checkRealSpaceOptions(optimisationOptions);

try
    InputDataset = InputUser.realSpaceInputDataset;
catch
    InputDataset = 'RTMPCA';
end
OutputDataset = 'RTMRefined';
%%
% while 
try
    PC = h5saveRead(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
        tileNum));
catch
    PC = h5saveRead(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCStatic", ...
        tileNum));
end
if strcmp(InputUser.inputEulerMode, "commercial")
    G_SO3 = h5saveRead(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercialDetector",tileNum));
    % else
    % G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/gmatrixCommercialDetector",tileNum));
end 
if ~isfield(rtmConfigRefine.options,'TestPseudoSymmetry')
    rtmConfigRefine.options.TestPseudoSymmetry =0;
end
if rtmConfigRefine.options.TestPseudoSymmetry==1
    try
        numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
        rotPseudo = rtmConfigRefine.options.rotPseudo;
    catch
        rtmConfigRefine.options.TestPseudoSymmetry=0;
        warning('Option Pseudosymmetry specified, but no matrix given, continuing without Pseudosymmetry')
    end
end
%  
%    PC = v
% end
h5name = InputUser.h5name;
% Pattern = h5read(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
orgSize = rtmConfigRefine.options.initPatSize;
% clear Pattern;
%%
numPatterns = size(PC,2); %% need to pass as
for phaseNumber = 1:numel(InputUser.Phase_Input)
    %%
    exportString = InputUser.Phase_Input{phaseNumber};
    if isfield(InputUser,"paramString")
        exportString = sprintf("%s/%s", exportString,InputUser.paramString);
        
    end
    
    if rtmConfigRefine.options.TestPseudoSymmetry==1
        gPseudo = zeros(3,3,numPatterns,numPseudo);
        phPseudo =  zeros(1, numPatterns,numPseudo) ;
    end
    
    Detector_tilt = rtmConfigRefine.RTM.Rx(InputUser.MicroscopeData.TotalTilt);
    gNew = zeros(3,3,numPatterns);
    PCBestList = zeros(3,numPatterns);
    NCCBestList = zeros(1,numPatterns);
    try
        if strcmp(InputDataset,"RTMCoarse")
             G_RefinedList = h5saveRead(InputUser.h5name,sprintf('/RTM/%s/%s/%i/G_SO3',exportString,InputDataset,tileNum));
        else
            G_RefinedList = h5saveRead(InputUser.h5name,sprintf('/RTM/%s/%s/%i/G_Out',exportString,InputDataset,tileNum));
        end
    catch
        G_RefinedList = G_SO3;
    end
    PC_bestList = PC;
    %%
    if optimisationOptions.parFor==0
        for patternNumber = 1: numPatterns
            try
                if InputUser.PCA_Mode == "Full_Map"
                    Pat_exp = bReadEBSP(InputUser.MapInfo.EBSPData,patternNumber);
                else
                    try
                        Pat_exp = h5saveRead(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum),[1, 1,patternNumber],[orgSize(1), orgSize(2),1]);
                    catch
                        if patternNumber ==1
                            Pat_exp = h5saveRead(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum),[1, 1],[orgSize(1), orgSize(2)]);
                        % end
                        else
                        warning("Patterns could not be loaded, currently we are not catching this error. Please look into the scource code to fix.")
                        end
                    end
                end
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
                % imwrite(Pat_sim,fullfile(InputUser.h5folder, sprintf("Pat_%i.png",patternNumber)));
    
                % disp(PC_best);
                % disp(G_Refined);
                eangs = conv_G_to_EA(G_Refined);
                % disp(eangs);
                % Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
                % Pat_exp = Patterns(:,:,patternNumber);
                % disp(rtmConfigRefine.Settings_Cor);
                [ Pat_Ref ] = EBSP_BGCor( Pat_exp,rtmConfigRefine.Settings_Cor);
                % imwrite(Pat_Ref,fullfile(InputUser.h5folder, sprintf("Patexp_%i.png",patternNumber)));
                % error('end')
                tic
                %%
                if ~any(isnan(eangs))
                  
                    %% pc refinement
             
                    % gmatrix_out = conv_EA_to_G(eangs_out');
                    % test pseudo Symmetry
                    if rtmConfigRefine.options.TestPseudoSymmetry ==1
                        % eangs = conv_G_to_EA(gmatrix_out);
                        [G_RefinedPseudo, NCCPseudo,pC_Out] = testPseudoSymmetryReal(rtmConfigRefine,MultiPhaseStruct(phaseNumber),...
                        Pat_Ref, PC_best, ...
                        G_Refined, optimisationOptions,varargs);
        
                        % NCCBest = regout_RPseudo(4,:);
                        [maxPH,maxPos]=max(NCCPseudo(:));
                        G_Refined = G_RefinedPseudo(:, : ,maxPos);
                        PC_best = pC_Out(maxPos,:);
                        NCCBest = maxPH;
                        % gPseudo(:,:,patternNumber,:) = gmatrix_out;
                        % phPseudo(:,patternNumber,:) = PHPseudo;
    
                        % gmatrix_out = 
                        
                        %  [PC_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                        % eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                        % MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                    end
                    
                    eangs = conv_G_to_EA(G_Refined);
                    [pc_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                        eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                        MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                    % end
                    try
                        if ~ InputUser.PCA_Mode == "Full_Map"
                        writeStructToH5(optiparams,InputUser.h5name,sprintf("/RTM/%s/%s/%i/params/%i", ...
                            exportString,OutputDataset,tileNum,patternNumber))
                        end
                    end
                else
                    % toc
                    pc_best = [0 0 0];
                    eangs_out = [NaN NaN NaN];
                    NCCBest = 0;
                    gmatrix_out = conv_EA_to_G(eangs_out');
                end
                
                
                toc
                %%
                gmatrix_out = conv_EA_to_G(eangs_out');
                gTemp = gmatrix_out*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
                gNew(:,:,patternNumber) = gTemp';
                NCCBestList(:,patternNumber) = NCCBest;
                PCBestList(:,patternNumber) = pc_best;
                G_Out(:,:,patternNumber)=gmatrix_out;
                %%
                catch e %e is an MException struct
                    %% basically is just here to catch "wrong pc handling and thus wrong number of patterns, might be a better way to handle this of course.
                    fprintf(1,'The identifier was:\n%s',e.identifier);
                    fprintf(1,'There was an error! The message was:\n%s',e.message);
                    % more error handling...
                    break
                end
        end
        % end

    else %%
        parfor patternNumber = 1: numPatterns
            if mod(patternNumber,100)==0

                tic
            end
           if InputUser.PCA_Mode == "Full_Map"
                Pat_exp = bReadEBSP(InputUser.MapInfo.EBSPData,patternNumber);
            else
                Pat_exp = h5saveRead(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum),[1, 1,patternNumber],[orgSize(1), orgSize(2),1]);
            end
            %%% Pat_exp = squeeze(Pat_exp);
            % size(Pat_exp)
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
            % imwrite(Pat_sim,fullfile(InputUser.h5folder, sprintf("Pat_%i.png",patternNumber)));

            % if patternNumber ==1
            %     disp('PCBest:')
            %     disp(PC_best);
            %     disp('G_Refined:')
            %     disp(G_Refined);
            %     disp(eangs);
            % end

            % Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber);
            % Pat_exp = Patterns(:,:,patternNumber);
            % imwrite(Pat_exp,fullfile(InputUser.h5folder, sprintf("Patexp_%i.png",patternNumber)))
            [ Pat_Ref ] = EBSP_BGCor( Pat_exp,rtmConfigRefine.Settings_Cor);
            % disp(rtmConfigRefine.Settings_Cor);
            % imwrite(Pat_Ref,fullfile(InputUser.h5folder, sprintf("PatexpCor_%i.png",patternNumber)));
            % error('end')
             if ~any(isnan(eangs))
              
                %% pc refinement
         
                % gmatrix_out = conv_EA_to_G(eangs_out');
                % test pseudo Symmetry
                if rtmConfigRefine.options.TestPseudoSymmetry ==1
                    % eangs = conv_G_to_EA(gmatrix_out);
                    [G_RefinedPseudo, NCCPseudo,pC_Out] = testPseudoSymmetryReal(rtmConfigRefine,MultiPhaseStruct(phaseNumber),...
                    Pat_Ref, PC_best, ...
                    G_Refined, optimisationOptions,varargs);
    
                    % NCCBest = regout_RPseudo(4,:);
                    [maxPH,maxPos]=max(NCCPseudo(:));
                    G_Refined = G_RefinedPseudo(:, : ,maxPos);
                    PC_best = pC_Out(maxPos,:);
                    NCCBest = maxPH;
                    % gPseudo(:,:,patternNumber,:) = gmatrix_out;
                    % phPseudo(:,patternNumber,:) = PHPseudo;

                    % gmatrix_out = 
                    
                    %  [PC_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                    % eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                    % MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                end
                
                eangs = conv_G_to_EA(G_Refined);
                [pc_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                    eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                    MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                % end
                %% TODO: make this optional
                 try
                        if ~ InputUser.PCA_Mode == "Full_Map"
                        writeStructToH5(optiparams,InputUser.h5name,sprintf("/RTM/%s/%s/%i/params/%i", ...
                            exportString,OutputDataset,tileNum,patternNumber))
                        end
                    end

            else
                % toc
                pc_best = [0 0 0];
                eangs_out = [NaN NaN NaN];
                NCCBest = 0;
                gmatrix_out = conv_EA_to_G(eangs_out');
            end
            %%
            gmatrix_out = conv_EA_to_G(eangs_out);
            gTemp = gmatrix_out*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
            gNew(:,:,patternNumber) = gTemp';
            NCCBestList(:,patternNumber) = NCCBest;
            PCBestList(:,patternNumber) = pc_best;
            G_Out(:,:,patternNumber)=gmatrix_out;
            if mod(patternNumber,100)==0

                toc
                disp(sprintf("Running pattern number: %i",patternNumber))
            end
            % toc
        end
    end
    %%
    % profile off
    % error('end')
    % eangs to mtex
    eangs_out = rotation.byMatrix(gNew);
    % the quaternions where causing problems witht he parfor...
    eangs_out_export = zeros(3,numel(eangs_out));
    for i = 1:size(eangs_out,2)
        eangs_out_export(1,i) =     eangs_out(i).phi1';
        eangs_out_export(2,i) =     eangs_out(i).Phi';
        eangs_out_export(3,i) =     eangs_out(i).phi2';
    end
    %%
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/RefinedNCC", ...
        exportString,OutputDataset,tileNum), NCCBestList);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/eangs_out", ...
        exportString,OutputDataset,tileNum), eangs_out_export);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/G_Out", ...
        exportString,OutputDataset,tileNum), G_Out);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/%s/%s/%i/PC", ...
        exportString,OutputDataset,tileNum), PCBestList);
    
    
    
    
end
end
function EBSP2=fix_mean(EBSP)
%zero mean & fix stdev
% EBSP2=(EBSP-mean(EBSP(:),"omitnan"))/std(EBSP(:),"omitnan");
EBSP2=(EBSP-min(EBSP(:)))/max(EBSP(:));
%make positive
% EBSP2=EBSP2-min(EBSP2(:))+1;
end