function realSpaceOptimisation(InputUser,rtmConfigRefine,MultiPhaseStruct,tileNum,optimisationOptions)
% numPatterns = numel(InputUser.PatternList)
if nargin<=4
    optimisationOptions.Display = 'iter';
end
if ~isfield(optimisationOptions,'fminsearch')
    optimisationOptions.fminsearch=0;
end

if ~isfield(optimisationOptions,'refinePC')
    optimisationOptions.refinePC=1;
    % this is actually supposed to be one as default
end
if ~isfield(optimisationOptions,'genetic')
    optimisationOptions.genetic=0;
end
if ~isfield(optimisationOptions,'particleSwarm')
    optimisationOptions.particleSwarm=0;
end
if ~isfield(optimisationOptions,'patternsearch')
    optimisationOptions.patternsearch=0;
end
if ~isfield(optimisationOptions,'snobfit')
    optimisationOptions.snobfit=0;
end
if ~isfield(optimisationOptions,'gradients')
    optimisationOptions.gradients = 0;
end
if ~isfield(optimisationOptions,['doPlot' ...
        '' ...
        ''])
    optimisationOptions.doPlot = 0;
end
if ~isfield(optimisationOptions, 'parFor')
    optimisationOptions.parFor = 1; % parallel default;
    optimisationOptions.UseParallel = 0; %do no tuse parallel in inner loops.               
end
% optimisationOptions.refinePC ==1

InputDataset = 'RTMPCA';
OutputDataset = 'RTMRefined';

% set up variable argumetns
varargNum = 1;
% varargs = {""}
if optimisationOptions.fminsearch ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'fminsearch';
    % varargNum = varargNum+
end
% this way the first option will be overwritten, by the later one.
if optimisationOptions.genetic ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'genetic';
end
if optimisationOptions.snobfit ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'snobfit';
end
if optimisationOptions.particleSwarm ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'particleSwarm';
end
if optimisationOptions.patternsearch ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'patternsearch';
end
varargNum = varargNum+2;


if optimisationOptions.gradients ==1
    varargs{varargNum} = 'gradients';
    varargNum = varargNum+1;
end
if optimisationOptions.refinePC ==1
    varargs{varargNum} = 'refinePC';
    varargNum = varargNum+1;
end
if optimisationOptions.doPlot ==1
    varargs{varargNum} ='doPlot'
    varargNum = varargNum+1;
end
if ~exist("varargs")
    varargs = "";
end





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
%     PC = v
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
    
    Detector_tilt = rtmConfigRefine.RTM.Rx(InputUser.MicroscopeData.TotalTilt);
    gNew = zeros(3,3,numPatterns);
    PCBestList = zeros(3,numPatterns);
    NCCBestList = zeros(1,numPatterns);
    try
        G_RefinedList = h5saveRead(InputUser.h5name,sprintf('/RTM/%s/%s/%i/G_Out',exportString,InputDataset,tileNum));
    catch
        G_RefinedList = G_SO3;
    end
    PC_bestList = PC;
    %%
    if optimisationOptions.parFor==0
        for patternNumber = 1: numPatterns
            Pat_exp = h5saveRead(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum),[1, 1,patternNumber],[orgSize(1), orgSize(2),1]);
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
            if ~any(isnan(eangs))
                [pc_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                    eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                    MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                %% pc refinement
                writeStructToH5(optiparams,InputUser.h5name,sprintf("/RTM/%s/%s/%i/params/%i", ...
                    exportString,OutputDataset,tileNum,patternNumber))
            else
                pc_best = [0 0 0];
                eangs_out = [NaN NaN NaN];
                NCCBest = 0;
            end
            %%
            gmatrix_out = conv_EA_to_G(eangs_out);
            gTemp = gmatrix_out*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
            gNew(:,:,patternNumber) = gTemp';
            NCCBestList(:,patternNumber) = NCCBest;
            PCBestList(:,patternNumber) = pc_best;
            G_Out(:,:,patternNumber)=gmatrix_out;
        end

    else %%
        parfor patternNumber = 1: numPatterns
            Pat_exp = h5saveRead(h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum),[1, 1,patternNumber],[orgSize(1), orgSize(2),1]);
            % Pat_exp = squeeze(Pat_exp);
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
                [pc_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC_best,...
                    eangs,Pat_Ref,rtmConfigRefine.Settings_Cor,...
                    MultiPhaseStruct(phaseNumber).screen_int,optimisationOptions,varargs{:});
                %% pc refinement
                writeStructToH5(optiparams,InputUser.h5name,sprintf("/RTM/%s/%s/%i/params/%i", ...
                    exportString,OutputDataset,tileNum,patternNumber))
            else
                pc_best = [0 0 0];
                eangs_out = [NaN NaN NaN];
                NCCBest = 0;
            end
            %%
            gmatrix_out = conv_EA_to_G(eangs_out);
            gTemp = gmatrix_out*inv(rtmConfigRefine.RTM.EdaxConversionMatrix)*inv(Detector_tilt);
            gNew(:,:,patternNumber) = gTemp';
            NCCBestList(:,patternNumber) = NCCBest;
            PCBestList(:,patternNumber) = pc_best;
            G_Out(:,:,patternNumber)=gmatrix_out;
        end
    end
    %%
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
