function [PCBest,eangs_out,phbest,eangs] = fitPatternCentreMultiPointNonPar(MapInfo,MicroscopeData,InputUser,Settings_Cor,screen_int,patternsRefine,pnum,options)
    % Rx  = InputUser.RTM.Rx
    % Ry = InputUser.RTM.Ry
    RTM = InputUser.RTM;
    PC = [InputUser.calibration.pcx,InputUser.calibration.pcy,InputUser.calibration.pcz];
    ebsd_map = MapInfo.EBSPData.ebsd_map;
    if ~isfield(options,'gradients')
        options.gradients = 0
    end
    if ~isfield(options,'fminsearch')
        options.fminsearch = 0
    end
    % cmap_div = returnCmapDiv;
    %%
    % InputUser
    %  ebsd_size = [MicroscopeData.NCOLS,MicroscopeData.NROWS]
    % [col,row] = ind2sub(ebsd_size,pnum)
    % pnumTransformed = sub2ind((ebsd_size), row,col)  
    %
    h5name = InputUser.h5name;
    % try
    
    
    % try
    %     writeStructToH5(options, h5name,sprintf("%s/options",InputUser.hdf5DataSet))
    % catch
    %  warning("Could not write options to h5")
    % end
 
 
    % writeStructToHDF5(options, h5name)
    %
    eangs = ebsd_map(pnum).orientations;
    eangsList = [eangs.phi1, eangs.Phi,eangs.phi2]';
    pcbest = zeros(3,size(eangs,1));
    phbest = zeros(1,size(eangs,1));
    eangs_out = zeros(3,size(eangs,1));
    gNew = zeros(3,3,size(eangs,1));
    % profile on
    gmatrix_in=eangs.matrix;
    %%
    patternSim = zeros(Settings_Cor.size(1),Settings_Cor.size(2),size(eangs,1));
    patternsCorrected =  zeros(Settings_Cor.size(1),Settings_Cor.size(2),size(eangs,1));
    %%
    % tic
    % disp(size(pnum))
    for i = 1:size(pnum(:),1)
        tic
         
        Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
        gmatrix = gmatrix_in(:,:,i);
        gmatrix = gmatrix'*Detector_tilt*InputUser.EdaxConversionMatrix ; %% only for EDAX: *RTM.Ry(180*degree)
        eangsIn=conv_G_to_EA(gmatrix);
        %
        % options = struct()
        [patternIn,PatternInfo] = EBSP_BGCor(patternsRefine(:,:,i),Settings_Cor);
        patternsCorrected(:,:,i) = patternIn;
        %%
        %%
        if ~options.fminsearch ==1
            if options.gradients == 1
                [PCBest(:,i),eangs_o(:,i),phbest(:,i),~] = optimize_pc_ori_rtm(PC,eangsIn,patternIn,Settings_Cor,screen_int,options,'Mode','fmincon','refinePC','gradients');%,'LPRatio');
            else
                            [PCBest(:,i),eangs_o(:,i),phbest(:,i),~] = optimize_pc_ori_rtm(PC,eangsIn,patternIn,Settings_Cor,screen_int,options,'Mode','fmincon','refinePC');%,'LPRatio');
            end
            %% we need to find a bettwer way of implementing
        else
         if options.gradients == 1
            [PCBest(:,i),eangs_o(:,i),phbest(:,i),~] = optimize_pc_ori_rtm(PC,eangsIn,patternIn,Settings_Cor,screen_int,options,'Mode','fminsearch','refinePC','gradients');%,'LPRatio');
        else
                        [PCBest(:,i),eangs_o(:,i),phbest(:,i),~] = optimize_pc_ori_rtm(PC,eangsIn,patternIn,Settings_Cor,screen_int,options,'Mode','fminsearch','refinePC');%,'LPRatio');
         end
        end
            toc
        % outputs
%%
        gmatrix_out = conv_EA_to_G(eangs_o(:,i));
        gTemp = gmatrix_out*inv(InputUser.EdaxConversionMatrix)*inv(Detector_tilt);
        gNew(:,:,i) = gTemp';
        % eangs_Converted = conv_G_to_EA(gNew);
 
%% 
        % writeStructToHDF5(options, h5name)
        %%
        % eangs_export = rotation.byMatrix(gNew(:,:,i))
        EBSD_geometry = EBSP_Gnom( PatternInfo,PCBest(:,i));
        patternSim(:,:,i) = EBSP_gen(EBSD_geometry,conv_EA_to_G(eangs_o(:,i)),screen_int,0);
        
        
    end
    %%
    eangs_out = rotation.byMatrix(gNew);
    % the quaternions where causing problems witht he parfor...
    eangs_out_export = zeros(3,size(eangs_out,1));
    for i = 1:size(eangs_out,2)
       eangs_out_export(i,1) =     eangs_out(i).phi1;
       eangs_out_export(i,2) =     eangs_out(i).Phi;
       eangs_out_export(i,3) =     eangs_out(i).phi2;
    end
eangs_out_export = eangs_out_export'
    %% replicate some data before saving
      writeDataToH5(h5name,sprintf("%s/Patterns/Eangs_out (radians)",...
                    InputUser.hdf5DataSet), eangs_out_export)
    PC_IN = repmat(PC, size(phbest'))';
    %%
    
        % writeDataToH5(h5name,sprintf("%s/Patterns/PatternsEXP",InputUser.hdf5DataSet), patternsCorrected);
        
   
        % writeDataToH5(h5name,sprintf("%s/Patterns/SimulatedPattern",...
        %             InputUser.hdf5DataSet), patternSim);
        
        % writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/Eangs_out (radians)",...
        %             InputUser.hdf5DataSet,i), [eangs_export.phi1 eangs_export.Phi eangs_export.phi2])
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Eangs_in (radians)",...
                    InputUser.hdf5DataSet), [eangsList]);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/PC_in (pcx,pcxy,pcz)", ...
                    InputUser.hdf5DataSet), PC_IN);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/NCC", ...
                    InputUser.hdf5DataSet),phbest);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/PC_out (pcx,pcy,pcz)",...
                    InputUser.hdf5DataSet),PCBest)
    % profile off
    % profile viewer
end