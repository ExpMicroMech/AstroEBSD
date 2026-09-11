function [pcbest,eangs_out,phbest,eangs] = fitPatternCentreMultiPoint(MapInfo,MicroscopeData,InputUser,Settings_Cor,screen_int,patternsRefine,pnum,options)
    % Rx  = InputUser.RTM.Rx
    % Ry = InputUser.RTM.Ry
    RTM = InputUser.RTM;
    PC = [InputUser.calibration.pcx,InputUser.calibration.pcy,InputUser.calibration.pcz];
    if isfield(MapInfo.MapData,"pcx")
        PC = [MapInfo.MapData.pcx(:)'; MapInfo.MapData.pcy(:)'; MapInfo.MapData.pcz(:)'];
        mappedPC = 1;
    else
        pcSel = PC;
    end
    %%
    ebsd_map = MapInfo.EBSPData.ebsd_map;
    % cmap_div = returnCmapDiv;
    if ~isfield(options,'gradients')
        options.gradients = 0
    end
    if ~isfield(options,'fminsearch')
        options.fminsearch = 0
    end
    if ~isfield(options,'refinePC')
    options.refinePC=1
    % this is actually supposed to be one as default
    end
    if ~isfield(options,'genetic')
        options.genetic=0
    end
    if ~isfield(options,'gradients')
    options.gradients = 0;
    end
    varargNum = 1
    % varargs = {""}
    if options.fminsearch ==1
        varargs{varargNum} = 'Mode';
        varargs{varargNum+1} = 'fminsearch';
        % varargNum = varargNum+
    end
    % this way the first option will be overwritten, by the later one.
    if options.genetic ==1
        varargs{varargNum} = 'Mode';
        varargs{varargNum+1} = 'genetic';
    end
    varargNum = varargNum+2


    if options.gradients ==1
        varargs{varargNum} = 'gradients';
        varargNum = varargNum+1;
    end
    if options.refinePC ==1
        varargs{varargNum} = 'refinePC';
        varargNum = varargNum+1;
    end
    %%
    % InputUser
    %  ebsd_size = [MicroscopeData.NCOLS,MicroscopeData.NROWS]
    % [col,row] = ind2sub(ebsd_size,pnum)
    % pnumTransformed = sub2ind((ebsd_size), row,col)  
    %
    h5name = InputUser.h5name;
    % try
    
    
    try
        writeStructToH5(options, h5name,sprintf("%s/options",InputUser.hdf5DataSet))
    catch
     warning("Could not write options to h5")
    end

    if ~isfield(options,'gradients')
        options.gradients == 0;
    end
 
 
    % writeStructToHDF5(options, h5name)
    %
    eangs = ebsd_map(pnum).orientations;
    eangsList = [eangs.phi1, eangs.Phi,eangs.phi2]';
    pcbest = zeros(3,size(eangs,1));
    phbest = zeros(1,size(eangs,1));
    eangs_o = zeros(3,size(eangs,1));
    gNew = zeros(3,3,size(eangs,1));
    % profile on
    gmatrix_in=eangs.matrix;
    %%
    patternSim = zeros(Settings_Cor.size(1),Settings_Cor.size(2),size(eangs,1));
    patternsCorrected =  zeros(Settings_Cor.size(1),Settings_Cor.size(2),size(eangs,1));
    %%
    % tic
    % disp(size(pnum))
    parfor patternNum = 1:numel(pnum(:))
        tic
         
        Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
        gmatrix = gmatrix_in(:,:,patternNum);
        gmatrix = (gmatrix'*Detector_tilt)*InputUser.EdaxConversionMatrix ; %% only for EDAX: *RTM.Ry(180*degree)
        eangsIn=conv_G_to_EA(gmatrix);
        %
        % options = struct()
        [patternIn,PatternInfo] = EBSP_BGCor(patternsRefine(:,:,patternNum),Settings_Cor);
        patternsCorrected(:,:,patternNum) = patternIn;
        if mappedPC == 1
            pcSel = PC(:,pnum(patternNum));
        end
        %%
        %%
        %   if ~any(isnan(eangs))
            [pcbest(:,patternNum),eangs_o(:,patternNum),phbest(:,patternNum),~] = optimize_pc_ori_rtm(pcSel,eangsIn,patternIn,Settings_Cor,screen_int,options,varargs{:});%,'LPRatio');
            %% pc refinement
            % writeStructToH5(optiparams,InputUser.h5name,sprintf("/RTM/%s/%s/%i/params/%i", ...
            %     exportString,OutputDataset,tileNum,patternNumber))
            % else
            % % else
            % % else
            %     PCBest = [0 0 0];
            %     eangs_out = [NaN NaN NaN];
            %     NCCBest = 0;
            % end
            % outputs
    %%
        gmatrix_out = conv_EA_to_G(eangs_o(:,patternNum));
        gTemp = gmatrix_out*inv(InputUser.EdaxConversionMatrix)*inv(Detector_tilt);
        gNew(:,:,patternNum) = gTemp';
        % eangs_Converted = conv_G_to_EA(gNew);
 
%% 
        % writeStructToHDF5(options, h5name)
        %%
        % eangs_export = rotation.byMatrix(gNew(:,:,i))
        EBSD_geometry = EBSP_Gnom( PatternInfo,pcbest(:,patternNum));
        patternSim = EBSP_gen(EBSD_geometry,conv_EA_to_G(eangs_o(:,patternNum)),screen_int,0);
        patternSim = EBSP_BGCor(patternSim,Settings_Cor)
        
    end
    %%
    eangs_out = rotation.byMatrix(gNew);
    % the quaternions where causing problems witht he parfor...
    eangs_out_export = zeros(3,size(eangs_out,1));
    for eangsNum = 1:size(eangs_out,2)
       eangs_out_export(eangsNum,1) =     eangs_out(eangsNum).phi1;
       eangs_out_export(eangsNum,2) =     eangs_out(eangsNum).Phi;
       eangs_out_export(eangsNum,3) =     eangs_out(eangsNum).phi2;
    end
    eangs_out_export = eangs_out_export'
    %% replicate some data before saving
      writeDataToH5(h5name,sprintf("%s/Patterns/Eangs_out (radians)",...
                    InputUser.hdf5DataSet), eangs_out_export)
    % PC_IN = repmat(PC, size(phbest'))';
    PC_IN = PC
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
                    InputUser.hdf5DataSet),pcbest)
    % profile off
    % profile viewer
end