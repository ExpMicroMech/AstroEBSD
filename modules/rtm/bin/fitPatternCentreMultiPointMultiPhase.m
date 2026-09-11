function [PCBest,eangs_out,phbest,eangs] = fitPatternCentreMultiPoint(MapInfo,MicroscopeData,InputUser,Settings_Cor,screen_int,patternsRefine,pnum,options)
    % Rx  = InputUser.RTM.Rx
    % Ry = InputUser.RTM.Ry
    RTM = InputUser.RTM;
    PC = [InputUser.calibration.pcx,InputUser.calibration.pcy,InputUser.calibration.pcz];
    ebsd_map = MapInfo.EBSPData.ebsd_map;
   % cmap_div = returnCmapDiv;
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
 
 
    % writeStructToHDF5(options, h5name)
    %
    eangs = ebsd_map(pnum).orientations;
    eangsList = [eangs.phi1, eangs.Phi,eangs.phi2];
    pcbest = zeros(size(eangs,1),3);
    phbest = zeros(size(eangs,1),1);
    eangs_out = zeros(size(eangs,1),3);
    gNew = zeros(3,3,size(eangs,1));
    % profile on
    gmatrix_in=eangs.matrix;
    %%
    % tic
    % disp(size(pnum))
    parfor i = 1:size(pnum(:),1)
        tic
        
        Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
        gmatrix = gmatrix_in(:,:,i);
        gmatrix = gmatrix'*Detector_tilt*InputUser.EdaxConversionMatrix ; %% only for EDAX: *RTM.Ry(180*degree)
        eangsIn=conv_G_to_EA(gmatrix);
        %
        % options = struct()
        [patternIn,PatternInfo] = EBSP_BGCor(patternsRefine(:,:,i),Settings_Cor);
        [PCBest(i,:),eangs_o(i,:),phbest(i),~] = optimize_pc_ori(PC,eangsIn,patternIn,Settings_Cor,screen_int,options,'Mode','fmincon','refinePC');%,'LPRatio');
        toc
        % outputs
%%
        gmatrix_out = conv_EA_to_G(eangs_o(i,:));
        gTemp = gmatrix_out*inv(InputUser.EdaxConversionMatrix)*inv(Detector_tilt);
        gNew(:,:,i) = gTemp';
        % eangs_Converted = conv_G_to_EA(gNew);
 
%% 
        % writeStructToHDF5(options, h5name)
        %%
        % eangs_export = rotation.byMatrix(gNew(:,:,i))
        EBSD_geometry = EBSP_Gnom( PatternInfo,PCBest(i,:));
        patternSim = EBSP_gen(EBSD_geometry,conv_EA_to_G(eangs_o(i,:)),screen_int,0);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/ExperimentalPattern",InputUser.hdf5DataSet,i), flipud(patternIn))
        
   
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/SimulatedPattern",...
                    InputUser.hdf5DataSet,i), flipud(patternSim));
        
        % writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/Eangs_out (radians)",...
        %             InputUser.hdf5DataSet,i), [eangs_export.phi1 eangs_export.Phi eangs_export.phi2])
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/Eangs_in (radians)",...
                    InputUser.hdf5DataSet,i), [eangsList(i,:)]);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/PC_in pcx,pcxy,pcz", ...
                    InputUser.hdf5DataSet,i), PC);
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/NCC", ...
                    InputUser.hdf5DataSet,i),phbest(i));
        
        writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/PC_out (pcx,pcy,pcz)",...
                    InputUser.hdf5DataSet,i),PCBest(i,:))
    end
    eangs_out = rotation.byMatrix(gNew);
    % the quaternions where causing problems witht he parfor...
    for i = 1:size(eangs_out,2)
                writeDataToH5(h5name,sprintf("%s/Patterns/Pattern_%04i/Eangs_out (radians)",...
                    InputUser.hdf5DataSet,i), [eangs_out(i).phi1 eangs_out(i).Phi eangs_out(i).phi2])
    end

    % profile off
    % profile viewer
end