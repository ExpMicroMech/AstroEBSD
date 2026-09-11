function getOrisPCAVMOut(VMOutput,tile,rtmConfigRefine,InputUser)
    numPatterns =VMOutput.NumPCA;
    MapInfo = tile.MapInfo;
    tileNum =tile.tileNum; %hardcoded externally
    [YSearchGrid ,XSearchGrid ]= meshgrid(tile.rowstart(tileNum):1:tile.rowfin(tileNum),tile.colstart(tileNum):1:tile.colfin(tileNum));
    % Plotting starts here)
    % numPatterns = size(numPatterns)
    tileNum = 1; % handcarfte, thus VMOutput only contains 1 tile anyways
    numPhasesEBSD = unique(MapInfo.EBSPData.ebsd_map.phase);
    gmatrix_conv = zeros(3,3,numPatterns);
    gmatrix_in = zeros(3,3,numPatterns);
    PC = zeros(3,numPatterns);
    Phase = zeros(1,numPatterns);
    eangs_in_export = zeros(3,numPatterns);
    try
        numPatterns=tile.numPatterns;
    end
    for pcaNum =1:numPatterns% numPatterns
        %%
        selection = VMOutput.filtered_map'==pcaNum;
        [pnum,~]=find(MapInfo.MapData.XBeam==XSearchGrid(selection)' & MapInfo.MapData.YBeam==YSearchGrid(selection)');
        [pnumOrg,~]=find(MapInfo.MapData.XBeam==XSearchGrid(:)' & MapInfo.MapData.YBeam==YSearchGrid(:)');
        ebsdMapSelection = MapInfo.EBSPData.ebsd_map(pnum);
        %%
        [u, ~, uidx] = unique(ebsdMapSelection.phase);
        counts = accumarray(uidx, 1);
        phaseId = find(max(uidx(counts)));
        try
        Phase(:,pcaNum) = phaseId;
        catch
            warning("phase Id could not be found keeping zero")
        end
        %%
        ebsdMapSelection = ebsdMapSelection(ebsdMapSelection.phase==phaseId);
        pnumNew = ebsdMapSelection.id;

        meanOri = mean(ebsdMapSelection.orientations);
        meanX = mean(MapInfo.MapData.XBeam(pnumNew));
        meanY = mean(MapInfo.MapData.YBeam(pnumNew));
        rtmConfigRefine.Detector_tilt = rtmConfigRefine.RTM.Rx(tile.MicroscopeData.TotalTilt);
        gmatrix_in(:,:,pcaNum) = meanOri.matrix;
        eangs_in_export(1,pcaNum) = meanOri.phi1;
        eangs_in_export(2,pcaNum) = meanOri.Phi;
        eangs_in_export(3,pcaNum) = meanOri.phi2;
        % Patterns_export(:,:,pcaNum) = VMOutput.PCA_VM_Pat_n(:,:,pcaNum);
        gmatrix = gmatrix_in(:,:,pcaNum);
        gmatrix_conv(:,:,pcaNum) = gmatrix'*rtmConfigRefine.Detector_tilt *InputUser.EdaxConversionMatrix ;
        try
           PC(:,pcaNum)  = getPCFromPlane(meanX,meanY,rtmConfigRefine.RTM.PCX_mod,...
            rtmConfigRefine.RTM.PCY_mod,rtmConfigRefine.RTM.PCZ_mod);
        catch
            PC(:,pcaNum)=InputUser.PC
        end
    end
    %%
    Patterns_export = VMOutput.PCA_VM_Pat_n(:,:,1:numPatterns);
    % for i = 1:size(eangs_out,2)
    %    eangs_in_export(i,1) =     eangs_out(i).phi1;
    %    eangs_in_export(i,2) =     eangs_out(i).Phi;
    %    eangs_in_export(i,3) =     eangs_out(i).phi2;
    % end
% for i = 1:
    % for i = 1:numPhasesEBSD
    % writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercial", ...
    %             InputUser.Phase_Input{phaseNumber},tileNum), gmatrix_in);
    %%
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercial", ...
                tileNum), gmatrix_in);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns", ...
                tileNum), single(Patterns_export),"single");
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/EangsCommercial", ...
               tileNum), single(eangs_in_export));
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercialDetector", ...
                tileNum), gmatrix_conv);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
                    tileNum), PC);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/Phase", ...
                    tileNum), Phase,"uint16");
%%
end
