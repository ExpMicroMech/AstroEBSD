function preparePCARefineTile(VMOutput,tile,rtmConfigRefine,InputUser)
    numPatterns =VMOutput.NumPCA;
    MapInfo = tile.MapInfo;
    tileNum =tile.tileNum; %hardcoded externally
    [YSearchGrid ,XSearchGrid ]= meshgrid(tile.rowstart(tileNum):1:tile.rowfin(tileNum),tile.colstart(tileNum):1:tile.colfin(tileNum));
    % Plotting starts here)
    % numPatterns = size(numPatterns)
    % tileNum = 1; % handcarfte, thus VMOutput only contains 1 tile anyways
    numPhasesEBSD = unique(MapInfo.EBSPData.ebsd_map.phase);
    % gmatrix_conv = zeros(3,3,numPatterns);
    % gmatrix_in = zeros(3,3,numPatterns);
    % PC = zeros(3,numPatterns);
    % % Phase = zeros(1,numPatterns);
    % eangs_in_export = zeros(3,numPatterns);
    try
        numPatterns=tile.numPatterns;
    end
    for pcaNum =1:numPatterns% numPatterns
        selection = VMOutput.filtered_map'==pcaNum;
        [pnum,~]=find(MapInfo.MapData.XBeam==XSearchGrid(selection)' & MapInfo.MapData.YBeam==YSearchGrid(selection)');
        [pnumOrg,~]=find(MapInfo.MapData.XBeam==XSearchGrid(:)' & MapInfo.MapData.YBeam==YSearchGrid(:)');
        ebsdMapSelection = MapInfo.EBSPData.ebsd_map(pnum);
        %%
        [u, ~, uidx] = unique(ebsdMapSelection.phase);
        counts = accumarray(uidx, 1);
        phaseId = find(max(uidx(counts)));
        Phase(:,pcaNum) = phaseId;
        %%
        ebsdMapSelection = ebsdMapSelection(ebsdMapSelection.phase==phaseId);
        pnumNew = ebsdMapSelection.id;

        meanOri = mean(ebsdMapSelection.orientations);
        meanX = mean(MapInfo.MapData.XBeam(pnumNew));
        meanY = mean(MapInfo.MapData.YBeam(pnumNew));
        try
           PC(:,pcaNum)  = getPCFromPlane(meanX,meanY,rtmConfigRefine.RTM.PCX_mod,...
            rtmConfigRefine.RTM.PCY_mod,rtmConfigRefine.RTM.PCZ_mod);
        catch
            PC(:,pcaNum)=InputUser.PC
        end
    end
    PC = repmat(InputUser.PC,numPatterns)
    Patterns_export = VMOutput.PCA_VM_Pat_n(:,:,1:numPatterns);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns", ...
                tileNum), single(Patterns_export),"single");
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCStatic", ...
                tileNum), single(Patterns_export),"single");
    % for pcaNum =1:n
%%
end
