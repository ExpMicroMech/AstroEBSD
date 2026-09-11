function preparePCARefineTile(tile,tileName,rtmConfigRefine,InputUser)
    numPatterns =tile.(tileName).NumPCA;
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
    % try
        % numPatterns=tile.numPatterns;
    % end
    Patterns_export = tile.(tileName).PCA_VM_Pat_n(:,:,1:numPatterns);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns", ...
                tileNum), single(Patterns_export),"single");
    % for pcaNum =1:n
%%
end
