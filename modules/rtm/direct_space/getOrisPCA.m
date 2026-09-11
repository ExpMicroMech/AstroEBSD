function getOrisPCA(tile)
    numPatterns =VMOutput.NumPCA
    % Plotting starts here)
    % numPatterns = size(numPatterns)
    for pcaNum =1: numPatterns
        selection = VMOutput.filtered_map'==pcaNum;
        
        [pnum,~]=find(MapInfo.MapData.XBeam==XSearchGrid(selection)' & MapInfo.MapData.YBeam==YSearchGrid(selection)');
        [pnumOrg,~]=find(MapInfo.MapData.XBeam==XSearchGrid(:)' & MapInfo.MapData.YBeam==YSearchGrid(:)');
        ebsdMapSelection = MapInfo.EBSPData.ebsd_map(pnum);
        ebsdMapSelection = ebsdMapSelection(ebsdMapSelection.phase==1);
        pnumNew = ebsdMapSelection.id;

        meanOri = mean(ebsdMapSelection.orientations)
        meanX = mean(MapInfo.MapData.XBeam(pnumNew))
        meanY = mean(MapInfo.MapData.YBeam(pnumNew))
        Detector_tilt = rtmConfigRefine.RTM.Rx(MicroscopeData.TotalTilt);
        gmatrix_in=meanOri.matrix;
        gmatrix = gmatrix_in(:,:,i);
        gmatrix = gmatrix'*Detector_tilt*InputUser.EdaxConversionMatrix ;
        rtmConfigRefine.Settings_PCin.start  = getPCFromPlane(meanX,meanY,PCX_mod,PCY_mod,PCZ_mod)
    end
end