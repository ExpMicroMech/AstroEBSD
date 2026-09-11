function getOrisFullMap(rtmConfigRefine,InputUser,MapInfo,MicroscopeData)
    numPatterns = numel(MapInfo.MapData.XBeam)
    tileNum=1
    % MapInfo = tile.MapInfo;
    % tileNum =tile.tileNum; %hardcoded externally
    % [YSearchGrid ,XSearchGrid ]= meshgrid(tile.rowstart(tileNum):1:tile.rowfin(tileNum),tile.colstart(tileNum):1:tile.colfin(tileNum));
    % Plotting starts here)
    % numPatterns = size(numPatterns)
    % tileNum = 1; % handcarfte, thus VMOutput only contains 1 tile anyways
    numPhasesEBSD = unique(MapInfo.EBSPData.ebsd_map.phase);
    gmatrix_conv = zeros(3,3,numPatterns);
    gmatrix_in = zeros(3,3,numPatterns);
    PC = zeros(3,numPatterns);
    Phase = zeros(1,numPatterns);
    eangs_in_export = zeros(3,numPatterns);
    try
        [PCX,PCY,PCZ] = getPatternCentreProjective(MapInfo.MapData.XSampleCor,MapInfo.MapData.YSampleCor,...
        rtmConfigRefine.RTM.ProjectiveTransform );
        PC =[PCX(:)';PCY(:)';PCZ(:)']
    catch
        warning('loading PC from map data')
        try
            PC = [MapInfo.MapData.pcx'; MapInfo.MapData.pcy';MapInfo.MapData.pcz'];
        catch
            PC = repmat(InputUser.PC',1,numPatterns);
            warning(' PC from mapdata could not be loaed. setting static PC.')
        end
    end
%%
    rtmConfigRefine.Detector_tilt = rtmConfigRefine.RTM.Rx(MicroscopeData.TotalTilt);
    % for numPhasesEBSD
    EBSD_Map = MapInfo.EBSPData.ebsd_map;
    for pNum = 0: numPhasesEBSD
        comparisonBool = EBSD_Map.phase==pNum;
        EBSD_MapLoc = EBSD_Map(comparisonBool);
        ori =  EBSD_MapLoc.orientations;
        gmatrix_in(:,:,comparisonBool) =ori.matrix;
        eangs_in_export(1,comparisonBool) = ori.phi1;
        eangs_in_export(2,comparisonBool) = ori.Phi;
        eangs_in_export(3,comparisonBool) = ori.phi2;
        % Patterns_export(:,:,pcaNum) = VMOutput.PCA_VM_Pat_n(:,:,pcaNum);
        gmatrix =gmatrix_in(:,:,comparisonBool);
        gmatrix_conv(:,:,comparisonBool) = pagemtimes( pagemtimes( pagetranspose(gmatrix),rtmConfigRefine.Detector_tilt) ...
                                            , InputUser.EdaxConversionMatrix )  ;


    end
    %%
    % Patterns_export = tile.(tileName).PCA_VM_Pat_n(:,:,1:numPatterns);
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
    % writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns", ...
    %             tileNum), single(Patterns_export),"single");
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/EangsCommercial", ...
               tileNum), eangs_in_export);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercialDetector", ...
                tileNum), gmatrix_conv);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
                    tileNum), PC);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/Phase", ...
                    tileNum), Phase,"uint16");
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/XVals", ...
                    tileNum), MapInfo.MapData.XSampleCor);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/YVals", ...
                    tileNum), MapInfo.MapData.YSampleCor);
end
%%
% end
