function [pcx,pcy, pcz] = patternCentreGeometric(calibration,MapData)
%     MapData = MapInfo.MapData
%     calibration = InputUser.calibration;
    % detPixelSize = 55;
    % detNumPixels = 512;
    detSize = calibration.DetPXSize * calibration.DetNumPixels;
    % we assume no sample roation only til
    pcy = -sin(calibration.SampleTilt*degree)*MapData.YSampleCor/detSize+calibration.pcy;
    pcx = -MapData.XSampleCor/detSize+calibration.pcx;
    pcz =  cos(calibration.SampleTilt*degree)*MapData.YSampleCor/detSize+calibration.pcz;
end