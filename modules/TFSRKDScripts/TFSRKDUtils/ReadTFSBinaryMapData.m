function mapData = ReadTFSBinaryMapData(mapPath, mapName, TFSJson, dataType)
%ReadTFSBinaryMap reads a map data from a TFS EBSD binary file (usually .bin)
% and converts it automatically to the dimension of the map.
%
% Inputs:
% mapPath: folder path to the map
% mapName: file name for the map (e.g. counts.bin)
% TFSJson: the metadata of the current EBSD project, deserialized 
% from "site.json" of the project. This can be obtained
% using the ReadTFSJson() function. Map width and height are extracted.
%
% dataType: type of data in the binary file, e.g. 'int16', 'int32'. If 
% not specified, it defaults to single precision float point number. A
% warning message will be printed if the wrong datatype is used, which will
% cause a mismatch of the size between the data read from the binary file
% and the physical size of the map.
%
% Outputs:
% mapData: reshaped map data. An empty object is returned if the wrong
% datatype is provided.


if ~exist('dataType', 'var')
    dataType = 'single';
    warning("No data format specified - default to single precision.")
end

fileID = fopen(fullfile(mapPath, mapName));
mapVectorData = fread(fileID, dataType);

try
mapData = reshape(mapVectorData, [TFSJson.Info.Width TFSJson.Info.Height])';
catch
warning("Size mismatch - your input datatype is likely wrong.")
mapData = [];
end

end