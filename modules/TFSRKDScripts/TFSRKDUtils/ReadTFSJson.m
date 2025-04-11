%% this function is simply a json reader

function [TfsStruct] = ReadTFSJson(path, file)
jsonFile = fullfile(path, file);

jsonText = fileread(jsonFile);

TfsStruct = jsondecode(jsonText);

end