function [Map_Out]=bLoadTFS_Bin(bin_location,MapSizeY,MapSizeX)
% Loads TFS bin files to create map data that can be added to the EBSD
% containers for MTEX
%
% INPUTS - 
% bin_location - the path location of the bin file
% MapSizeY - map y size in pts
% MapSizeX - map x side in pts
% 
% OUTPUTS - 
% Map_Out = double array reshaped from the input bin file
%
% TBB 2025

if exist(bin_location) == 0
    error(['The file to load does not exist: ' bin_location)
end

    fid = fopen(bin_location, 'rb');                   % Open file for reading in binary mode
    fseek(fid, 0, 'bof');                       % Offset in bytes (0 means beginning of file)
    x = fread(fid, inf, 'single');              % Read all data as 32-bit floats (single precision)
    fclose(fid);   
    Map_Out=double(reshape(x,[MapSizeY,MapSizeX]));
end