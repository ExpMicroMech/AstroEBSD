function [single_pattern,row_patterns] = rEBSD_TFS_Cor(y_val,x_val,TFS_DataLoc)
%read patterns from TFS data - will give you back a single pattern and the row values

% pattern_file='row_00000.bin';
pattern_file=['row_' sprintf('%05i',y_val) '.bin']; %matlab is 1 index, TFS is 0 index

pattern_full=fullfile(TFS_DataLoc,'\patterns',pattern_file);

% Open the file for reading in binary mode
fileID = fopen(pattern_full, 'rb');

% Read data from the file as unsigned 8-bit integers
data = fread(fileID, inf, 'uint8');

% Close the file
fclose(fileID);

num_pats = length(data) / 229441; % Calculate the 
row_patterns = rot90(reshape(data, [479, 479, num_pats]));

single_pattern=(row_patterns(:,:,x_val+1));
end