LoadProcsOpts.path_patterns = "D:\RKD_Data\2023-11-16T12_36_34_WSe2_orig\2023-11-16T13_36_03\patterns\";
LoadProcsOpts.path_outputs = "D:\RKD_Data\RKD_python_Scripts\Script_Test";

LoadProcsOpts.patternSize = 479;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;
Mapsize_x = 100;
Mapsize_y = 45;

count = 0;

all_processed_data = zeros(LoadProcsOpts.patternSize, LoadProcsOpts.patternSize, Mapsize_x * Mapsize_y);

for i=1:Mapsize_y
    paddedRowNumber = sprintf( '%05d', (i-1) );
    binFileName = strcat("row_", paddedRowNumber, ".bin");
    fullBinFileName = fullfile(LoadProcsOpts.path_patterns, binFileName);
    
    currentRowMinIndex = 1 + (i-1) * Mapsize_x;
    currentRowMaxIndex = i * Mapsize_x;

    fileID = fopen(fullBinFileName);
    data =  fread(fileID);

    if mod(size(data), LoadProcsOpts.magicNumber) ~=0
        error("length of pattern file not dividable by size of the pattern");
    end
    
    dataReshaped = reshape(data, [LoadProcsOpts.patternSize LoadProcsOpts.patternSize size(data,1)/LoadProcsOpts.magicNumber]);
    dataReshaped = rot90(dataReshaped);
    all_processed_data(:,:,[currentRowMinIndex:currentRowMaxIndex]) = dataReshaped;
    fclose(fileID);
end

clear data dataReshaped;
clear LoadProcsOpts;