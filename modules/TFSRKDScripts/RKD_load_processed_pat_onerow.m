LoadProcsOpts.path_patterns = "D:\RKD_Data\2025-04-24T12_30_27_Wse2_sample2_quantofoil2 - Copy\2025-04-24T13_28_38\patterns\";
% LoadProcsOpts.path_outputs = "D:\RKD_Data\RKD_python_Scripts\Script_Test";

LoadProcsOpts.patternSize = 479;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;
Mapsize_x = 241;

mid_y = 64;
y_span = 3;
start_y = mid_y - y_span;
end_y = mid_y + y_span;

Mapsize_y = y_span * 2 + 1;

% count = 0;

all_processed_data = zeros(LoadProcsOpts.patternSize, LoadProcsOpts.patternSize, Mapsize_x * Mapsize_y);

% for i=1:Mapsize_y
for i=start_y:end_y
    paddedRowNumber = sprintf( '%05d', (i-1) );
    binFileName = strcat("row_", paddedRowNumber, ".bin");
    fullBinFileName = fullfile(LoadProcsOpts.path_patterns, binFileName);
    
    currentRowMinIndex = 1 + (i-start_y) * Mapsize_x;
    currentRowMaxIndex = (i-start_y + 1) * Mapsize_x;

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