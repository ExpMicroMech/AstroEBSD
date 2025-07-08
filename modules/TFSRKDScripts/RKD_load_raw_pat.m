LoadRawOpts.path_raw_patterns = "C:\Users\billy\Documents\MATLAB\2023-11-16T12_36_34_WSe2_orig\2023-11-16T13_36_03\raw\2023-11-16T13_36_30.392241091";
% path_patterns = "C:\Users\billy\Documents\RKD_Data\2023-11-16T12_36_34_WSe2_orig\2023-11-16T13_36_03\patterns";
% path_outputs = "C:\Users\billy\Documents\RKD_Data\RKD_python_Scripts\Script_Test";

LoadRawOpts.quadrantSize = 256;
LoadRawOpts.quadrantpatternSize = LoadRawOpts.quadrantSize ^ 2;
LoadRawOpts.fullRawPatternSize = LoadRawOpts.quadrantSize ^ 2 * 4;

LoadRawOpts.edgeoffset = 35;

LoadRawOpts.rawPatternSize = 549;
LoadRawOpts.patternSize = 479;
LoadRawOpts.magicNumber = patternSize * patternSize;

Mapsize_x = 100; % to be replaced by json
Mapsize_y = 62; % to be replaced by json

%%
cd(LoadRawOpts.path_raw_patterns);
LoadRawOpts.infList = ls('*.inf');
LoadRawOpts.pakList = ls('*.pak');

LoadRawOpts.totalPakFiles = size(LoadRawOpts.pakList,1);
LoadRawOpts.rowIndexToPakLUT = zeros(1,LoadRawOpts.totalPakFiles); % max row index for each file


for i=1:LoadRawOpts.totalPakFiles
    thisInfFile = LoadRawOpts.infList(i,:);
    
    % rowIndices = readTFSRawPatInf(thisInfFile);
    LoadRawOpts.rowIndexToPakLUT(i) = max(ReadTFSRawPatInf(thisInfFile));
    
    % the corresponding pak file has the same index as 
end

clear infList rowIndices

%% read from .pak
all_raw_data = zeros(LoadRawOpts.rawPatternSize, LoadRawOpts.rawPatternSize, Mapsize_x * Mapsize_y);

pattern_number = 1;
pakStopAt = find(LoadRawOpts.rowIndexToPakLUT <= Mapsize_y, 1, 'last' );

% for i=1:totalPakFiles
for i=1:LoadRawOpts.pakStopAt
thisPakFile = LoadRawOpts.pakList(i,:);

% read the big thing only once per file!
fileID = fopen(thisPakFile);
all_x = fread(fileID,'uint16');
fclose(fileID);

if i == 1
    thisPakMinRowIndex = 0;
else
    thisPakMinRowIndex = rowIndexToPakLUT(i-1)+1;
end
thisPakMaxRowIndex = rowIndexToPakLUT(i);

for j = thisPakMinRowIndex:thisPakMaxRowIndex

    for k = 1:(Mapsize_x)
        x_slice = all_x([1:LoadRawOpts.fullRawPatternSize]+(k-1)*LoadRawOpts.fullRawPatternSize);
        all_raw_data(:,:,pattern_number) = TileQuadPatterns(x_slice);

        pattern_number = pattern_number +1;
    end

end

end

clear x_slice

%%

% bg = mean(all_raw_data,3);
% bg = bg(36:512,36:512);
% 
%%
% pat2 = all_raw_data(36:512,36:512,4151);
% pat2cor = pat2 ./ imgaussfilt(pat2, 5);
% 
% figure;
% imagesc(pat2cor); axis xy; axis image; colormap('gray');
