%% EBSD_load_processed
% Example script of loading processed patterns from xTalView to a MATLAB
% workspace variable
% Tianbi Zhang and T. Ben Britton, April 2025

%% Directories
InputUser.HDF5_folder = 'D:\ApreoTruePixData\Ni_Mapping\2025-04-01T12_15_51';

LoadProcsOpts.path_patterns = fullfile(InputUser.HDF5_folder, 'patterns/');
LoadProcsOpts.path_outputs = fullfile(InputUser.HDF5_folder, 'patterns/');

LoadProcsOpts.patternSize = 256;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;

TFSJson = ReadTFSJson(InputUser.HDF5_folder, 'site.json');

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height; % change this
% Mapsize_y = 1 % you can edit this to only load first x rows.

pattern_number = 1;

disp('Pattern conversion started');

all_procd_data = zeros(LoadProcsOpts.patternSize, LoadProcsOpts.patternSize, Mapsize_y * Mapsize_x);

%%
% each row has a .bin file which stores the patterns
for i=1:Mapsize_y
    paddedRowNumber = sprintf( '%05d', (i-1) );
    binFileName = strcat("row_", paddedRowNumber, ".bin");
    fullBinFileName = fullfile(LoadProcsOpts.path_patterns, binFileName);
    
    currentRowMinIndex = i + (i-1) * Mapsize_x;
    currentRowMaxIndex = i * Mapsize_x;

    fileID = fopen(fullBinFileName, 'r');
    data =  fread(fileID);
    fclose(fileID);

    if mod(size(data,1), LoadProcsOpts.magicNumber) ~=0
        error("length of pattern file not dividable by pattern size");
    end
    
    dataReshaped = reshape(data, [LoadProcsOpts.patternSize LoadProcsOpts.patternSize size(data,1)/LoadProcsOpts.magicNumber]);
    
    for j=1:Mapsize_x
    
    thisProcdPattern = rot90(flipud(uint16(dataReshaped(:,:,j))),3);
    EBSDPat_re = imresize(thisProcdPattern,1);
    EBSDPat_re = uint8(EBSDPat_re);
        
    all_procd_data(:,:,pattern_number) = EBSDPat_re;

    pattern_number = pattern_number + 1;

    end
    
end

disp('Pattern conversion completed');


%% View a pattern

figure;
imagesc(all_procd_data(:,:,i)); axis image; colormap('gray'); axis off;
