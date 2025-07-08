%% EBSD_load_processed_for_h5.m
% Example script of loading processed patterns from xTalView to a HDF (.h5)
% file.
% Tianbi Zhang and T. Ben Britton, April 2025

% N.B. This script is called by bCreateHDF5FromTFSData.m and should not be
% run on its own.

%% Directories

LoadProcsOpts.patternSize = 256;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height;

pattern_number = 1;

disp('Pattern conversion started');

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
    
    thisProcdPattern = flipud(uint16(dataReshaped(:,:,j)))';
    EBSDPat_re = imresize(thisProcdPattern,1/binning);
    EBSDPat_re = imresize(thisProcdPattern,1);
    EBSDPat_re = uint16(EBSDPat_re);

        if pattern_number == 1
            EBSPData.made=0;
            
            EBSPData.PW=double(size(EBSDPat_re,2));
            EBSPData.PH=double(size(EBSDPat_re,1));

            h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',EBSPData.PH);
            h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',EBSPData.PW);
        else
            EBSPData.made=1;
            
        end
        bWriteEBSP(EBSPData,pattern_number,EBSDPat_re, 'uint16')
        
        if 1000*round(pattern_number/1000) == pattern_number
            disp(['Pattern ' int2str(pattern_number) ' of ' int2str(EBSPData.numpats) ' patterns converted']);
        end

        pattern_number = pattern_number + 1;

    end
    
end

disp('Pattern conversion completed');


clear data dataReshaped;
clear LoadProcsOpts;
