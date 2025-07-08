LoadProcsOpts.path_patterns = fullfile(InputUser.HDF5_folder, 'patterns/');
% LoadProcsOpts.path_outputs = "D:\RKD_Data\RKD_python_Scripts\Script_Test";

LoadProcsOpts.patternSize = 479;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height;
% Mapsize_y = 1;

pattern_number = 1;

disp('Pattern conversion started');


% each row has a .bin file which stores the patterns
for i=1:Mapsize_y
    paddedRowNumber = sprintf( '%05d', (i-1) );
    binFileName = strcat("row_", paddedRowNumber, ".bin");
    fullBinFileName = LoadProcsOpts.path_patterns + binFileName;
    
    currentRowMinIndex = i + (i-1) * Mapsize_x;
    currentRowMaxIndex = i * Mapsize_x;

    fileID = fopen(fullBinFileName);
    data =  fread(fileID);

    if mod(size(data), LoadProcsOpts.magicNumber) ~=0
        error("length of pattern file not dividable by 229441");
    end
    
    dataReshaped = reshape(data, [LoadProcsOpts.patternSize LoadProcsOpts.patternSize size(data,1)/LoadProcsOpts.magicNumber]);
    dataReshaped = rot90(dataReshaped);
    fclose(fileID);

    for j=1:Mapsize_x
    % all_processed_data(:,:,[currentRowMinIndex:currentRowMaxIndex]) = dataReshaped;
    
    thisProcdPattern = uint16(dataReshaped(:,:,j));
    EBSDPat_re = imresize(thisProcdPattern,1/binning);
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

% clear data dataReshaped;
% clear LoadProcsOpts;9