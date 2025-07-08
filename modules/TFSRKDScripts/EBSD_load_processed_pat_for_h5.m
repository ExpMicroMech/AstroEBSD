% LoadProcsOpts.path_patterns = fullfile(InputUser.HDF5_folder, 'patterns/');
% LoadProcsOpts.path_outputs = "D:\RKD_Data\RKD_python_Scripts\Script_Test";

LoadProcsOpts.patternSize = 256;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height;
% Mapsize_y = 20;

pattern_number = 1;

disp('Pattern conversion started');

% all_procd_data = zeros(LoadProcsOpts.patternSize, LoadProcsOpts.patternSize, Mapsize_y * Mapsize_x);

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
    % EBSDPat_re = imresize(thisProcdPattern,1);
    EBSDPat_re = flipud(uint16(EBSDPat_re));

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
        
        % all_procd_data(:,:,pattern_number) = EBSDPat_re;

        if 1000*round(pattern_number/1000) == pattern_number
            disp(['Pattern ' int2str(pattern_number) ' of ' int2str(EBSPData.numpats) ' patterns converted']);
        end

        pattern_number = pattern_number + 1;

    end
    
end

disp('Pattern conversion completed');


% clear data dataReshaped;
% clear LoadProcsOpts;9

%%
% pat2 = all_procd_data(:,:,744);
% 
% figure;
% imagesc(pat2); axis xy; axis image; colormap('gray');