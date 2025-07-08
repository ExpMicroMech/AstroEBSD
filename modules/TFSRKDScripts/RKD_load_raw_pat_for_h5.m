% LoadRawOpts.path_raw_patterns = fullfile(InputUser.HDF5_folder, 'raw/2023-11-16T13_36_30.392241091');
% EBSD patterns have different dimensions, refer to

LoadRawOpts.quadrantSize = 256;
LoadRawOpts.quadrantpatternSize = LoadRawOpts.quadrantSize ^ 2;
LoadRawOpts.fullRawPatternSize = LoadRawOpts.quadrantSize ^ 2 * 4;

LoadRawOpts.edgeoffset = 35;

LoadRawOpts.rawPatternSize = 549;
LoadRawOpts.patternSize = 479;
% LoadRawOpts.magicNumber = LoadRawOpts.patternSize * LoadRawOpts.patternSize;

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height;

if exist('binning') == 0
    binning = 1;
end

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

clear infList 

%% read from .pak

pattern_number = 1;

disp('Pattern conversion started');
disp(['Loading from ' num2str(LoadRawOpts.totalPakFiles) ' raw files']);

for i=1:size(LoadRawOpts.pakList,1)

thisPakFile = LoadRawOpts.pakList(i,:);
thisInfFile = LoadRawOpts.infList(i,:);

[thisPakRowIndices, thisPakColumnIndices] = ReadTFSRawPatInf(thisInfFile);

% read the big thing only once per file!
fileID = fopen(thisPakFile);
all_x = fread(fileID, dataType); % pay attention to this- might be different!!!
fclose(fileID);

thisPakIndices = thisPakColumnIndices + thisPakRowIndices .* Mapsize_x + 1;

    for j = 1:size(thisPakIndices,1)
          
        x_slice = all_x([1:LoadRawOpts.fullRawPatternSize] + (j-1) * LoadRawOpts.fullRawPatternSize );

        EBSDPat_1 = TileQuadPatterns(x_slice);

        EBSDPat_re=uint16(EBSDPat_1);
        
        if binning ~= 1
            EBSDPat_re = imresize(EBSDPat_1,1/binning);
        else
            EBSDPat_re = imresize(EBSDPat_1,1);
        end

        EBSDPat_re=uint16(EBSDPat_re);

        if pattern_number == 1
        EBSPData.PW=double(size(EBSDPat_re,2));
        EBSPData.PH=double(size(EBSDPat_re,1));

        h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',EBSPData.PH);
        h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',EBSPData.PW);

        EBSPData.made=0;
        else
        EBSPData.made=1;
        end
        
        bWriteEBSP(EBSPData, thisPakIndices(j), EBSDPat_re, 'uint16')

        % allRawPatterns(:,:,thisPakIndices(j)) = EBSDPat_re;

        if mod(thisPakIndices(j), 1000) == 0
            disp(['Pattern ' int2str(thisPakIndices(j)) ' of ' int2str(Mapsize_full) ' patterns converted']);
        end

    end

end

disp('Pattern conversion finished')

clear x_slice all_x