% raw patterns are more annoying to deal with since the packaging is
% different.

% LoadRawOpts.path_raw_patterns = fullfile(InputUser.HDF5_folder, 'raw/');
LoadRawOpts.path_raw_patterns = fullfile("D:\RKD_Data\2024-11-04T16_22_44\raw\2024-11-04T14_32_19.743335136");

TFSJson = ReadTFSJson('D:\RKD_Data\2024-11-04T16_22_44', 'site.json');

LoadRawOpts.patternSize = 256;
LoadRawOpts.fullRawPatternSize = LoadRawOpts.patternSize * LoadRawOpts.patternSize;

Mapsize_x = TFSJson.Info.Width;
Mapsize_y = TFSJson.Info.Height;
Mapsize_full = TFSJson.Info.TotalNumberOfPoints;

% pre-allocate an array for the patterns
allRawPatterns = zeros(LoadRawOpts.patternSize, LoadRawOpts.patternSize, Mapsize_full);

%%
cd(LoadRawOpts.path_raw_patterns);
LoadRawOpts.infList = ls('*.inf');
LoadRawOpts.pakList = ls('*.pak');

LoadRawOpts.totalPakFiles = size(LoadRawOpts.pakList,1);

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

        EBSDPat_1 = reshape(x_slice, [LoadRawOpts.patternSize, LoadRawOpts.patternSize])';

        EBSDPat_re=uint16(EBSDPat_1);
        
        allRawPatterns(:,:,thisPakIndices(j)) = EBSDPat_re;

        if mod(thisPakIndices(j), 1000) == 0
            disp(['Pattern ' int2str(thisPakIndices(j)) ' of ' int2str(Mapsize_full) ' patterns converted']);
        end

    end

end

disp('Pattern conversion finished')

clear x_slice all_x