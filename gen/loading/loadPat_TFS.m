function [data_out] = loadPat_TFS(TFS_DataLoc,pattern_array,pat_type)
%LOADPAT_TFS Load processed patterns from TFS xTalView Data (works with TruePix data)
%
% [patterns_out] = loadPat_TFS(maplocation,pattern_array,ebsd_header)
%
% maplocation = file location for the data
% - should a folder that has a sub directory 'patterns'
%
% pattern_array - a [xvals, yvals] column vector,
% with 0 based indexing (TFS indexing)
%
% pat_type = '1', '2', '3'
% determine the type of pattern you want to read
%
% 1 = processed
% 2 = raw
% 3 = raw + processed

LoadProcsOpts.patternSize = 256;
LoadProcsOpts.magicNumber = LoadProcsOpts.patternSize * LoadProcsOpts.patternSize;
data_out=struct;

x_locations=pattern_array(:,1);
y_locations=pattern_array(:,2);
num_to_read=size(pattern_array,1); %number of patterns

if pat_type == 1 || pat_type == 3

    LoadProcsOpts.path_patterns = fullfile(TFS_DataLoc, 'patterns/');
    LoadProcsOpts.path_outputs = fullfile(TFS_DataLoc, 'patterns/');


    data_proc=zeros(LoadProcsOpts.patternSize,LoadProcsOpts.patternSize,num_to_read);
    %
    % x_locations=pattern_array(:,1);
    % y_locations=pattern_array(:,2);

    rows_unique=unique(y_locations); %find the unique rows
    num_rows=numel(rows_unique);

    %pre process the array data

    for row_num=1:num_rows
        % read a row
        row_counter=rows_unique(row_num);

        %find the counter values for the X patterns
        [x_counts]=find(y_locations == row_counter);
        x_pts=x_locations(x_counts);

        %create the filename of this row
        paddedRowNumber = sprintf( '%05d', (row_counter) );
        binFileName = strcat("row_", paddedRowNumber, ".bin");
        fullBinFileName = fullfile(LoadProcsOpts.path_patterns, binFileName);

        %open the pattern stack for this row
        fileID = fopen(fullBinFileName, 'r');

        for col_counter=1:numel(x_pts)
            n=x_pts(col_counter);

            % data =  fread(fileID);

            %move to the pattern in the array, and read it
            status=fseek(fileID,(LoadProcsOpts.magicNumber)*(n),'bof');
            if status == -1
                error(['The pattern location ' int2str(n) ' was not loaded properly'])
            end

            data =  fread(fileID,[256,256]);
            data90=rot90(data);

            %take this pattern and put it into the array of patterns
            data_proc(:,:,x_counts(col_counter))=data90;
        end

        fclose(fileID);

        %rotate it to give this in the image coords properly


        % if mod(size(data,1), LoadProcsOpts.magicNumber) ~=0
        %     error("length of pattern file not dividable by pattern size");
        % end
        %
        % dataReshaped = reshape(data, [LoadProcsOpts.patternSize LoadProcsOpts.patternSize size(data,1)/LoadProcsOpts.magicNumber]);
        %
        % thisProcdPattern = rot90(flipud(uint16(dataReshaped(:,:,j))),3);

        % patterns_out(:,:,1)=thisProcdPattern;
    end

    data_out.pat_proc=data_proc;
end

if pat_type == 2 || pat_type == 3
    dataType = 'uint16'; datastep=2;
    % dataType = 'uint32'; datastep=4;

    data_raw=zeros(LoadProcsOpts.patternSize,LoadProcsOpts.patternSize,num_to_read,dataType);
    % cd(LoadRawOpts.path_raw_patterns);
    LoadProcsOpts.path_raw = fullfile(TFS_DataLoc, 'raw/');
    % LoadProcsOpts.path_outputs = fullfile(TFS_DataLoc, 'patterns/');

    LoadRawOpts.infList = ls([LoadProcsOpts.path_raw '*.inf']);
    % LoadRawOpts.pakList = ls([LoadProcsOpts.path_raw '*.pak']);
    LoadRawOpts.ff = ls([LoadProcsOpts.path_raw '*.f32']);
    %read the inf list to get a list locations

    % PakIndices=[];

    %read the locations of all the points in an inf file
    for i=1:size(LoadRawOpts.infList,1)
        %read the inf file
        thisInfFile = fullfile(LoadProcsOpts.path_raw,LoadRawOpts.infList(i,:));
        [thisPakRowIndices, thisPakColumnIndices] = ReadTFSRawPatInf(thisInfFile);

        %check if there are any points from the point list which match this pak

        ok_list=false(num_to_read,1);
        ok_val=zeros(num_to_read,1);
        patlist=1:num_to_read;
        patlist=patlist';
        for s=1:num_to_read
            [xi]=find(thisPakColumnIndices == x_locations(s) & thisPakRowIndices == y_locations(s));
            
                % [thisPakColumnIndices(xi),thisPakRowIndices(xi)]
            if ~isempty(xi)
                ok_val(s)=xi;
                ok_list(s)=true;
            end
        end
        ok_short=ok_val(ok_list);
        patlist_short=patlist(ok_list);

        if ~isempty(ok_short) & sum(ok_short) > 0
            %open the pakFile
            % thisPakFile = fullfile(LoadProcsOpts.path_raw,LoadRawOpts.pakList(i,:));
            thisPakFile = [thisInfFile(1:end-3) 'pak'];
            fileID = fopen(thisPakFile);
            
            for ss=1:numel(ok_short)
                % 
                status=fseek(fileID,(LoadProcsOpts.magicNumber)*datastep*(ok_short(ss)-1),'bof'); %the 4 is because the data is in uint32
                
                if status == -1
                    error(['The pattern location ' int2str(ok_short(ss)-1) ' was not loaded properly'])
                end
                
                data =  fread(fileID,[256,256],dataType);
             
                data90=rot90(data,3);
                % figure; imagesc(data90);
                data_raw(:,:,patlist_short(ss))=data90;
            end


        end
    end
    data_out.pat_raw=data_raw;

    %read the flat field
    thisFFFile = fullfile(LoadProcsOpts.path_raw,LoadRawOpts.ff);

    fileID = fopen(thisFFFile);
    data =  fread(fileID,[256,256],'single');
    data90=rot90(data,3);
    data_out.pat_ff=  data90;
end


