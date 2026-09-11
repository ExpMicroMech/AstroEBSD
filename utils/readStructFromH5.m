function DataStruct = readStructFromH5(AnalysisHDF5File,varargin)
    % Get HDF5 structure information
    
    h5structure = h5info(AnalysisHDF5File);
    
    % Initialize output structure
    DataStruct = struct();
    
    % Iterate over all groups in the HDF5 file
    for groupNumber = 1:numel(h5structure.Groups)
        groupInfo = h5structure.Groups(groupNumber);

        % DataStruct.(fName) = h5read(fileName, datasetReadName);
        
        % Read data from this group
        % DataStruct(groupNumber).GroupName = groupInfo.Name;
         datasetName= groupInfo.Name;

        % Clean dataset name
        fName = regexprep(datasetName, ' |\.|-|,|:|%|~|#|', '_'); 
        fName = regexprep(fName, '\(.*\)', '');
        fName = strsplit(fName,'/');
        fName = fName{end};
        DataStruct.(fName) = readDatasetsFromGroup(AnalysisHDF5File, groupInfo,varargin{:});
        %%
    end

end

function datasetsStruct = readDatasetsFromGroup(fileName, groupInfo,varargin)
    skipParams = 0;
    skipPatterns = 0;
    for i = numel(varargin)
        if strcmp(varargin{i}, "skipPatterns")
            skipPatterns = 1;
        end
        if strcmp(varargin{i},"skipParams")
            skipParams = 1;
        end
    end


    % Initialize a structure to hold dataset information
    datasetsStruct = struct();
    
    for numDataset = 1:numel(groupInfo.Datasets)
       %%
        datasetName= groupInfo.Datasets(numDataset).Name;

        % Clean dataset name
        fName = regexprep(datasetName, ' |-|,|:|%|~|#|', '_'); 
        fName = regexprep(fName, '\(.*\)', '');  % Remove everything in parentheses
        if skipPatterns==1
            if contains(fName, "Pattern")
                continue;
            end
        end
        DataSize = groupInfo.Datasets(numDataset).Dataspace.Size;
        
        try
            % try
                % warning(sprintf('%s is not a proper fieldname, skipping!'),fName)
        % end
            datasetsStruct.(fName) = zeros(DataSize); 
            datasetReadName = sprintf("%s/%s", groupInfo.Name, datasetName);
            datasetsStruct.(fName) = h5read(fileName, datasetReadName);
           

        catch ME
             try
           %%
                 warning(sprintf('%s is not a proper fieldname, skipping!'),fName)
             % datasetsStruct((fName)) = zeros(DataSize); 
            datasetReadName = sprintf("%s/%s", groupInfo.Name, datasetName);
            fName2 = regexprep(groupInfo.Name, ' |-|,|:|%|~|#|\.|', '');
            fName2 = regexprep(fName2, '\(.*\)', '');
            fName2 = strsplit(fName2,'/');
            fName2 = fName2{end};%(groupInfo.Name)            
            % 
            if numel(datasetsStruct)==0
                datasetsStruct.(fName2)=""
            end
                readVal = h5read(fileName, datasetReadName);
            
            datasetsStruct(str2num(fName)).(fName2) = readVal;
           


            % Assign empty if there's an error
             catch  
                 % datasetsStruct.(fName) = [];
                    warning('Error reading dataset "%s" in group "%s": %s', fName, groupInfo.Name, ME.message);
            end
        end
        
    end
    

    for groupNumber= 1:numel(groupInfo.Groups)
      %%
      
        groupInfo2 = groupInfo.Groups(groupNumber);
       
        fName = regexprep(groupInfo2.Name, ' |-|,|:|%|~|#|\.|', '');
        fName = regexprep(fName, '\(.*\)', '');
        fName = strsplit(fName,'/');
        fName = fName{end};
        % disp(fName)
        if ~isnan(str2double(fName))
            fName = sprintf("t%s",fName);
        end
        %%
        if skipPatterns==1
            if contains(groupInfo2.Name, "Patterns")
                continue;
            end
        end
        if skipParams==1
            if contains(groupInfo2.Name, "params")
                continue;
            end
        end
        % ':)'
        % Read data from this group
        % datasetsStruct.(fName).GroupName = groupInfo2.Name;
 % % DataStruct(groupNumber).GroupName = groupInfo.Name;
        % groupName = groupInfo.Datasets(numGroups).Name;
        
        datasetsStruct.(fName) = readDatasetsFromGroup(fileName, groupInfo2,varargin);
    end
    %%
end

