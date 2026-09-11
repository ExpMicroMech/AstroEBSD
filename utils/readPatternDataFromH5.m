function PatternData=readPatternDataFromH5(AnalysisHDF5File,varargin)

skipPatterns = 0;
for i = numel(varargin)
    if strcmp(varargin{i}, "skipPatterns")
        skipPatterns = 1;
    end
end
h5strucutre = h5info(AnalysisHDF5File);%,'/Patterns','min')
%%
for grNumber = 1:size(h5strucutre.Groups)
    
    if strcmp(h5strucutre.Groups(grNumber).Name, '/Patterns')
        PatternGroup = grNumber;
        break;
    end
end
%%
patternPaths = {};
patternPaths = findGroupPaths(h5strucutre.Groups(PatternGroup), '/','Patterns',patternPaths);
patternPaths = patternPaths(2:end); %drop the first one, as this is the "super group"

% for patternPath =1 : numel(patternPaths{1})
    % DataInfo = h5info(AnalysisHDF5File,patternPaths{1});
% end


% DataName = '';

%%
PatternData = struct();
% 
try
    % This is hard coded for now. In the future iterate over the h5
    % iteratively.
    PatternData.InputUser.X =   h5read(AnalysisHDF5File,'/UserInfo/XGrid');
    PatternData.InputUser.Y =   h5read(AnalysisHDF5File,'/UserInfo/YGrid');
catch
    warning('X and Y grid not read')
end
for groupNumbers = 1 :numel(patternPaths)
    DataInfo = h5info(AnalysisHDF5File,patternPaths{groupNumbers});
    % NumberGroups = size(DataInfo.Groups);

   
    DataName =[];
    
    patternPathSplit = strsplit(patternPaths{groupNumbers},'/');
   
    for i = 3:numel(patternPathSplit)-1
        DataName = [DataName  '_'  patternPathSplit{i}];
    end
    PatternData(groupNumbers).DataName = DataName;
    for numDataset = 1:size(DataInfo.Datasets)
        datasetName= DataInfo.Datasets(numDataset).Name;
        fName = regexprep(DataInfo.Datasets(numDataset).Name,' |-|,|:|%|~|#| |','_');
        fName = regexprep(fName,'\(.*\)','');
        if skipPatterns==1
            if contains(fName, "Pattern")
                continue;
            end
        end
        DataSize = DataInfo.Datasets(numDataset).Dataspace.Size;
        try
            PatternData(groupNumbers).(fName) = zeros(DataSize(1),DataSize(2),DataSize(3));
        catch
            PatternData(groupNumbers).(fName) = zeros(DataSize(1),DataSize(2));
        end
        % for pNum = 1:NumberGroups(1)
            try

                %%
                datasetReadName = sprintf("%s/%s",DataInfo.Name,datasetName);
                % pNum;

                
                PatternData(groupNumbers).(fName)= h5read(AnalysisHDF5File,datasetReadName);
            %%
            catch
                warning(sprintf('empty DAtaset detected in %s , %s', DataName, fName))
            end
        % end

    end
end
% Fi
    
end
    function patternPaths = findGroupPaths(groupInfo,currentPath,patternName,patternPaths)
    %% Update the current path with the group name
    if isempty(currentPath)
        currentPath = [currentPath '/' groupInfo.Name];
    else
        currentPath = groupInfo.Name;
    end
    %
    % Check if this group is named 'Pattern'
    currentGroupName = strsplit(groupInfo.Name,'/');
    currentGroupName = currentGroupName{end};

    if strcmp(currentGroupName, [patternName])
        patternPaths{end+1} = currentPath; % Add path to results
        fprintf('Found Pattern Group: %s\n', currentPath); % Display found path
    end
    
    % Loop through all subgroups in this group recursively
    for i = 1:length(groupInfo.Groups)
        subgroup = groupInfo.Groups(i);
        try
        patternPaths = findGroupPaths(subgroup,currentPath,patternName,patternPaths); % Recursive call
        catch
            continue
        end
    end
end
