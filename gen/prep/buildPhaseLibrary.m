function Phases=buildPhaseLibrary(phaseFolder)
    phaseFiles = dir(fullfile(phaseFolder,'phasefiles','*.pha'))
    phaseNames = {}
    Phases = struct();
    for phaseNum  = 1:numel(phaseFiles)
        phaseNames{end+1}= strrep(phaseFiles(phaseNum).name,".pha","");
        phaseFiles(phaseNum).name;
        [phaseKeys, phaseVals] = readAllPhaKeys(fullfile(phaseFiles(phaseNum).folder,phaseFiles(phaseNum).name));
        try
        for keyNum = 1: numel(phaseKeys)
            fieldName = strrep(phaseKeys(keyNum),'$',"");
            Phases(phaseNum).(fieldName) = phaseVals(keyNum);
        end
        catch
            warning('invalid field in Phase: %s,\n field: %s \n skiping phase',...
                phaseFiles(phaseNum).name,fieldName)
            
        end
    end

    %% build table
    fNames = fieldnames(Phases); 
    fTypes = repmat(["string"],numel(fNames),1);
    allPhaseTable = table('Size',[numel(Phases),numel(fNames)],'VariableNames',fNames,'VariableTypes',fTypes);
    %%
    fieldClasses = []
    for fieldNum = 1:numel(fNames)
        % fieldNum
        field = fNames{fieldNum};
        isNumeric = 0;
        isString = 0;
        isArray = 0;
        % fieldsNumeric = []
        % fieldsString = []
        for rowNum = 1:numel(Phases)
            if ~isempty(Phases(rowNum).(field))
                fieldNumeric = str2num(Phases(rowNum).(field));
                if ~isempty(fieldNumeric)
                    isNumeric = isNumeric+1;
                    if numel(fieldNumeric)>1
                        isArray = isArray+1;
                    end
                else
                    isString = isString+1;
                end
            end
        end
        %%
        if isNumeric>isString
            if isArray>=1
               allPhaseTable=convertvars(allPhaseTable,field,'cell');  %maybe not the bes for some ints. but so it is ( we can check modulo later on)  
            else
                allPhaseTable=convertvars(allPhaseTable,field,'double');
            end
            for rowNum = 1:numel(Phases)
            if ~isempty(Phases(rowNum).(field))
                try
                allPhaseTable.(field)(rowNum) = str2num(Phases(rowNum).(field));
                catch
                    try
                 allPhaseTable.(field)(rowNum)= {str2num(Phases(rowNum).(field))};
                    catch
                        warning('faulty phase file, not generating entry')
                    end
                end
            end
            end
        % end
        else
            for rowNum = 1:numel(Phases)
                if ~isempty(Phases(rowNum).(field))
                    try
                    allPhaseTable.(field)(rowNum)=Phases(rowNum).(field);
                    catch
                        try
                     allPhaseTable.(field)(rowNum)={Phases(rowNum).(field)};
                     catch
                        warning('faulty phase file, not generating entry')
                    end
                    end
                end
            end
        end
    end
% end
    %%
    % allPhaseTable
    %%
    writetable(allPhaseTable,fullfile(phaseFolder,"PhaseLibrary.dat"));
    writetable(allPhaseTable,fullfile(phaseFolder,"PhaseLibrary.txt"));
    save(fullfile(phaseFolder,"PhaseLibrary.mat"),"allPhaseTable");
    % allPhaseTable = struct2table(Phases)
    %%
    % [numPhases, RTM_info] = Phase_Builder_RTMSimple(phaseNames,phaseFolder);
%%
end