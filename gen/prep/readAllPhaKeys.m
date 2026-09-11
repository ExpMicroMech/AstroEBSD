function [keys, vals] = readAllPhaKeys(Phase_File)
    % Read all keys from the pha file, to check, what they contain
    % This way we do not overwrite
    % from Phase_Builder_RTM
    delimiter = '';
    formatSpec = '%s%[^\n\r]';
    fileID = fopen(Phase_File,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
    fclose(fileID);
    PhaseRaw = [dataArray{1:end-1}];
    array_num = 1:size(PhaseRaw,1);

    parPosition =  strfind(PhaseRaw,'$');
    parLines = find(~cellfun(@isempty,parPosition));
    keys = [];
    vals = [];
    for parNum = 1:numel(parLines)
        keys = [keys string( PhaseRaw{parLines(parNum)})];
        vals = [vals string(PhaseRaw{parLines(parNum)+1})];
    end
    %%
end