function [PatternData,Phase_Names,MaxIters,PatSizes] =GetPatternIndices_PCRefinement(PatternData)
PatSizes = [];
MaxIters = [];
Phase_Names = [];
for RefinementNumber = 1:numel(PatternData)
    patSizeMatch = regexp(PatternData(RefinementNumber).DataName, '(?<=_PatSize_)\d+', 'match');
    if ~isempty(patSizeMatch)
        PatSizeValue = str2double(patSizeMatch{1});
    else
        PatSizeValue = []; % Handle case where no match is found
    end
    
       pattern = '(?:(.*Iter_\d+_)|joubert|_\d{1,2}_|_\d{1,2}kV_|\d{4}_\d{2}_\d{2}__\d{2}_\d{2}_\d{2}|dmin-[\d.-]+|npx-\d+)';
    
    % Remove unwanted substrings
    cleanedString = regexprep(PatternData(RefinementNumber).DataName, pattern, '');
    cleanedString = regexprep(cleanedString, "ehistmin", 'kv');
    % Clean up any extra underscores or spaces left behind
    cleanedString = string(strtrim(regexprep(cleanedString, '_+', '_'))); % Replace multiple underscores with a single one

    
    maxIterMatch = regexp(PatternData(RefinementNumber).DataName, '(?<=_MaxIter_)\d+', 'match');
    if ~isempty(maxIterMatch)
        MaxIterValue = str2double(maxIterMatch{1});
    else
        MaxIterValue = []; % Handle case where no match is found
    end
    PatSizes = [PatSizes PatSizeValue];
    MaxIters = [MaxIters MaxIterValue];
    Phase_Names = [Phase_Names cleanedString];
end
[~,index] = sort(PatSizes.');
PatternData = PatternData(index);
Phase_Names = Phase_Names(index);
MaxIters = MaxIters(index);
PatSizes = (PatSizes(index));
end