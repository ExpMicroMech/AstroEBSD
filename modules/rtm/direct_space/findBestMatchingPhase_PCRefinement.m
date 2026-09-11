function [PatternData] = findBestMatchingPhase_PCRefinement(PatternData,patSize)
    [PatternData,Phase_Names,MaxIters,PatSizes] = GetPatternIndices_PCRefinement(PatternData);
    maxNCC = 0;
    maxNCCIndex = NaN;
    PatIndex = find(PatSizes == patSize);
    %%
    for i =1:numel(PatIndex)
        meanNCC = mean(PatternData(PatIndex(i)).NCC);
       if meanNCC>maxNCC
           maxNCC = meanNCC;
           maxNCCIndex = PatIndex(i);
       end
    end
    %%
    PatternData = PatternData(maxNCCIndex)
end
