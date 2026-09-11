function outputs = getFHKLAndD(InputUser,cutoff,pNum)
    % InputUser.
    [m_unique,f_hkl] = calcStructureFactorNonCorr(InputUser.Phase_Input(pNum),InputUser.maxHKL,InputUser);
    %%
    dSpacing = m_unique.dspacing;
    FDotD = f_hkl.*dSpacing;
    [FDotDSorted,idsFHKL] = sort(abs(FDotD),'descend');
    m_uniqueSorted = m_unique(idsFHKL);
    f_hklUniqueSorted = f_hkl(idsFHKL);
    dSorted = dSpacing(idsFHKL);
    %%
    MillerChar = [];
    for i = 1:size(m_uniqueSorted)
        MillerChar = [MillerChar, sprintf("%s",char(m_uniqueSorted(i)))];
    end
    MillerChar = MillerChar';
    selection = FDotDSorted>cutoff*max(FDotDSorted(:));
    OutputTableFHKL = table(MillerChar(selection), abs(f_hklUniqueSorted(selection)), FDotDSorted(selection),dSorted(selection),...
        'VariableNames',{'Miller','F_HKL','F_HKL * D', 'D'} )
    writetable(OutputTableFHKL, fullfile(InputUser.FigureFolder , ...
                                    sprintf('%s_FHKLList.csv',InputUser.Phase_Input{pNum}) ) ) ;
    writetable(OutputTableFHKL, fullfile(InputUser.FigureFolder , ...
                                    sprintf('%s_FHKLList.xls',InputUser.Phase_Input{pNum}) ) ); 
    
    %%
    m_unique_norm = normalize(m_uniqueSorted);
    Crossproducts = norm(cross(m_unique_norm, m_unique_norm'));
    indicesFHKLCleaned = zeros(size(m_uniqueSorted));
    
    for i = 1: size(m_uniqueSorted,1)
        equalMillerIndices = find(Crossproducts(i,:)==0);
        % [~,maxIndex] = max(f_hklSorted(equalMillerIndices));
        % because we sorted before, we can just select the first one, as it
        % will be the top most
        indicesFHKLCleaned(i) = equalMillerIndices(1); 
    end
    %%
    mCleaned = m_uniqueSorted(unique(indicesFHKLCleaned))
    fhklCleaned = f_hklUniqueSorted(unique(indicesFHKLCleaned))
    FDotDSortedCleaned = FDotDSorted(unique(indicesFHKLCleaned))
    dSortedCleaned = dSorted(unique(indicesFHKLCleaned))
    %%
    % figure
    % imagesc(Crossproducts==0)
    % colorbar
    %%
    MillerChar = [];
    
    selection = FDotDSortedCleaned>cutoff*max(FDotDSortedCleaned(:));
    mSelected  = mCleaned(selection)
    for i = 1:size(mSelected)
        MillerChar = [MillerChar, sprintf("%s",char(mSelected(i)))];
    end
    MillerChar = MillerChar';
    OutputTableFHKLCleaned = table(MillerChar(selection), abs(fhklCleaned(selection)),FDotDSortedCleaned(selection),dSortedCleaned(selection), 100*FDotDSortedCleaned(selection)/max(FDotDSortedCleaned(:)),...
        'VariableNames',{'Miller','F_HKL','F_HKL * D', 'D','% max(F_hkl*d)'} )
    writetable(OutputTableFHKLCleaned, fullfile(InputUser.FigureFolder , ...
                                    sprintf('%s_FHKLList_Cleaned.csv',InputUser.Phase_Input{pNum}) ) ) ;
    writetable(OutputTableFHKLCleaned, fullfile(InputUser.FigureFolder , ...
                                    sprintf('%s_FHKLList_Cleaned.xls',InputUser.Phase_Input{pNum}) ) );
    outputs.miller = mSelected;
    outputs.fhkl = fhklCleaned(selection);
    outputs.FDotD = FDotDSortedCleaned(selection);
    outputs.dSpacing = dSorted(selection);
end