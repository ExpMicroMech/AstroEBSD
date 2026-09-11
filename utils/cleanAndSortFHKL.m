function [m,I] = cleanAndSortFHKL(m,f_hkl,cutoff)
    I = f_hkl.*conj(f_hkl);
    I = abs(f_hkl);
    I_norm = I./max(I(:));
    I = I_norm(I_norm>cutoff);
    [I,indexSorted] = sort(I,'descend');
    
    m = m(I_norm>cutoff);
    m = m(indexSorted);
end