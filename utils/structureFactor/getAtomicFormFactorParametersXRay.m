function [a, b, c] = getAtomicFormFactorParameters(element)
    %
    %   return the atomic form factors fromhttps://lampz.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    %  [a1 b1 a2 b2 a3 b3 a4 b4 c]
    %   f(g) = sum(i=1:4) a_i * exp( -b_i *(G/4pi)^2) +c)
    %
    %
    %
    %
    allFactors = readtable('xRayFormFactors.csv','ReadRowNames',1);
    %%
    try
        factorParams = allFactors{[element],:};
        %[a1 b1 a2 b2 a3 b3 a4 b4 c]
        a = [factorParams(:,[1,3,5,7])];
        b = [factorParams(:,[2,4,6,8])];
        c = [factorParams(:,9)];
    catch
        warning('Element not specified, or mistyped')
        a = [1 1 1 1];
        b = [1 1 1 1];
        c = [1];
    end
    %%
end