
function I = mutualInfo(A, B)
    % MUTUALINFO Calculates the mutual information between two images.
    %
    %   I = MUTUALINFO(A, B)        Mutual information of images A and B of
    %   identical size, using 256 bins for histograms.
    %
    %   It is assumed that 0*log2(0) = 0.
    %
    %   The result of MUTUALINFO(X, X) is nearly equal to that of
    %   entropy(X) for an arbitrary image X. Theoretically, they should be
    %   exactly equal, but are not in practice due to floating-point error.
    %   
    %   See also ENTROPY.
    %
    %   Adapted from <a href="matlab:web('https://www.mathworks.com/matlabcentral/fileexchange/13289-fast-mutual-information-of-two-images-or-signals')">MI</a> by Jose Delpiano.
    %
    %   References:
    %   T. M. Cover and J. A. Thomas, "Entropy, Relative Entropy, and
    %       Mutual Information," in Elements of Information Theory, 2nd ed.
    %       Hoboken, NJ: Wiley-Interscience, 2006, pp. 13–55.
    %   F. Maes, D. Loeckx, D. Vandermeulen, and P. Suetens, "Image
    %       registration using mutual information," in Handbook of
    %       Biomedical Imaging: Methodologies and Clinical Research,
    %       N. Paragios, J. Duncan, and N. Ayache, Eds., Boston, MA:
    %       Springer US, 2015, pp. 295–308. doi: 10.1007/978-0-387-09749-7_
    %       16.
    % Handle variable type same way as entropy()
    if ~islogical(A)
        A = im2uint8(A);
        B = im2uint8(B);
        L = 256;
    else
        L = 2;
    end
    
    % Calculate gray-level probability distributions
    pa = histcounts(A, L, ...
        BinLimits=[0 L-1], ...
        BinMethod='integers', ...
        Normalization='probability');
    pb = histcounts(B, L, ...
        BinLimits=[0 L-1], ...
        BinMethod='integers', ...
        Normalization='probability');
    pab = histcounts2(A, B, L, ...
        XBinLimits=[0 L-1], ...
        YBinLimits=[0 L-1], ...
        BinMethod='integers', ...
        Normalization='probability');
    % Calculate mutual information
    papb = pa'*pb;
    idx = papb > 1e-12 & pab > 1e-12;   % Screen out singularities
    I = sum(pab(idx).*log2(pab(idx) ./ papb(idx)));
end
