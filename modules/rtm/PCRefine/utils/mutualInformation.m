
function MutInf= mutualInformation(patSim,patExp)
    %% C. Zhu et al Ultramicroscopy 233 (2022) 113407
    % % patSim = imresize(patSim,2);
    % patExp =  imresize(patExp,2);
   pA =pixelprobability(patSim);
   pB =pixelprobability(patExp);
   pJoint = jointPixProbability(patSim,patExp);
   % pJoint = pixelprobability(pA.*pB')
   MutInf = -sum(pA(:).*log2(pA(:)),'omitnan') -sum(pB(:).*log2(pB(:)),'omitnan') + sum(pJoint(:).*log2(pJoint(:)),'omitnan')

end
% optimisationLandscapeSynthetic
function p = pixelprobability(img)
    img = uint8(rescale(img,0,255));
    h = histcounts(img,256);
    p = h/numel(img(:));
end

function jointProb1DNoZero = jointPixProbability(im1,im2)
    % im1 =doulbe
    im1 = uint8(rescale(im1,0,255));
    im2 = uint8(rescale(im2,0,255));
    indrow = double(im1(:)) + 1;
    indcol = double(im2(:)) + 1; %// Should be the same size as indrow
    jointHistogram = accumarray([indrow indcol], 1);
    jointProb = jointHistogram / numel(indrow);
    indNoZero = jointHistogram ~= 0;
    jointProb1DNoZero = jointProb(indNoZero);
    % jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));
end