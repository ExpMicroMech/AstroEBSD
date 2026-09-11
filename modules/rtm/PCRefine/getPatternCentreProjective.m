function [PCX,PCY,PCZ] = getPatternCentreProjective(x,y,P)
    %%


    % check dimensionallity is correct
    if size(x,2)==1 & size(x,1)>size(x,2)
        x = x';
        y = y';
    end
    
%%
    PCs = P(1:3,:)*[x;y;ones(size(x))];
    %
    scale = P(4,1)*x + P(4,2)*y +1;
    PCs = PCs./scale;
    PCX = PCs(1,:);
    PCY = PCs(2,:);   
    PCZ = PCs(3,:);
end