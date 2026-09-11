function [n,x1] = getSurfaceNormalFromPC(PCXGrid,PCYGrid,PCZGrid)
    n = [0 0 0];
    maxcounts = numel(PCXGrid);
    breakCount = 2
    while sum(n)==0
        if breakCount>= maxcounts
            break
        end
    x0 = [PCXGrid(1), PCYGrid(1),PCZGrid(1)];
    x1 = [PCXGrid(breakCount+1), PCYGrid(breakCount), PCZGrid(2)];
    x2 = [PCXGrid(end), PCYGrid(end), PCZGrid(end)];
    n = cross([x0-x1],[x0-x2]);
    breakCount = breakCount+1;
    end
        %%
    % deltaPC = dot(([PCX',PCY',PCZ']-x1)',repmat(n',size(PCX)))./norm(n)
end