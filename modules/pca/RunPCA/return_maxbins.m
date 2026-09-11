function [minind,maxind]=return_maxbins(N)
    [ms,is]=maxk(N,2);
    if ms(2)/ms(1)<0.2  %% if the hist bin is much higher than the rest use this one.
        [~,minind]=max(N);
        maxind=minind+1;
    elseif is(1)-is(2)==1 %% bins are neighboured if this is 1, the 1st index is larger, hence has to be the right boundary
        minind=is(2);
        maxind=is(1)+1;
    elseif is(2)-is(1)==1
        minind=is(1);
        maxind=is(2)+1;
    else
       warning('proably bimodal orientation distribution')
       % return
       [~,minind]=max(N);
       maxind=minind+1; 
    end
    %%
end