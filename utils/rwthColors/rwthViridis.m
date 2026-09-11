function cmap = rwthViridis(numValues)
    cmap(1,:) = rwthViolett(100)
    cmap(2,:) = rwthTuerkis(100)
    % cmap(3,:) = rwthTuerkis(100)%// color first row - red
    cmap(3,:) = rwthGruen(100)   %// color 25th row - green
    cmap(4,:) = rwthGelb(100);   %// color 50th row - blue
    
    [X,Y] = meshgrid([1:3],[1:256]);  %// mesh of indices
    %
    pos1 = 130/256
    pos2 = 175/256
    % if 
       % [X,Y] = meshgrid([1:3],[1:256]);  %// mesh of indices
    %%
    cmap = interp2(X([1,130,175,256],:),Y([1,130,175,256],:),cmap,X,Y);%// interpolate colormap   % colormap(cmap)
    cmap = cmap(uint16(linspace(1,256,numValues)),:)
    %%
end