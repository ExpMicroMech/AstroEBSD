function ebsd_pat = kernelAveraging(x,y,kernelSize,MapInfo)
    %%
    % TODO: Implement kernel
    if mod(kernelSize,2)==0
        kernelSize = kernelSize+1;
    end
    % kernel = [0.5 0.5 0.5; 0.5 1 0.5; 0.5 0.5 0.5]
    [x_cord,y_cord] = (meshgrid(1:kernelSize,1:kernelSize));
    x_cord =uint32(x_cord) +x- ceil(kernelSize/2);
    % x_cord;
    y_cord = uint32(y_cord) +y -ceil(kernelSize/2);
    [pnum,~]=find(MapInfo.MapData.XBeam==x_cord(:)' & MapInfo.MapData.YBeam==y_cord(:)');
    %% for more efficient pattern reading
    pnum_sort = sort(pnum');
    jumpList = pnum_sort(2:end)-pnum_sort(1:end-1);
    jumps = find(jumpList>1);
    jumps = [0 jumps numel(pnum_sort)];
    % patterns = zeros(size())
    tempPatList = [];
    for i = 1:size(jumps,2)-1
        index = [pnum_sort(jumps(i)+1),jumps(i+1)-(jumps(i))  ];
        % pnum_sort(index)
        tempPat = bReadEBSP(MapInfo.EBSPData,index);%,Settings_Cor;
        tempPatList = cat(3,tempPatList,tempPat);
    end
    %%
    % for i = 1:size(pnum,1)
    %     []
    % 
    % end

    ebsd_pat = mean(tempPatList,3);
end