function plotPCRefinementResults(PCX,PCY,PCZ,PCXGrid_unfiltered,PCYGrid_unfiltered,PCZGrid_unfiltered,x_unfiltered,nDelta)
    figure
    scatter3(PCX,PCY,PCZ)
    hold on
    % scatter3(PCXGrid_unfiltered,PCYGrid_unfiltered,PCZGrid_unfiltered)
    PCXGrid_r = reshape(PCXGrid_unfiltered,size(x_unfiltered));
    PCYGrid_r = reshape(PCYGrid_unfiltered,size(x_unfiltered));
    PCZGrid_r = reshape(PCZGrid_unfiltered,size(x_unfiltered));
    mesh(PCXGrid_r,PCYGrid_r,PCZGrid_r);
    % figure
    quiver3(PCX,PCY,PCZ,-nDelta(1,:),-nDelta(2,:),-nDelta(3,:),0);
    view(-45,10)
    axis equal
    % exportgraphics(gcf,fullfile(HDF5Folder,sprintf('PCPlane%sFit.png',(PatternData(RefinementNumber).DataName))))
end
% end