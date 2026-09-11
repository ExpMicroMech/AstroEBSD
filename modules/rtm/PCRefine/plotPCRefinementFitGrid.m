function plotPCRefinementFitGrid(PCXGrid,PCYGrid,PCZGrid,PCX,PCY,PCZ)
%
%plotPCRefinementFitGrid(PCXGrid,PCYGrid,PCZGrid,PCX,PCY,PCZ)
% - Input Parameters :
% PCXGrid, PCYGrid, PCZGrid -> Pattern Centre coordinates after fitting a
% "plane"
% (using either projective or affine method
% PCX, PCY, PCZ -> Projection center coordinates from e.g. pattern matching
% used in the fit to obtain PCXGrid .. initially.
    f=figure
    f.Position = [100 100 1200 400]
    subplot(1,3,1)
    scatter(PCXGrid(:),PCYGrid(:),30,[0.5 0.5 0.5], 'filled','DisplayName','PC, Plane fitted')
    hold on
    scatter(PCX,PCY,10,[0 1 1],'filled','DisplayName','PC, input ','MarkerEdgeColor',[0.3 0.3 0.3])
    xlabel('PCX')
    ylabel('PCY')
    axis equal
    legend()
    
    % f=figure
    % f.Position = [100 100 1200 400]
    subplot(1,3,2)
    scatter(PCXGrid(:),PCZGrid(:),30,[0.5 0.5 0.5], 'filled','DisplayName','PC, Plane fitted')
    hold on
    scatter(PCX,PCZ,10,[0 1 1],'filled','DisplayName','PC, input ','MarkerEdgeColor',[0.3 0.3 0.3])
    xlabel('PCX')
    ylabel('PCZ')
    axis equal
    % f=figure
    % f.Position = [100 100 1200 400]
    subplot(1,3,3)
    scatter(PCZGrid(:),PCYGrid(:),30,[0.5 0.5 0.5], 'filled','DisplayName','PC, Plane fitted')
    hold on
    scatter(PCZ,PCY,10,[0 1 1],'filled','DisplayName','PC, input ','MarkerEdgeColor',[0.3 0.3 0.3])
    xlabel('PCZ')
    ylabel('PCY')
    % end
    
    % close all
end