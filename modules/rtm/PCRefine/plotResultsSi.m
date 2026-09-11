%%
function plotResultsSi(HDF5Folder,PatSizes,MaxIters,Phase_Names,DeltaPC_AllfromPoints,DeltaPC_AllfromPointsClean,PatternData,DeltaX_allClean,DeltaY_allClean,DeltaZ_allClean,PCXGrid_allClean,PCYGrid_allClean,PCZGrid_allClean)

Markers = ['^','v',"<",">","s","d","+","*","p","o"]
colors = tab20
uniquePatSizes = unique(PatSizes);
uniuqePhase_Names = unique(Phase_Names)
uniqueIters = unique(MaxIters);
LineWidth =2
MarkerLineWidth =1;
MarkerSize = 40
    % for maxIter in 
f=figure
f.Position = [100 100 1200 1300]
hold on
h = tiledlayout(3,3, 'TileSpacing', 'tight', 'Padding', 'none');
colors = tab20
% bad variable naming but this makes switching the
% loop easier if we rather want to plot the number of iterations
phaseNumber = 1
% for patSize =uniquePatSizes
%%
for iters = uniqueIters
    
    % patSizeIndex = PatSizes == patSize;
   phaseIndex = iters ==MaxIters;
    nexttile(1)
    hold on
    l=plot(PatSizes(phaseIndex)', mean(DeltaPC_AllfromPoints(phaseIndex,:), 2), 'handleVisibility', 'off', 'LineWidth', LineWidth);
    
    
    s=scatter(PatSizes(phaseIndex), mean(DeltaPC_AllfromPoints(phaseIndex,:), 2), MarkerSize,...
        Markers(phaseNumber),'filled', 'DisplayName', string(iters),'MarkerEdgeColor','black');
    l.Color = colors(phaseNumber,:)
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:)
   

   %%
    
    nexttile(2)
    l=plot(PatSizes(phaseIndex), mean(DeltaPC_AllfromPointsClean(phaseIndex,:), 2), 'DisplayName', string(iters), ...
        'LineWidth', LineWidth);
     l.Color = colors(phaseNumber,:)
    
    s=scatter(PatSizes(phaseIndex), mean(DeltaPC_AllfromPointsClean(phaseIndex,:), 2),...
        MarkerSize,Markers(phaseNumber),'filled','DisplayName', string(iters),'MarkerEdgeColor','black' );
    hold on
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:)



    % s.Filled =1;
    nexttile(3)
    l=plot(PatSizes(phaseIndex), mean(cat(1, PatternData(phaseIndex).NCC), 2),'DisplayName', string(iters), ...
        'LineWidth', LineWidth);
     l.Color = colors(phaseNumber,:)

    s=scatter(PatSizes(phaseIndex), mean(cat(1, PatternData(phaseIndex).NCC), 2), ...
        MarkerSize,Markers(phaseNumber),'filled', 'DisplayName', string(iters), ...
        'Marker', Markers(phaseNumber),'MarkerEdgeColor','black' );
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:)
    % l.EdgeColor = colors(phaseNumber,:)
    hold on



    nexttile(4)
    l=plot(PatSizes(phaseIndex), DeltaX_allClean(phaseIndex),'handleVisibility','off','LineWidth', LineWidth);
    l.Color = colors(phaseNumber,:);

    s=scatter(PatSizes(phaseIndex), DeltaX_allClean(phaseIndex), MarkerSize, Markers(phaseNumber),'filled',...
        'DisplayName', string(iters),'MarkerEdgeColor','black');
    s.MarkerFaceColor = colors(phaseNumber,:);
    hold on
    
    nexttile(5)
    l=plot(PatSizes(phaseIndex), DeltaY_allClean(phaseIndex),'handleVisibility','off','LineWidth', LineWidth);
    l.Color = colors(phaseNumber,:);

    s=scatter(PatSizes(phaseIndex), DeltaY_allClean(phaseIndex),MarkerSize,...
        Markers(phaseNumber),'DisplayName', string(iters),'MarkerEdgeColor','black' );
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:);
    hold on
    
    nexttile (6 )
    l=plot ( PatSizes ( phaseIndex ) , DeltaZ_allClean ( phaseIndex ) ,  'DisplayName' ,string(iters)  ,...
        'LineWidth' , LineWidth );
     l.Color = colors(phaseNumber,:);
    
    s=scatter ( PatSizes ( phaseIndex ) , DeltaZ_allClean ( phaseIndex ), ...
        MarkerSize ,Markers(phaseNumber),'filled', 'DisplayName' ,string(iters) ,'MarkerEdgeColor','black');
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:);




    
    % l.Color = colors(phaseNumber,:)
    hold on
%%
    meanPCX = [];
    meanPCY = [];
    meanPCZ = [];
    selectedPhases = find(phaseIndex==1)
    for numSelPhases = 1:numel(selectedPhases)
    findZeroPoints = PCXGrid_allClean(selectedPhases(numSelPhases),:)==0;
    % find(phaseIndex==1)
    %%
    meanPCX = [meanPCX mean(PCXGrid_allClean(selectedPhases(numSelPhases),~findZeroPoints))];
    meanPCY = [meanPCY mean(PCYGrid_allClean(selectedPhases(numSelPhases),~findZeroPoints))];
    meanPCZ = [meanPCZ mean(PCZGrid_allClean(selectedPhases(numSelPhases),~findZeroPoints))];
    % temp=pcxPlot(findZeroPoints)
    end  
    meanPCX = meanPCX;
    meanPCY = meanPCY;
    meanPCZ = meanPCZ;
    %
     nexttile(7)
    l=plot(PatSizes(phaseIndex),meanPCX,'handleVisibility','off','LineWidth', LineWidth);
    l.Color = colors(phaseNumber,:)

    s=scatter(PatSizes(phaseIndex), meanPCX, MarkerSize, Markers(phaseNumber),'filled',...
        'DisplayName', string(iters),'MarkerEdgeColor','black');
    s.MarkerFaceColor = colors(phaseNumber,:);
    hold on
    
    nexttile(8)
    l=plot(PatSizes(phaseIndex), meanPCY,'handleVisibility','off','LineWidth', LineWidth);
    l.Color = colors(phaseNumber,:);

    s=scatter(PatSizes(phaseIndex),meanPCY,MarkerSize,...
        Markers(phaseNumber),'DisplayName', string(iters),'MarkerEdgeColor','black' );
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:);
    hold on
    
    nexttile (9 )
    l=plot ( PatSizes ( phaseIndex ) ,meanPCZ ,  'DisplayName' ,string(iters)  ,...
        'LineWidth' , LineWidth );
     l.Color = colors(phaseNumber,:);
    
    s=scatter ( PatSizes ( phaseIndex ) ,meanPCZ, ...
        MarkerSize ,Markers(phaseNumber),'filled', 'DisplayName' , string(iters),'MarkerEdgeColor','black');
    s.LineWidth = MarkerLineWidth;
    s.MarkerFaceColor = colors(phaseNumber,:);
    % l.Color = colors(phaseNumber,:)
    hold on


    phaseNumber = phaseNumber+1

end
nexttile(1)
ylabel('Mean Deviation from Plane')

title('Plane deviation before cleaning')
nexttile(2)
title('Plane deviation after cleaning')
nexttile(3)
title('mean normalised cross correlation')
nexttile(4)
%

ylabel('rot Angle, °')
title('angle x-axis')
nexttile(5)

title('angle y-axis')
nexttile(6)

title('angle z-axis')
nexttile(7)
title('mean PCx')

nexttile(8)
title('mean PCy')
legend('location','southoutside','NumColumns',3)

nexttile(9)
title('mean PCz')


% axis off
exportgraphics(gcf,fullfile(HDF5Folder,'PatternFittingResulsByPhase.png'))
end