function plotPattern(EBSD_geometry,Pattern1,currentAxis)
    imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,Pattern1)
    axis equal;
    axis xy;
    colormap gray;
    ylim([EBSD_geometry.y_gn_min,EBSD_geometry.y_gn_max])
    clim([-3.3,3.3])
    % xlim([EBSD_geometry.x_gn_min,EBSD_geometry.x_gn_max])
    axis off
    colormap(currentAxis,'gray')
    % colorbar()
end