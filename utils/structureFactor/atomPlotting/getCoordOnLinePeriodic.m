function distance =  getCoordOnLine(origin,direction,queryPoint)
      direction = normalize(direction);

      corners = get_lattice_corners(direction,'n_cells',1);
      bounds=alphaShape(corners(1,:)',corners(2,:)',corners(3,:)',Inf);
    distanceQuery = -1/max(corners(:)):0.05:1.8*max(corners(:));
      distance = minus(dot(queryPoint,direction)' , dot(origin,direction))
      %%
end