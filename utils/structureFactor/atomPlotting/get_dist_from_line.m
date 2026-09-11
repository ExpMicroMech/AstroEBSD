function distance=get_dist_from_line(normal,origin,points)
    normal=normal./norm(normal);
    originNew = vector3d(origin.x,origin.y,origin.z);
    pointsNew = vector3d(points.x, points.y, points.z);
    % i=4
    % dist = 
    distance =  norm((originNew - pointsNew') - (dot((originNew-pointsNew'),normal) * normal) );
    % figure
    % imagesc(distance)
    % colorbar
    %%
    %%https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane
    %%
    
end
