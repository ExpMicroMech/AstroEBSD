function dist=get_dist_from_plane(normal,origin,points)
    normal=normal./norm(normal);
    a=normal(1);
    b=normal(2);
    c=normal(3);

    x=points.x;
    y=points.y;
    z=points.z;

    px=origin(1);
    py=origin(2);
    pz=origin(3);
    D = sqrt(px^2+py^2+pz^2);
    d = -(a*px+b*py+c*pz);
    dist= (a.*(x-px)+b.*(y-py)+c.*(z-pz))
    %%
    %%https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane
    %%
    
end
