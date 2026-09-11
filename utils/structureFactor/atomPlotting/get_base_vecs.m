function [x_ax,y_ax,z_ax]=get_base_vecs(h)
    base_vecs=h.CS.axes;
    % base_vecs;
    x_ax=[base_vecs(1).x base_vecs(1).y base_vecs(1).z];
    y_ax=[base_vecs(2).x base_vecs(2).y base_vecs(2).z];
    z_ax=[base_vecs(3).x base_vecs(3).y base_vecs(3).z];

end