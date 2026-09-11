function corners=get_lattice_corners(h,varargin)
    n_cells=1;
    extend_neg = 0;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
        elseif strcmp(varargin{i},'extend_neg')
            extend_neg = 1;
        else
            continue;
        end
    end
    %%
    [x_ax,y_ax,z_ax]=get_base_vecs(h);
    x_ax=x_ax*n_cells;
    y_ax=y_ax*n_cells;
    z_ax=z_ax*n_cells;
    if extend_neg == 0
        corners= [0 0 0; x_ax; y_ax; z_ax;
        x_ax+y_ax; x_ax+z_ax; y_ax+z_ax;
        x_ax+y_ax+z_ax;]';
    else
        corners= [-(x_ax+y_ax+z_ax); -(x_ax+y_ax)+z_ax; -(x_ax+z_ax)+y_ax; -(y_ax+z_ax)+x_ax;
        x_ax+y_ax-z_ax; x_ax+z_ax-y_ax; y_ax+z_ax-x_ax;
        x_ax+y_ax+z_ax;]';
    end
end

