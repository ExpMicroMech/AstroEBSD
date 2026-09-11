function edges=get_lattice_edges(h,varargin)
    n_cells=1;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
        else
            continue;
        end
    end
    [x_ax,y_ax,z_ax]=get_base_vecs(h);
    x_ax=x_ax*n_cells;
    y_ax=y_ax*n_cells;
    z_ax=z_ax*n_cells;
    edges=cat(3,[0 0 0; x_ax]',[0 0 0; y_ax]', [x_ax;x_ax+y_ax]',[y_ax;x_ax+y_ax]',... %% "bottom"
        [z_ax; z_ax+x_ax]',[z_ax;z_ax+y_ax]', [z_ax+x_ax;z_ax+x_ax+y_ax]',[z_ax+y_ax;z_ax+x_ax+y_ax]',...%% "top"
         [0 0 0; z_ax]',[y_ax; z_ax+y_ax]', [x_ax;z_ax+x_ax]',[x_ax+y_ax;x_ax+y_ax+z_ax]'); %% sides   
% %         [0 0 0; z_ax];[] y_ax; z_ax;
%                 x_ax+y_ax;x_ax+z_ax;y_ax+z_ax;
%                 x_ax+y_ax+z_ax;]';get_lattice_corners(h)

end
