function plot_lattice_plane(outputs,h,dists,varargin)
    plot_edges=0;
    n_cells=1;
    plot_planes=true;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'plot_edges')
            plot_edges=1;
        elseif strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
        elseif strcmp(varargin{i},'plot_plane')
            plot_planes=varargin{i+1};
        else
            continue;
        end
    end
    outputs=rep_atoms(outputs,n_cells,h,varargin{:});
    all_atoms=outputs.all_atoms;
    cs=h.CS;
    colors=outputs.colors_elements;
%     figure
    scatter3(all_atoms.x,all_atoms.y,all_atoms.z,200,colors(:,1:3),'filled','MarkerEdgeColor','black');
    
    
%     hold on
%     quiver3(0,0,0,cs.aAxis.x,cs.aAxis.y,cs.aAxis.z,'black','linewidth',2)
%     quiver3(0,0,0,cs.bAxis.x,cs.bAxis.y,cs.bAxis.z,'black','linewidth',2)
%     quiver3(0,0,0,cs.cAxis.x,cs.cAxis.y,cs.cAxis.z,'black','linewidth',2)
    axis equal
    hold on
    %%
    %%
    if plot_planes
        color=get_colors_from_list(size(h,2));
        for i=1:size(h,2)
%             disp(i)
            hold on
            col=color{i};
            if size(dists,2)==1
                intersections=get_plane_intersections(h(i),dists,varargin{:});
                p=fill3(intersections(1,:),intersections(2,:),intersections(3,:),col(1:3));
                p.FaceAlpha = 0.5;
            else
                intersections=get_plane_intersections(h(i),dists(i),varargin{:});
                p=fill3(intersections(1,:),intersections(2,:),intersections(3,:),col(1:3));
                p.FaceAlpha = 0.5;
            end
        end
        
    end
    if plot_edges==1
            edges=get_lattice_edges(h(1),varargin{:});
            for i=1:size(edges,3)
                plot3(edges(1,:,i),edges(2,:,i),edges(3,:,i),'black');
            end
            [x_ax,y_ax,z_ax]=get_base_vecs(h);
            plot3([0 x_ax(1)], [0 x_ax(2)], [0 x_ax(3)],'black', 'linewidth', 2)
            plot3([0 y_ax(1)], [0 y_ax(2)], [0 y_ax(3)], 'black', 'linewidth', 2)
            plot3([0 z_ax(1)], [0 z_ax(2)], [0 z_ax(3)], 'black', 'linewidth', 2)
    end
end
