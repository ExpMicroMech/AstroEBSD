function at_out=rep_atoms(at,n_cells,h,varargin)
    
    extend_negative=0;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'extend_negative')
            extend_negative=1;
        else
            continue;
        end
        
    end
    
    if isscalar(n_cells)
        n_cells_x=n_cells;
        n_cells_y=n_cells;
        n_cells_z=n_cells;
    else
        n_cells_x=n_cells(1);
        n_cells_y=n_cells(2);
        n_cells_z=n_cells(3);
    end

    if extend_negative==1
        start_cell_x=-n_cells_x;
        start_cell_y=-n_cells_y;
        start_cell_z=-n_cells_z;
    else
        start_cell_x=0;
        start_cell_y=0;
        start_cell_z=0;
    end

    elements=[];
    labels=[];
    colors_label=[];
    colors_elements=[];
    all_atoms=[];
    all_rep=[];
    radii = [];

    [x_ax,y_ax,z_ax]=get_base_vecs(h);
      [x_ax,y_ax,z_ax]=get_base_vecs(h);
%     mults=zeros(n_cells_x,n_cells_y,n_cells_z,3);
    at_cord=[at.all_atoms.x at.all_atoms.y at.all_atoms.z];
    % at_miller = [at]
    for i=start_cell_x:n_cells_x-1
        for j=start_cell_y:n_cells_y-1
            for k=start_cell_z:n_cells_z-1
    %             mults(i+1,j+1,k+1,:)=[i*x_ax + j*y_ax + k*z_ax]
                all_rep=vertcat(all_rep,at_cord+[i*x_ax + j*y_ax + k*z_ax]);
                labels=vertcat(labels, at.labels);
                elements=vertcat(elements, at.elements);
                colors_label=vertcat(colors_label, at.colors_label);
                colors_elements=vertcat(colors_elements, at.colors_elements);
                if isfield(at,'radii')
                    radii = vertcat(radii, at.radii);
                end
%                 all_atoms=[all_atoms at.all_atoms];

            end
        end
    end
%%
%     all_rep=Miller(all_rep(:,1),all_rep(:,2),all_rep(:,3),'UVW',h.CS);
    at_out.all_atoms.x=all_rep(:,1);
    at_out.all_atoms.y=all_rep(:,2);
    at_out.all_atoms.z=all_rep(:,3);
    at_out.labels=labels;
    at_out.elements=elements;
    at_out.colors_label=colors_label;
    at_out.colors_elements=colors_elements;
    if exist('radii')
        at_out.radii = radii;
    end
    
    %%
    
end

