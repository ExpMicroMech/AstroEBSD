function [atoms_return,layers]=get_atoms_layerwise(atoms_in,h,varargin)
%     
    n_cells=1;

%     h=Miller(1,0,-1,10,cs);
%     corners=get_lattice_corners(h,'n_cells',n_cells)
    n=[h.x;h.y;h.z];  %% we use column vectors here
%     n=[-1;1;10]
    n=n./norm(n);
    [x_ax,y_ax,z_ax]=get_base_vecs(h);
    l=[x_ax;y_ax;z_ax]';
    %
%     nLen_non=[x_ax*n, y_ax*n,z_ax*n]
%     norm_n=norm(n)
%         temp=norm([norm(x_ax)*n,norm(y_ax)*n,norm(z_ax)*n])
%     n=n./norm(n);
    diag=x_ax+y_ax+z_ax;   
    nLen=norm([diag*n]);%%gives the lattice paramters if e.g. n=[0;0;1], so it is the "extend" of the Miller in the unit cell
    %%
    slices=(nLen/h.dspacing)*n_cells; % meaning one unit cell
    circular=0;
    dSpacing = h.dspacing
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'slices')
            slices=varargin{i+1};
        elseif strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
            slices=(nLen/h.dspacing)*n_cells;
            if ~mod(slices,1)==0
                slices = ceil(slices);
                dSpacing = nLen/2;
            end
        elseif strcmp(varargin{i},'circular')
            circular=1;
        else
            continue;
        end
        
    end
    %%
%     figure 
%     hold on
%     axis equal;
%     
%%
    
    if circular==1
        atoms=rep_atoms(atoms_in,ceil([max(norm(h.CS.axes))./norm(h.CS.axes)]),h,'extend_negative');
    else
        atoms=rep_atoms(atoms_in,n_cells+1,h,'extend_negative');
    end
%     all_atoms=atoms.all_atoms;
%
    atom_cords=[atoms.all_atoms.x,atoms.all_atoms.y,atoms.all_atoms.z];
%     %%
    color=get_colors_from_list(uint8(slices+1));
%     for i=1:slices
%         disp(i)
%         hold on
    edges=get_lattice_edges(h,varargin{:});
    %%
    corners=get_lattice_corners(h,varargin{:});
%

%         for i=1:size(edges,3)
%             plot3(edges(1,:,i),edges(2,:,i),edges(3,:,i),'black')
%         end
    count=1;
    atoms_return=struct();
    %c
%     figure
%     hold on
% figure
%     for i=0:slices+1
%         
% %         layers(count).intersections_0=edge_ints_0;
% %         layers(count).intersections_1=edge_ints_1;
% %         edge_ints_1=get_plane_intersections(h,nLen*(i+1)/(slices));
%     end

    for i=0:slices-1
%         disp(i)
%%
%         if i<slices
        edge_center=get_plane_intersections(h,(i+0.5)*dSpacing,varargin{:});      
        edge_left=get_plane_intersections(h,(i)*dSpacing,varargin{:}); 
        edge_right=get_plane_intersections(h,(i+1)*dSpacing,varargin{:}); 

        edge_center_plot=get_plane_intersections(h,(i+0.5)*dSpacing);     
        edge_left_plot=get_plane_intersections(h,(i)*dSpacing); 
        edge_right_plot=get_plane_intersections(h,(i+1)*dSpacing); 


        col=color{abs(uint8(i))+1};
%         edge_ints_both=[edge_ints_0 edge_ints_1];
%             bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf);
%%      
        if ~isempty(edge_left)
        if circular==1
%             n_rot=1
% figure
            n_rot=[1:12]/12;
            rand_ort=normalize(orth(h));
            rot=rotation.byAxisAngle(h,n_rot*2*pi);
            outside=rot*rand_ort*norm(norm(h.CS.axes));
            intersections_0=[outside.x outside.y outside.z] +(i)*dSpacing*n';
            intersections_1=[outside.x outside.y outside.z]+ (i+1)*dSpacing*n';
            intersections_both=[intersections_0' intersections_1'];
            intersections_both=intersections_both+(intersections_both-mean(intersections_both,2))*0.03;
            bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf); 
%             plot(bounds)
%%
        else
            intersections_0=get_plane_intersections(h,i*dSpacing,varargin{:});
            intersections_1=get_plane_intersections(h,(i+1)*dSpacing,varargin{:});
            intersections_both=[intersections_0 intersections_1];
            intersections_both=intersections_both+(intersections_both-mean(intersections_both,2))*0.03;
            bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf); 
        end

        try
            %%
            selection=bounds.inShape(atom_cords);
            temp=atom_cords(selection,:);
            m_a=max(get_lattice_corners(h,varargin{:}),[],2);
            m_i=min(get_lattice_corners(h,varargin{:}),[],2);
            if all(temp(:,3)>=m_a(3)) || all(temp(:,3)<m_i(3))
                error('z to high or low')
                disp('z_error')
            end
    %             scatter3(temp(:,1),temp(:,2),temp(:,3))
    %             p=fill3(intersections_0(1,:),intersections_0(2,:),intersections_0(3,:),col(1:3));
    %             p=fill3(intersections_1(1,:),intersections_1(2,:),intersections_1(3,:),col(1:3));
            %
%              if circular==1
%              else

%              end
%             if i<slices
%             if i<=0
%             else
%             end
%             if ~negative_control_var==1
                layers(count).intersections_center=edge_center;
    %             end
                layers(count).intersections_right=edge_left;
                layers(count).intersections_left=edge_right;

                layers(count).intersections_center_plot=edge_center_plot;
    %             end
                layers(count).intersections_right_plot=edge_left_plot;
                layers(count).intersections_left_plot=edge_right_plot;


                atoms_return(count).labels=atoms.labels(selection);
                atoms_return(count).elements=atoms.elements(selection);
                
                atoms_return(count).colors_elements=atoms.colors_elements(selection,:);
                atoms_return(count).colors_label=atoms.colors_label(selection,:);
        %             scatter3(intersections_both(1,:),intersections_both(2,:),intersections_both(3,:))
                atoms_return(count).all_atoms.x=atoms.all_atoms.x(selection);
                atoms_return(count).all_atoms.y=atoms.all_atoms.y(selection);
                atoms_return(count).all_atoms.z=atoms.all_atoms.z(selection);
                atoms_return(count).dSpacing = dSpacing;
                if isfield(atoms, 'radii')
                    atoms_return(count).radii = atoms.radii(selection);
                end
                count=count+1;
%             end
            %%
        catch
            warning('slice probably outside of requested unit cell')
        end
        %%
        end

        %%
        % atoms_return.dSpacing = dSpacing;
    end
    
    %%
% p0=n*dist;
% l0=corners;
% l=[base_vecs.x;base_vecs.y;base_vecs.z];
    
    
end
% function replicate_atoms(atoms,n_unit_cells)
    