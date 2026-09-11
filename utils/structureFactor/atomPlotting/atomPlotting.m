

%%
% x_ax=[base_vecs(1).x base_vecs(1).y base_vecs(1).z];
% y_ax=[base_vecs(2).x base_vecs(2).y base_vecs(2).z];
% z_ax=[base_vecs(3).x base_vecs(3).y base_vecs(3).z];
% corners= [0 0 0; x_ax; y_ax; z_ax;
%             x_ax+y_ax;x_ax+z_ax;y_ax+z_ax;
%             x_ax+y_ax+z_ax;]';
% n_cells=1;
% for i=1:size(varargin(:))
%         if strcmp(varargin{i},'n_cells')
%             n_cells=varargin{i+1};
%         else
%             continue;
%         end
%     end
[x_ax,y_ax,z_ax]=get_base_vecs(h);
corners=get_lattice_corners(h,varargin{:});


dot_direction=round(dot(repmat(h.xyz',[1 size(corners,2)]),corners),3);
    %%
index_min=find(min(dot_direction)==dot_direction);
index_min=index_min(1);

% [~,index_min]=min(dot(repmat(h.xyz',[1 size(corners,2)]),corners));
start_corner=corners(:,index_min);

n=[h.x;h.y;h.z];  %% we use column vectors here
n_non=n;
n=n./norm(n);
p0=n*dist+start_corner;
% [~,ind]=min(n);
if ~ n(3)==0
    ind=3;
else
    inf_e=0;
end
l0=corners;
l=[x_ax;y_ax;z_ax]';

% bounds=alphaShape(corners(1,:)',corners(2,:)',corners(3,:)',Inf);
% bounds.Alpha=bounds.Alpha*1.2;
corners_bound=corners+(corners-mean(corners,2))*0.03;

% corners
bounds=alphaShape(corners_bound(1,:)',corners_bound(2,:)',corners_bound(3,:)',Inf);
if extend~=1
    
end
%%

if inf_e==1 & extend~=1
    corners_bound=corners+(corners-mean(corners,2))*extend;
    corners_bound(ind,:)=corners_bound(ind,:)+(corners_bound(ind,:)-mean(corners_bound(ind,:),2))*5;
    bounds=alphaShape(corners_bound(1,:)',corners_bound(2,:)',corners_bound(3,:)',Inf);
elseif inf_e==1
    corners_bound(ind,:)=corners_bound(ind,:)+(corners_bound(ind,:)-mean(corners_bound(ind,:),2))*5;
    bounds=alphaShape(corners_bound(1,:)',corners_bound(2,:)',corners_bound(3,:)',Inf);
end
%
% figure
% hold on
intersections=[];
for i=1:size(l0,2)
    for j=1:size(l)
        %%
%         to_plot=[l0(:,i) l0(:,i)+l(:,j)]
%         plot(to_plot)
%         plot3(to_plot(1,:),to_plot(2,:),to_plot(3,:))
       if dot(p0-l0(:,i),n)==0
            d=0;
            intersections=[intersections l0(:,i)];
       elseif ~dot(n,l(:,j))==0
            d=dot(p0-l0(:,i),n)/dot(n,l(:,j));
            intersections=[intersections l0(:,i)+l(:,j)*d];
%             disp(d)
        
        end
    end
end
%%
% figure
% scattrer3(intersections(1,:),intersections(2,:),intersections(3,:),'black');
% pts=bounds.Points;
% % f
% hold on
% % plot(bounds)
% hold on
% % plot3([0 n(1)],[0 n(2)], [0 n(3)]);
% % scatter3(pts(1,:),pts(2,:),pts(3,:),10)
% axis equal;
%%
    intersections=intersections(:,bounds.inShape(intersections'));
    %%
    intersections=unique(round(intersections,4)','rows')'; % this needs to be transposed
    center=mean(intersections,2);
    std_var=std(intersections,0,2);
    rs=center-intersections;
%     intersections;
% scatter3(intersections(:,1),intersections(:,2),intersections(:,3))
    %% intersections
    [min_min,ind_min]=min(std_var); %% e.g. if we have vertical problems, the ordering  does not work as well.
   %%
    if ind_min ==1
        angle=atan2(rs(2,:),rs(3,:));
        [~,ind]=sort(angle);
        intersections=intersections(:,ind);
    elseif ind_min==2
        angle=atan2(rs(1,:),rs(3,:));
        [~,ind]=sort(angle);
        intersections=intersections(:,ind);
    elseif ind_min==3
        angle=atan2(rs(1,:),rs(2,:));
        [~,ind]=sort(angle);
        intersections=intersections(:,ind);
    end
    %%
end
% end

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
    scatter3(all_atoms.x,all_atoms.y,all_atoms.z,100,colors(:,1:3),'filled','MarkerEdgeColor','black');
    
    
%     hold on
%     quiver3(0,0,0,cs.aAxis.x,cs.aAxis.y,cs.aAxis.z,'black','linewidth',2)
%     quiver3(0,0,0,cs.bAxis.x,cs.bAxis.y,cs.bAxis.z,'black','linewidth',2)
%     quiver3(0,0,0,cs.cAxis.x,cs.cAxis.y,cs.cAxis.z,'black','linewidth',2)
    axis equal
    hold on
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

function corners=get_lattice_corners(h,varargin)
    n_cells=1;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
        else
            continue;
        end
    end
    %%
    [x_ax,y_ax,z_ax]=get_base_vecs(h);
    x_ax=x_ax*n_cells;
    y_ax=y_ax*n_cells;
    z_ax=z_ax*n_cells;
    corners= [0 0 0; x_ax; y_ax; z_ax;
                x_ax+y_ax;x_ax+z_ax;y_ax+z_ax;
                x_ax+y_ax+z_ax;]';
end
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


function [x_u,y_u,z_u]=restrict_to_unit_cell(x,y,z)
    x_u=x-floor(x)
    y_u=y-floor(y)
    z_u=z-floor(z)
end
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
    slices=nLen/h.dspacing; % meaning one unit cell
    circular=0;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'slices')
            slices=varargin{i+1};
        elseif strcmp(varargin{i},'n_cells')
            n_cells=varargin{i+1};
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
    edges=get_lattice_edges(h);
    
    corners=get_lattice_corners(h);
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
        edge_center=get_plane_intersections(h,(i+0.5)*h.dspacing,varargin{:});      
        edge_left=get_plane_intersections(h,(i)*h.dspacing,varargin{:}); 
        edge_right=get_plane_intersections(h,(i+1)*h.dspacing,varargin{:}); 

        edge_center_plot=get_plane_intersections(h,(i+0.5)*h.dspacing);     
        edge_left_plot=get_plane_intersections(h,(i)*h.dspacing); 
        edge_right_plot=get_plane_intersections(h,(i+1)*h.dspacing); 


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
            intersections_0=[outside.x outside.y outside.z] +(i)*h.dspacing*n';
            intersections_1=[outside.x outside.y outside.z]+ (i+1)*h.dspacing*n';
            intersections_both=[intersections_0' intersections_1'];
            intersections_both=intersections_both+(intersections_both-mean(intersections_both,2))*0.03;
            bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf); 
%             plot(bounds)
%%
        else
            intersections_0=get_plane_intersections(h,i*h.dspacing,varargin{:});
            intersections_1=get_plane_intersections(h,(i+1)*h.dspacing,varargin{:});
            intersections_both=[intersections_0 intersections_1];
            intersections_both=intersections_both+(intersections_both-mean(intersections_both,2))*0.03;
            bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf); 
        end

        try
            %%
            selection=bounds.inShape(atom_cords);
            temp=atom_cords(selection,:);
            m_a=max(get_lattice_corners(h),[],2);
            m_i=min(get_lattice_corners(h),[],2);
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
                count=count+1;
%             end
            %%
        catch
            warning('slice probably outside of requested unit cell')
        end
        %%
        end

        %%
    end
    %%
% p0=n*dist;
% l0=corners;
% l=[base_vecs.x;base_vecs.y;base_vecs.z];
    
    
end
% function replicate_atoms(atoms,n_unit_cells)
    
function colors=get_colors_from_list(n)
    color_list= {[0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0];...
%  [0.6823529411764706, 0.7803921568627451, 0.9098039215686274, 1.0];...
 [1.0, 0.4980392156862745, 0.054901960784313725, 1.0];...
%  [1.0, 0.7333333333333333, 0.47058823529411764, 1.0];...
 [0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0];...
%  [0.596078431372549, 0.8745098039215686, 0.5411764705882353, 1.0];...
 [0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0];...
%  [1.0, 0.596078431372549, 0.5882352941176471, 1.0];...
 [0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0];...
%  [0.7725490196078432, 0.6901960784313725, 0.8352941176470589, 1.0];...
 [0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0];...
%  [0.7686274509803922, 0.611764705882353, 0.5803921568627451, 1.0];...
 [0.8901960784313725, 0.4666666666666667, 0.7607843137254902, 1.0];...
 [0.9686274509803922, 0.7137254901960784, 0.8235294117647058, 1.0];...
 [0.4980392156862745, 0.4980392156862745, 0.4980392156862745, 1.0];...
 [0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0];...
 [0.7372549019607844, 0.7411764705882353, 0.13333333333333333, 1.0];...
 [0.8588235294117647, 0.8588235294117647, 0.5529411764705883, 1.0];...
 [0.09019607843137255, 0.7450980392156863, 0.8117647058823529, 1.0];...
 [0.6196078431372549, 0.8549019607843137, 0.8980392156862745, 1.0]};
    colors=color_list(1:n,:);
end

% function get_atoms_within_bound()
% end
function [x_ax,y_ax,z_ax]=get_base_vecs(h)
    base_vecs=h.CS.axes;
    % base_vecs;
    x_ax=[base_vecs(1).x base_vecs(1).y base_vecs(1).z];
    y_ax=[base_vecs(2).x base_vecs(2).y base_vecs(2).z];
    z_ax=[base_vecs(3).x base_vecs(3).y base_vecs(3).z];

end
%%
% function rot_mat=get_h_rot(h)
%     rot_x=rotation.byAxisAngle(vector3d.X,angle(vector3d.Z,h(1)));
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.Y,h(1)));
%     rot_mat=rot_x*rot_z;
% end
function rot_mat=get_h_rot(h)
% nach h.j. bunge 1988
    Phi= acos(h(1).l/sqrt(h(1).h^2+h(1).k^2+h(1).l^2));
    %%;
    phi2=acos(h(1).k/sqrt(h(1).h^2+h(1).k^2));
    %%
%     rot_x=rotation.byAxisAngle(vector3d.X,angle(vector3d.Z,h(1)));
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% %     rot_x=rotation.byAxisAngle(vector3d.X,55*degree);
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
%     rot_mat=rot_x*rot_z;
    if isnan(phi2)
        rot_mat=rotation.byEuler(0,Phi,0);
    else
        rot_mat=rotation.byEuler(0,Phi,phi2);
    end
    
end
% function plot_along_plane_normal(atoms,h)
% end

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

    [x_ax,y_ax,z_ax]=get_base_vecs(h);
      [x_ax,y_ax,z_ax]=get_base_vecs(h);
%     mults=zeros(n_cells_x,n_cells_y,n_cells_z,3);
    at_cord=[at.all_atoms.x at.all_atoms.y at.all_atoms.z];
   
    for i=start_cell_x:n_cells_x-1
        for j=start_cell_y:n_cells_y-1
            for k=start_cell_z:n_cells_z-1
    %             mults(i+1,j+1,k+1,:)=[i*x_ax + j*y_ax + k*z_ax]
                all_rep=vertcat(all_rep,at_cord+[i*x_ax + j*y_ax + k*z_ax]);
                labels=vertcat(labels, at.labels);
                elements=vertcat(elements, at.elements);
                colors_label=vertcat(colors_label, at.colors_label);
                colors_elements=vertcat(colors_elements, at.colors_elements);
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
    %%
    
end


% % function layers=get_intersections(atoms_in,h,varargin)
% %% this funciton is too full righ tnow, bu t i was lazy
% 
% %     
%     n_cells=1;
% %     h=Miller(1,0,-1,10,cs);
% %     corners=get_lattice_corners(h,'n_cells',n_cells)
%     n=[h.x;h.y;h.z];  %% we use column vectors here
% %     n=[-1;1;10]
%     n=n./norm(n);
%     [x_ax,y_ax,z_ax]=get_base_vecs(h);
%     l=[x_ax;y_ax;z_ax]';
%     %
% %     nLen_non=[x_ax*n, y_ax*n,z_ax*n]
% %     norm_n=norm(n)
% %         temp=norm([norm(x_ax)*n,norm(y_ax)*n,norm(z_ax)*n])
% %     n=n./norm(n);
%     diag=x_ax+y_ax+z_ax;   
%     nLen=norm([diag*n]);%%gives the lattice paramters if e.g. n=[0;0;1], so it is the "extend" of the Miller in the unit cell
%     %%
%     slices=5; % meaning one unit cell
%     for i=1:size(varargin(:))
%         if strcmp(varargin{i},'slices')
%             slices=varargin{i+1};
%         elseif strcmp(varargin{i},'n_cells')
%             n_cells=varargin{i+1};
%         else
%             continue;
%         end
%         
%     end
%     %%
% %     figure 
% %     hold on
% %     axis equal;
% %     
% %%
%     atoms=rep_atoms(atoms_in,n_cells,h);
% %     all_atoms=atoms.all_atoms;
%     %%
%     atom_cords=[atoms.all_atoms.x,atoms.all_atoms.y,atoms.all_atoms.z];
%     color=get_colors_from_list(slices+1);
% %     for i=1:slices
% %         disp(i)
% %         hold on
%     edges=get_lattice_edges(h(1));
% %         for i=1:size(edges,3)
% %             plot3(edges(1,:,i),edges(2,:,i),edges(3,:,i),'black')
% %         end
%     count=1;
%     atoms_return=struct();
%     %c
% %     figure
% %     hold on
%     for i=-slices:slices
%         col=color{abs(i)+1};
%         try
%             intersections_0=get_plane_intersections(h,nLen*i/(slices+1),varargin{:});
%             intersections_1=get_plane_intersections(h,nLen*(i+1)/(slices+1),varargin{:});
%             intersections_both=[intersections_0 intersections_1];
%             bounds=alphaShape(intersections_both(1,:)',intersections_both(2,:)',intersections_both(3,:)',Inf);
%             selection=bounds.inShape(atom_cords);
%             temp=atom_cords(selection,:);
% %             scatter3(temp(:,1),temp(:,2),temp(:,3))
% %             p=fill3(intersections_0(1,:),intersections_0(2,:),intersections_0(3,:),col(1:3));
% %             p=fill3(intersections_1(1,:),intersections_1(2,:),intersections_1(3,:),col(1:3));
%             %
%             layers(count).intersections_0=intersections_0;
%             layers(count).intersections_1=intersections_1;
% %             atoms_return(count).labels=atoms.labels(selection);
% %             atoms_return(count).elements=atoms.elements(selection);
% %             
% %             atoms_return(count).colors_elements=atoms.colors_elements(selection,:);
% %             atoms_return(count).colors_label=atoms.colors_label(selection,:);
% % %             scatter3(intersections_both(1,:),intersections_both(2,:),intersections_both(3,:))
% %             atoms_return(count).all_atoms.x=atoms.all_atoms.x(selection);
% %             atoms_return(count).all_atoms.y=atoms.all_atoms.y(selection);
% %             atoms_return(count).all_atoms.z=atoms.all_atoms.z(selection);
%             count=count+1;
%         catch
%             warning('slice probably outside of requested unit cell')
%         end
%         %%
%         
% 
%         %%
%     end
%     %%
% % p0=n*dist;
% % l0=corners;
% % l=[base_vecs.x;base_vecs.y;base_vecs.z];
%     
%     
% end
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

    dist= a.*(x-px)+b.*(y-py)+c.*(z-pz)/sqrt(a^2+b^2+c^2);
    %%
    %%https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane
    %%
    
end

function atoms_checked=check_zeros(atom_cords)
    if ~any(atom_cords==0)
        atoms_checked=atom_cords;
    elseif sum(atom_cords==0)==3
        atoms_checked=[[0;0;0;],[0;0;1],[0;1;0],[1;0;0],[1;1;0],[1;0;1],[0;1;1],[1;1;1]];
    elseif sum(atom_cords==0)==2
        atoms_checked=[]
        index=find(atom_cords==0);
        for i=0:1
            for j=0:1
                atom_cords_in=atom_cords;
                atom_cords_in(index(1))=i;
                atom_cords_in(index(2))=j;
                atoms_checked=[[atoms_checked], [atom_cords_in]];
            end
        end
    elseif sum(atom_cords==0)==1
        atoms_checked=[]
        index=find(atom_cords==0);
        for i=0:1
                atom_cords_in=atom_cords;
                atom_cords_in(index(1))=i;
                atoms_checked=[[atoms_checked], [atom_cords_in]];
        end

%         atoms_checked=[[0;0;0;],[0;0;1],[0;1;0],[1;0;0],[1;1;0],[1;0;1],[0;1;1],[1;1;1]]

    end
    %%
%     disp(atom_cords);
end