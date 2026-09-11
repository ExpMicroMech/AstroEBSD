function intersections = get_plane_intersections(h,dist,varargin) %%https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
% base_vecs=h.CS.axes;
% % base_vecs;


    extend=1;
    inf_e=0;
    for i=1:size(varargin(:))
        if strcmp(varargin{i},'extend')
            extend=varargin{i+1};
        elseif strcmp(varargin{i},'inf_e')
            inf_e=1;
            continue;
        end
    end
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
