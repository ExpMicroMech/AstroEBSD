function plotAtomLayersToScale(atoms,sections,millerPlot,col,sliceNumber,varargin)
   % col=color{i};
   % figure,ax = plt.subplots(2,2)
   % OFFSET in multiples of the lattice spacing
   offset = 0
   plotAtomPosition = 0
     for i=1:size(varargin(:))
        if strcmp(varargin{i},'offset')
            offset=varargin{i+1};
        end    
        if strcmp(varargin{i},'plotAtomPosition')
            plotAtomPosition = 1;
        end
    end
   figure
%    for i=1:sz
   % ax=nexttile(t3);
    hold on
    rot_mat=get_h_rot(millerPlot);
    % rot_mat = rotation.byEuler(0,90*degree,0*degree)
    % rot_mat=rotation.byEuler(180*degree,0,0)*rot_mat
    temp=vector3d(atoms.all_atoms.x,atoms.all_atoms.y, atoms.all_atoms.z);
    %
    % newTemp = rot_mat.matrix'*(temp.xyz)';

    % figure
    % scatter3(newTemp(1,:),newTemp(2,:),newTemp(3,:))
    % axis equal;
    %
    % % tempVec = [atoms.all_atoms.x,atoms.all_atoms.y, atoms.all_atoms.z]
    temp= rot_mat*temp;
    % figure
    % scatter3(temp.x,temp.y,temp.z)
    % axis equal;
    %
    % temp = vector3d(temp(1,:),temp(2,:,2),temp(:,3))
    axesRot = rot_mat*millerPlot.CS.axes;
    
    %
    [~,uniqueIndices,~] =  unique([atoms.all_atoms.x,atoms.all_atoms.y,atoms.all_atoms.z],'rows')
%     figure
    
    axis equal;
    
%     axis equal;
   
% end
    edges=get_lattice_edges(millerPlot);
    % [x_ax,y_ax,z_ax]=get_base_vecs(millerPlot); 
    x_ax=[axesRot(1).x axesRot(1).y axesRot(1).z];
    y_ax=[axesRot(2).x axesRot(2).y axesRot(2).z];
    z_ax=[axesRot(3).x axesRot(3).y axesRot(3).z];
 
    first=vector3d(squeeze(edges(:,1,:)));
    second=vector3d(squeeze(edges(:,2,:)));
    first=vector3d([0 0 0], [0 0 0], [0 0 0])
    second=vector3d(x_ax,y_ax,z_ax)

    %
    cornersTop=vector3d(sections.intersections_left);
    corners = vector3d(sections.intersections_center);
    cornersBottom=vector3d(sections.intersections_right);
    %
    %
    dot_direction=round(dot(repmat(millerPlot.xyz',[1 size(cornersTop,2)]),cornersTop.xyz'),3);
    %
    index_min=find(min(dot_direction)==dot_direction);
    index_min=index_min(1);
    % start_corner=corners(:,index_min).xyz;

    corners_plotTop=(rot_mat)*cornersTop;
    corners_plotBottom  = (rot_mat) * cornersBottom;
    corners_plot = (rot_mat) * corners;
    


    first=rot_mat*first;
    second=rot_mat*second;
    % col=color{i};
     p=patch(corners_plotTop.x,corners_plotTop.y,corners_plotTop.z,col(1:3),'LineStyle','--','EdgeColor',[0.5 0.5 0.5]);%col(1:3)
     p.FaceAlpha = 0;
     p=patch(corners_plotBottom.x,corners_plotBottom.y,corners_plotBottom.z,col(1:3),'LineStyle','--','EdgeColor',[0.5 0.5 0.5]);%col(1:3)
     % fill3(sections(i).intersections_center_plot(1,:),sections(i).intersections_center_plot(2,:),sections(i).intersections_center_plot(3,:),col(1:3),'LineStyle','--','EdgeColor',[0.3 0.3 0.3])
     p.FaceAlpha=0;
     p=patch(corners_plot.x,corners_plot.y,corners_plot.z,col(1:3),'LineStyle',':','EdgeColor',[0.5 0.5 0.5]);%col(1:3)
     p.FaceAlpha = 0.1;
     corners = get_lattice_corners(millerPlot,'n_cells',1);
    start_corner=corners(:,index_min);
     %
    plot3([0 x_ax(1)], [0 x_ax(2)], [0 x_ax(3)],'black', 'linewidth', 2)
    text(x_ax(1)+0.3,x_ax(2),x_ax(3), 'a','FontSize',14)
    plot3([0 y_ax(1)], [0 y_ax(2)], [0 y_ax(3)], 'black', 'linewidth', 2)
    text(y_ax(1)+0.3,y_ax(2),y_ax(3), 'b','FontSize',14)
    plot3([0 z_ax(1)], [0 z_ax(2)], [0 z_ax(3)], 'black', 'linewidth', 2)
    text(z_ax(1)+0.3,z_ax(2),z_ax(3), 'c','FontSize',14)
    % grid(ax,'off')
    grid off
    axis off
    
    h_cur=millerPlot;
    n_vec=[h_cur.x;h_cur.y;h_cur.z];
    n_vec=n_vec./norm(n_vec) ;

%     [x_ax,y_ax,z_ax]=get_base_vecs(h_cur);
%     diag=x_ax+y_ax+z_ax;   
%     nLen=norm([diag*n_vec]);%%gives the lattice paramters if e.g. n=[0;0;1], so it is the "extend" of the Miller in the unit cell
%     slices=nLen/h_cur.dspacing
%
    max_z_dev=0.5*atoms.dSpacing;
      
    % z_diff=get_dist_from_plane(n_vec,(sliceNumber-0.5+offset)*atoms.dSpacing*n_vec,  atoms.all_atoms);
    % z_diff = temp.z - (sliceNumber-0.5)*h_cur.dspacing;
    % z_diff = z_diff(:)
    z_diff = temp.z'-(sliceNumber-0.5+offset)*atoms.dSpacing;
    rInPlane = sqrt((atoms.radii).^2-(z_diff).^2);
    rInPlane = real(rInPlane);
    %
    rInPlane = rInPlane(uniqueIndices);
    plotCordsX  = temp.x(uniqueIndices);
    plotCordsY =  temp.y(uniqueIndices);
    plotCordsZ =  temp.z(uniqueIndices);
    rReal = atoms.radii;
    rReal = rReal(uniqueIndices);
    z_diff_sel = z_diff(uniqueIndices);
    col_new=rgb2hsv(atoms.colors_elements(:,1:3));
    numPoints = 20; 
    % radius = linspace(0,100,numrings+1);   %where the rings are to be drawn    
    theta = linspace(0,2*pi,numPoints);
    [THETA,RADIUS] = ndgrid(theta,rInPlane);
    [THETAREal,RADIUSReal] = ndgrid(theta,rReal);
    [X,Y] = pol2cart(THETA,RADIUS);
    [XR,YR] = pol2cart(THETAREal,RADIUSReal);
    % figure
    col_new(:,3)=0.7;
    col_new(:,2)=1;
    col_new(:,3)=(col_new(:,3) +0.25*(z_diff'/max_z_dev)');
    col_new(:,3) = (col_new(:,3))/max(col_new(:,3));
%     col_new(:,3)=col_new(:,3)./max(col_new(:,3));
    col_new(:,2)=(col_new(:,2) - 0.25*abs((z_diff'/max_z_dev)'));
    colSel = col_new(uniqueIndices,:)
    % colSel = rescale(colSel)
      
    filler = zeros(size(X(:,1)));
    for K = 1 : size(colSel)
       % xpoints = [X(:,K); flipud(X(:,K+1))];
       % ypoints = [Y(:,K); flipud(Y(:,K+1))];
       patch(X(:,K)+plotCordsX(K), Y(:,K)+plotCordsY(K),hsv2rgb(colSel(K,:)), 'EdgeColor', 'black','FaceAlpha',1);
       % [x,y,z]=sphere(10);
       % x=rReal(K)*x;
       % y = rReal(K)*y;
       % z = rReal(K)*z;
       % surf(x+plotCordsX(K),y+plotCordsY(K),z+plotCordsZ(K),'FaceColor',hsv2rgb(colSel(K,:)),'EdgeColor','None','FaceAlpha',0.5)
       hold on
       % p.faceAlpha = 1
    end
     for K = 1 : size(colSel)
       % xpoints = [X(:,K); flipud(X(:,K+1))];
       % ypoints = [Y(:,K); flipud(Y(:,K+1))];
       patch(XR(:,K)+plotCordsX(K), YR(:,K)+plotCordsY(K),repmat(z_diff_sel(K),size(X(:,K))), rescale(hsv2rgb(colSel(K,:))), 'EdgeColor', 'black','LineStyle','--','FaceAlpha',.15,'EdgeAlpha',.5);
       % p.FaceAlpha=0.5
       hold on
    end
    grid on
    axis on
    axis equal
    %%
    xN=    atoms.all_atoms.x(abs(z_diff)<0.6)
    yN=    atoms.all_atoms.y(abs(z_diff)<0.6)
    zN=    atoms.all_atoms.z(abs(z_diff)<0.6)
    % figure
    % scatter3(temp.x,temp.y,temp.z)
    % axis equal;
    % figure
    % scatter3(xN,yN,zN)
    % axis equal;
    %%
%     % min(col_new)
%     % max(col_new)
%     %
% 
% %     temp_colors=col_new(:,3)
% %     temp_colors(temp_colors<=0)=0;
% %     col_new(:,3)=temp_colors;
% %     temp_colors=col_new(:,2)
% %     temp_colors(temp_colors<=0)=0;
% %     col_new(:,2)=temp_colors;
%     % col_new(:,2)=col_new(:,3)./max(col_new(:,3));
%     % s=scatter(temp.x,temp.y,80,hsv2rgb(col_new),'filled','MarkerEdgeColor','black');
%     s = viscircles([temp.x',temp.y'],rInPlane)
% %     s.AlphaData=1-abs(z_diff/(millerPlot.dspacing*2));
% %     s.MarkerFaceAlpha='flat';
%     %
%     del_x=max(corners_plot.x)-min(corners_plot.x);
%     del_y=max(corners_plot.y)-min(corners_plot.y);
%     xl=[min(corners_plot.x)-del_x*0.1,max(corners_plot.x)+del_x*0.1];
%     yl=[min(corners_plot.y)-del_y*0.1,max(corners_plot.y)+del_y*0.1];
%     ax.XLim=xl;
%     ax.YLim=yl;
    %%
    % cols=unique(atoms.colors_elements(:,1:3),'rows');
% z_dists=-max_z_dev:0.01:max_z_dev;
% % nexttile
% 
% for i=1:size(cols,1)
% %     disp(i)
%     % ax=nexttile(t4);
%     box off
%     x=1:size(z_dists,2);
%     y=1:size(z_dists,2);
%     col_new=rgb2hsv(cols(i,:));
%     %
%     col_new(:,3)=0.7
%     col_new(:,2)=1;
%     z_col_3=col_new(3)+0.25*(z_dists/(max_z_dev)');
% %     z_col_3=z_col_3./max(z_col_3(:));%col_new(:,3)./max(col_new(:,3));
%     z_col_2=col_new(2)-0.25*abs((z_dists/(max_z_dev)'));
%     col_new=cat(1,repmat(col_new(1),size(z_dists,2),1)',z_col_2,z_col_3)';
%         % min(col_new)
%         % max(col_new)
%         %
% %         col_new(:,3)=(col_new(:,3)+0.25*(z_diff/(h(miller_number).dspacing))');
% %         col_new(:,3)=col_new(:,3)./max(col_new(:,3));
% %         col_new(:,2)=(col_new(:,2)-0.25*abs((z_diff/(h(miller_number).dspacing*4))'));
% %     s=scatter(x,y,30,hsv2rgb(col_new'),'filled')%,'MarkerSize',10);
%     custom_cmap=colormap(gca,hsv2rgb(col_new));
%     clim([min(z_dists),max(z_dists)]);
%     grid off
% % end
%     min(z_dists)*(millerPlot.dspacing)
%     c=colorbar()
%     set(gca,'XColor', 'none','YColor','none')
%     c.Location = 'southoutside';
% 
% % scatter(atoms_Out(i).colors_elements(:,1:3))
% 
%     axis off
%     elements=vertcat(atoms.labels);
%     keySet_element=unique(string(elements));
%     legend(keySet_element)
%     % text(0,0,keySet_element(i),'units','normalized')
%     clear ax;
%     % t4.Padding = 'tight';
%     % t4.TileSpacing = 'tight';
% 
%     if i==1
%         c.Label.String = 'z in Angstroms';
%     end
% end
end