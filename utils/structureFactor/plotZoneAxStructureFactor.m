function plotZoneAxStructureFactor(EBSD_geometry,m,r,InputData,varargin)
theta=(0:1:360)*degree;
plotColor = "";
for i = 1:size(varargin,1)
    if strcmp(varargin{i},'color')
        plotColor = varargin{i+1};
    end
end
if ~isfield(InputData,'maxBNum')
        InputData.maxBNum = 10
end
    if size(m,1) < 8
        colors = colorcube(9);
    else
        colors = colorcube(size(m,1));
    end
    % colors = colors(5:end,:)
    colors = colors*0.7;
    col_count = 0;
    % fix that we plot too many bands, and transparncy gets gone
    for i=1:size(m,1)
        %%
        m_sym=m(i).symmetrise();%'unique','antipodal');
        %%
        % col_count=1;
        isPlotted = 0;
        for j=1:size(m_sym,1)
            col_plot=colors(i,:);
            % Iplot = InputData.IPlot(i);
            % col_plot = [Iplot Iplot Iplot];
            % col_plot=[col_plot 1];
            % col_plot = [ F_AbsSquared(i) F_AbsSquared(i) F_AbsSquared(i)];
            % if j==1
            %     plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
            % end
            % lw=F_AbsSquared(i);
            lw =0.1;
            % for k=1:3
                %% depending on case, plot band center or band edges
                % plane=circlesMTEXModified(r*m_sym(j)./norm(r*m_sym(j)));%,theta,omega);
            plane = r*m_sym(j)./norm(r*m_sym(j));
            %%
            hold on
            % m_sym(j)
            gn_coords=[plane.x./plane.z, plane.y./plane.z];
            % only select coordinates that are in front
            gn_coords=gn_coords(plane.z>0.001,:);
            % get rid of everything that is bigger than the pattern
            gn_coords=find_limits(gn_coords,EBSD_geometry);
           
            %plot the bands
            if ~isempty(gn_coords)
                % col_plot=InputData.cset(i,:)*(1-(col_count/size(m_sym,1)));
                
                % col_plot=[col_plot 0.5];
                if ~ strcmp(plotColor,"")
                    col_plot = plotColor;
                end
                   scatter(gn_coords(:,1),gn_coords(:,2),150,'filled','MarkerEdgeColor','black','MarkerFaceColor',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
                    len = EBSD_geometry.x_gn_max-EBSD_geometry.x_gn_min;
                   % dim = [gn_coords(:,1)+0.02*len, gn_coords(:,2)-0.02*len,0.1*len,0.03*len]
                   text(gn_coords(:,1)+0.02*len, gn_coords(:,2)-0.02*len,...
                       char(m_sym(j)),'BackgroundColor','white');
            end
               hold on
               %%
        if col_count >InputData.maxBNum
            break;
        end
    %     disp([i j])
        end
    end
    % uistack(img, 'bottom');
    % legend('Location','eastoutside')
    xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
    ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
        % legend('FontSize',8)
        % axis off
end