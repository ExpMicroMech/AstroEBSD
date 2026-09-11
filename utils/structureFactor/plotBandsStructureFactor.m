function pattern =plotBandsStructureFactor(EBSD_geometry,m,r,InputData,varargin)
    theta=(0:0.5:360)*degree;
    plotColor = "";
    imagePattern=0;
    for i = 1:numel(varargin(:))
        if strcmp(varargin{i},'imagePattern')
            imagePattern = 1
        end
        if strcmp(varargin{i},'color')
            plotColor = varargin{i+1};
        end
    end
    if ~isfield(InputData,'maxBNum')
        InputData.maxBNum = 10

    end
    InputData.maxBNum = min(InputData.maxBNum,size(m,1))
    if size(m,1) < 8
        colors = colorcube(9);
    else
        colors = colorcube(size(m,1));
    end
    % colors = colors(5:end,:)
    colors = colors*0.7;
    col_count = 0;
    % fix that we plot too many bands, and transparncy gets gone
    pattern = zeros(size(EBSD_geometry.xpts_screen));
    % for i=1+(size(m,1)-InputData.maxBNum):size(m,1) %% cut of bands from the other side. this way the highest intensity band is plotted last and non top.
    for i = 1:InputData.maxBNum
    %%
        m_sym=m(i).symmetrise('unique','antipodal');
        %%
        % col_count=1;
        isPlotted = 0;
        for j=1:size(m_sym,1)
            col_plot=colors(i,:);
            Iplot = InputData.IPlot(i);
            col_plot = [Iplot Iplot Iplot];
            % col_plot=[col_plot 1];
            % col_plot = [ F_AbsSquared(i) F_AbsSquared(i) F_AbsSquared(i)];
            % if j==1
            %     plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
            % end
            % lw=F_AbsSquared(i);
            lw =0.1;
            for k=1:3
                %% depending on case, plot band center or band edges
                if k==1
                % omega=pi/2;
                    continue;
                elseif k==2
                    omega = pi/2-asin(InputData.lambda_p/(2*m_sym(j).dspacing));
                    omegaFirst = omega;
                else
                    omega = pi/2+asin(InputData.lambda_p/(2*m_sym(j).dspacing));
                    % omega
                end
                plane=circlesMTEXModified(r*m_sym(j)./norm(r*m_sym(j)),theta,omega);
            
            hold on
            % m_sym(j)
            gn_coords=[plane.x./plane.z, plane.y./plane.z];
            % only select coordinates that are in front
            gn_coords=gn_coords(plane.z>0.001,:);
            % get rid of everything that is bigger than the pattern
            gn_coords=find_limits(gn_coords,EBSD_geometry);
            difference = gn_coords(1:end-1,:)-gn_coords(2:end,:);
            [max1,ind1] = max(abs(difference(:,1)),[],1);
            [max2,ind2] = max(abs(difference(:,2)),[],1);
            if max1>max2
                ind = ind1;
            else
                ind = ind2;
            end
            try
                if max1>mean(abs(difference(:)))*3 || max2>mean(abs(difference(:)))*3
                    gn_coords = vertcat(gn_coords(ind+1:end,:),gn_coords(1:ind-1,:));
                end
            catch
                warning("comparison did not work");
                % print([max1, max2]);
            end
            %plot the bands
            if k ==2
                gn_coords_first = gn_coords;
            end
            if ~isempty(gn_coords)
                hold on
                % col_plot=InputData.cset(i,:)*(1-(col_count/size(m_sym,1)));
                
                % col_plot=[col_plot 0.5];
                if ~ strcmp(plotColor,"")
                    col_plot = plotColor;
                end
                if strcmp(plotColor,"fhkl2")
                    I_col = 1-InputData.IPlot(i);
                    col_plot =[I_col I_col I_col];
                end
                if strcmp(plotColor,"fhkl")
                    I_col = 1-sqrt(InputData.IPlot(i));
                    col_plot =[I_col I_col I_col];
                end
                if imagePattern == 1
                    I_col = 1-I_col;
                    if ~isPlotted == 1
                    % else
                        isPlotted =1;
                        col_count=col_count+1;
                    end

                    continue
                end
                if k==1
                   % plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
                elseif k==2
                      if isPlotted == 1
                        l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw,'HandleVisibility','off');
                      else
                        isPlotted =1;
                        l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m_sym(j).project2FundamentalRegion));
                        col_count=col_count+1;
                      end
                else
                    l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
                end

                % uistack(l, 'bottom');
%%
                % % l=patchline()
                % % l.lh=[0,0,0,0.5]
                %  if k ==2
                %     gn_old = gn_coords
                %  elseif k==3 
                %     s = patch([gn_old(:,1) gn_coords(:,1)], [gn_old(:,2) gn_coords(:,2)],col_plot)
                %     alpha(s, col_plot)
                %  end
             end
          %%
  
        %%
            end
                       hold on
        %%
        if strcmp(plotColor,"")
            continue
        end
         if ~isempty(gn_coords)
             %%
             % pattern = zeros(size(EBSD_geometry.xpts_screen));

            xSmaller1 = EBSD_geometry.xpts_screen(:) < gn_coords_first(:,1)';
            xSmaller2 = EBSD_geometry.xpts_screen(:) < gn_coords(:,1)';
            
            ySmaller1 = EBSD_geometry.ypts_screen(:) < gn_coords_first(:,2)';
            ySmaller2 = EBSD_geometry.ypts_screen(:) < gn_coords(:,2)';
            
            xLarger1 = ~xSmaller1;
            yLarger1 = ~ySmaller1;
            xLarger2 = ~xSmaller2;
            yLarger2 = ~ySmaller2;
            % xLarger1  = EBSD_geometry.xpts_screen(:) > gn_coords_first(:,1)';
            % xLarger2  = EBSD_geometry.xpts_screen(:) > gn_coords(:,1)';
            % 
            % yLarger1  = EBSD_geometry.ypts_screen(:) > gn_coords_first(:,2)';
            % yLarger2  = EBSD_geometry.ypts_screen(:) > gn_coords(:,2)';

            % This goes over the various possibilities where we have
            % diagonals running from top lef to bottom right (case 1 and2)
            % or top right to bottom lef (2 and 1). 
            % further which is further right/ left is distinguished.
            if mean(gn_coords_first(:,1)) >= mean(gn_coords(:,1)) & mean(gn_coords_first(:,2)) >= mean(gn_coords(:,2)) %Case1 
                selection1 = any( xLarger1   & yLarger1 ,2);
                selection2 = any( xSmaller2  & ySmaller2 ,2);


            elseif mean(gn_coords_first(:,1)) < mean(gn_coords(:,1)) & mean(gn_coords_first(:,2)) < mean(gn_coords(:,2)) %Case 2

                selection1 = any( xSmaller1  & ySmaller1 ,2);
                selection2 = any( xLarger2   & yLarger2 ,2);
            % else 
            %     continue
            % end
            elseif  mean(gn_coords_first(:,1)) < mean(gn_coords(:,1)) & mean(gn_coords_first(:,2)) >= mean(gn_coords(:,2)) %Case 3

                selection1 = any( xSmaller1  & yLarger1 ,2);
                selection2 = any( xLarger2   & ySmaller2 ,2);
            % 
            else
                % continue
            % 
                selection1 = any( xLarger1  & ySmaller1 ,2);
                selection2 = any( xSmaller2   & yLarger2 ,2);

            end
            
            % 
            % selection1 = any( (EBSD_geometry.xpts_screen(:) < gn_coords(:,1)'...
            %         & EBSD_geometry.ypts_screen(:) < gn_coords(:,2)') ,2);
            % selection2 =any( (EBSD_geometry.xpts_screen(:) > gn_coords_first(:,1)'...
            %         & EBSD_geometry.ypts_screen(:) > gn_coords_first(:,2)') ,2); % should work for bands running from top left to botoom right (first band left right

          %   selection3 =any( (EBSD_geometry.xpts_screen(:) > gn_coords(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) > gn_coords(:,2)') ,2);
          %   selection4 =any( (EBSD_geometry.xpts_screen(:) < gn_coords_first(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) < gn_coords_first(:,2)') ,2);
          % 
          % selection3 =any( (EBSD_geometry.xpts_screen(:) > gn_coords(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) > gn_coords(:,2)') ,2);
          %   selection4 =any( (EBSD_geometry.xpts_screen(:) < gn_coords_first(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) < gn_coords_first(:,2)') ,2);
          % 
          %     selection3 =any( (EBSD_geometry.xpts_screen(:) > gn_coords(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) > gn_coords(:,2)') ,2);
          %   selection4 =any( (EBSD_geometry.xpts_screen(:) < gn_coords_first(:,1)'...
          %           & EBSD_geometry.ypts_screen(:) < gn_coords_first(:,2)') ,2);
          %   % pattern()=0.5;
            % selection = any(
            % I_col
            pattern(~selection1 & ~ selection2) = pattern(~selection1 & ~ selection2) + I_col;%pattern(selection1)1;
            % pattern(selection2)= 1%pattern(selection2);
    
            % pattern(:,end)=1
            % pattern(end,:)=1
            % pattern(:,end)=1
            % pattern(end,:)=1
            
            % % figure
            % pPattern(pattern,EBSD_geometry);
            % % imagesc(pattern)
            % pattern = zeros(size(EBSD_geometry.xpts_screen));
            % pattern(~selection3 & ~selection4)= pattern(~selection3 & ~selection4) +I_col;
            %  figure
            % pPattern(pattern,EBSD_geometry);
            % %%
            % hold on
            % scatter(gn_coords(:,1),gn_coords(:,2))
            % scatter(gn_coords_first(:,1),gn_coords_first(:,2))
            % colormap gray
            % clim([0,1])
            %     % patch([gn_coords(:,1) gn_coords(:,2)])
                % 
                % try
                % catch:
                % end
    
                % if ~isempty(gn_coords) & isPlotted == 0;
                %     col_count=col_count+1;
                % end
         end
        end
        
        if col_count >InputData.maxBNum
            break;
        end
    %     disp([i j])
        
    end
    %%
    % uistack(img, 'bottom');
    % legend('Location','eastoutside')
    xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
    ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
        % legend('FontSize',8)
        % axis off
    figure
    pPattern(pattern,EBSD_geometry);
    colorbar()
    disp(i)
end