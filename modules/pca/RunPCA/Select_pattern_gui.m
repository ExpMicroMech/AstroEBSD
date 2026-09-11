function [OutputData]=Select_pattern_gui(InputData)
    %% get data from input string
    
    
    MapInfo=InputData.MapInfo;
%     PhaseData=InputData.PhaseData;
    InputUser=InputData.InputUser;
    ebsd=InputData.EBSD_map;
    MicroscopeData=InputData.MicroscopeData;


    %% make a figure and set up the axes
    
    f = figure('Visible','on','Position',[100,20,1200,750],'Color','w');
    set(gcf, 'units','normalized')
     axis off;
%      hold on
     for i=2:size(ebsd.mineralList,2)

        mineral_name=ebsd.mineralList{i}
        try
        ipfKey = ipfTSLKey(ebsd(mineral_name));
        ipfKey.inversePoleFigureDirection = vector3d.X;
        colors_X = ipfKey.orientation2color(ebsd(mineral_name).orientations);
        ipfKey.inversePoleFigureDirection = vector3d.Y;
        colors_Y = ipfKey.orientation2color(ebsd(mineral_name).orientations);
        ipfKey.inversePoleFigureDirection = vector3d.Z;
        colors_Z = ipfKey.orientation2color(ebsd(mineral_name).orientations);
        catch
            ipfKey = ipfTSLKey(ebsd(mineral_name));
            ipfKey.inversePoleFigureDirection = vector3d.X;
            colors_X = ipfKey.orientation2color(ebsd(mineral_name).orientations);
            ipfKey.inversePoleFigureDirection = vector3d.Y;
            colors_Y = ipfKey.orientation2color(ebsd(mineral_name).orientations);
            ipfKey.inversePoleFigureDirection = vector3d.Z;
            colors_Z = ipfKey.orientation2color(ebsd(mineral_name).orientations);
            warning('mineral names not loaded')
        end
        if i==2
            h_ebsd_map_ipfx = axes('Units','normalized','Position',[0.05 0.02 0.25 0.27],'Parent',f);
            h_ebsd_map_ipfy = axes('Units','normalized','Position',[0.27 0.02 0.25 0.27],'Parent',f);
            h_ebsd_map_ipfz = axes('Units','normalized','Position',[0.0 0.35 0.5 0.6],'Parent',f);
            h_ebsd_map_ipf_key = axes('Units','normalized','Position',[0.5 0.68 0.1 0.1],'Parent',f);
        end
            hold on 
            plot(ebsd(mineral_name),colors_X,'micronbar','off','parent',h_ebsd_map_ipfx);
            title('IPF X')
            plot(ebsd(mineral_name),colors_Z,'micronbar','off','parent',h_ebsd_map_ipfz);
            title('IPF Z')
        if isfield(ebsd.prop,'component')


            hold on
            plot(ipfKey,'parent',h_ebsd_map_ipf_key)
        else
            hold on
            plot(ebsd(mineral_name),colors_Y,'micronbar','off','parent',h_ebsd_map_ipfy);
            title('IPF Y')
%             hold on
        end

        if isfield(ebsd.prop,'component')
           plot(ebsd(mineral_name),ebsd(mineral_name).prop.component,'micronbar','off','parent',h_ebsd_map_ipfz);
           title('PCA components')
           colormap(h_ebsd_map_ipfz, "colorcube");
           [component_center_coordinates, component_names] = find_centerOfComponent(ebsd)
           text(component_center_coordinates(:,1), component_center_coordinates(:,2), int2str(component_names(:,1)));
        else
            % plot(ebsd(ebsd.phaseId==i),colors_Z,'micronbar','off','parent',h_ebsd_map_ipfz);
            plot(ebsd(ebsd.phaseId==i),ebsd(ebsd.phaseId==i).prop.iq,'micronbar','off','parent',h_ebsd_map_ipfy);
            colormap gray
            title('IQ')
%             hold on
try
            plot(ipfKey,'parent',h_ebsd_map_ipf_key)
catch
    warning('not able to plot ipf key')
end
        end
% <<<<<<< HEAD
     %    hold on
     %    % plot(ebsd(mineral_name),colors_Y,'micronbar','off','parent',h_ebsd_map_ipfy);
     %    plot(ebsd(mineral_name),ebsd(mineral_name).prop.imagequality,'micronbar','off','parent',h_ebsd_map_ipfy);
     %    title('IPF Y')
     %    hold on
     %    plot(ebsd(mineral_name),colors_Z,'micronbar','off','parent',h_ebsd_map_ipfz);
     %    title('IPF Z')
     %    hold on
     %    plot(ipfKey,'parent',h_ebsd_map_ipf_key)
     % end
% =======
% >>>>>>> eda632a50ef20e5012d186ec265bd0f5a25b88bf

     end
% 
%     plotR(ebsd,ebsd.or,'parent',h_ebsd_map)
    x_cord=10;
    y_cord=10;
    % [MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=breadh5oina(InputUser);
    pnum=find(MapInfo.MapData.XBeam==x_cord & MapInfo.MapData.YBeam==y_cord);
    disp(pnum)
    % components=EBSPData.tile.map_reshaped_original';
    % components=components(:);
    % pattern_index=components(pattern_number);
    %%
    pat=bReadEBSP(MapInfo.EBSPData,[pnum]);

    h_ebsd_pattern = axes('Units','normalized','Position',[0.5 0.1 0.5 0.5],'Parent',f);
%     hold on
    imagesc(pat)
    axis xy;
    axis equal;
    axis off;
    colormap(h_ebsd_pattern, 'gray');
    axis(h_ebsd_pattern,'tight');
    %% ToDo: check, which crytsal strucutre we have and make it dependent on this, which to plot
%     if 
    


%% set up our data cursor

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'off';
dcm.UpdateFcn = @displayCoordinates;
% nexttile
% hPlot.ButtonDownFcn = @(h,e) disp(e.IntersectionPoint);

xsep=80; xstart=290; xwid=70; ystart=100; ysep=-18; yhig=15;
% x_start=0.1
% xsep=0.1
% ystart=0.1
% xwid=0.1
% yhig=0.1

%% set up all the buttons
% h_dipslay_pattern=uicontrol('style','pushbutton','string', ...
%     'Display Pattern',...
%     'position',[xstart+xsep*5 ystart xwid yhig],'BackgroundColor','w');
% h_display_pattern.Units = 'normalized';
% h_display_pattern.Callback=@Update_pattern;

h_coords=uicontrol('style','text','string', ...
    sprintf('coords'),...
    'position',[xstart+xsep*0 ystart+ysep xwid yhig],'BackgroundColor','w','units','normalized');
% h_coord.Value=x_cord

h_coords_real=uicontrol('style','text','string', ...
    sprintf('coords, µm'),...
    'position',[xstart+xsep*0 ystart+ysep*2 xwid yhig],'BackgroundColor','w','units','normalized');
% h_coord.Value=x_cord

h_x_coord=uicontrol('style','text','string', ...
    sprintf('%.3f',round(x_cord./MicroscopeData.XSTEP)),...
    'position',[xstart+xsep*1 ystart+ysep xwid yhig],'BackgroundColor','w','units','normalized');
h_x_coord.Value=round(x_cord./MicroscopeData.XSTEP);

h_y_coord=uicontrol('style','text','string', ...
    sprintf('%.3f', round(y_cord./MicroscopeData.XSTEP)),...
    'position',[xstart+xsep*2 ystart+ysep xwid yhig],'BackgroundColor','w','units','normalized');

h_y_coord.Value=round(y_cord./MicroscopeData.YSTEP);


h_x_coord_real=uicontrol('style','text','string', ...
    sprintf('%.3f',x_cord),...
    'position',[xstart+xsep*1 ystart+ysep*2 xwid yhig],'BackgroundColor','w','units','normalized');
h_x_coord_real.Value=x_cord;

h_y_coord_real=uicontrol('style','text','string', ...
    sprintf('%.3f', y_cord),...
    'position',[xstart+xsep*2 ystart+ysep*2 xwid yhig],'BackgroundColor','w','units','normalized');

h_y_coord_real.Value=y_cord;

%% Pole figure plots

    h_pole_figure_1 = axes('Units','normalized','Position',[0.5 0.6 0.2 0.2],'Parent',f);
    h_pole_figure_2 = axes('Units','normalized','Position',[0.7 0.6 0.2 0.2],'Parent',f);
%     hold on
        mineral_name=ebsd.mineralList{2}
    ori=ebsd(mineral_name).orientations
    
    m1=Miller(0,0,1,ebsd(mineral_name).CS)
    plotPDF(ori,m1,'Parent',h_pole_figure_1)
    title(char(m1))
    axis(h_pole_figure_1,'tight');

    m2=Miller(1,0,0,ebsd(mineral_name).CS)
    plotPDF(ori, m2, 'Parent',h_pole_figure_2)
    title(char(m2))
    axis(h_pole_figure_1,'tight');


%%
function  [txt,x,y] = displayCoordinates(~,eventdata)
    x = eventdata.Position(1);
    y = eventdata.Position(2);
    h_x_coord_real.Value=x;
    h_y_coord_real.Value=y;
    
    Update_coordinates;
    Update_pattern;
    Update_PoleFigure;
    txt = ['(' num2str(x) ', ' num2str(y) ')'];
end

function Update_pattern(~,eventdata)
   set(f, 'CurrentAxes', h_ebsd_pattern)
% <<<<<<< HEAD
%     % disp([h_x_coord.Value,h_y_coord.Value])
%     pnum=find(MapInfo.MapData.XBeam==h_x_coord.Value & MapInfo.MapData.YBeam==h_y_coord.Value);
%     disp(pnum)
% =======
    pnum=find(MapInfo.MapData.XBeam==h_x_coord.Value & MapInfo.MapData.YBeam==h_y_coord.Value   );
% >>>>>>> eda632a50ef20e5012d186ec265bd0f5a25b88bf
    pat=bReadEBSP(MapInfo.EBSPData,[pnum]);
    if isfield(InputData,'SettingsCor')
      pat = EBSP_BGCor(pat,InputData.SettingsCor);
    end
    imagesc(pat)
    axis equal;
    axis xy;
    axis off;
    colormap(h_ebsd_pattern,gray);
    axis(h_ebsd_pattern,'tight');
    disp([h_x_coord.Value,h_y_coord.Value])

end

function Update_PoleFigure(~,eventdata)
    set( f, 'CurrentAxes', h_pole_figure_1);
    ori=ebsd(h_x_coord.Value,h_y_coord.Value).orientations;
%     delete(pf1)
    disp(h_pole_figure_1.Children);
    disp(h_pole_figure_1.Children(5).Children())
    if size(h_pole_figure_1.Children(5)>=2)
        delete(h_pole_figure_1.Children(5).Children(1))
    end
    hold off;
%     m1=Miller(0,0,1,ebsd(mineral_name).CS)
    plotPDF(ori, m1,'Parent',h_pole_figure_1)
    title(char(m1))
    axis(h_pole_figure_1,'tight');

    set( f, 'CurrentAxes', h_pole_figure_2);
%     m2=Miller(1,0,0,cs,ebsd(mineral_name).CS)
%     delete(pf2)
    hold off;
    if size(h_pole_figure_2.Children(5)>=2)
         delete(h_pole_figure_2.Children(5).Children(1))
    end
    plotPDF( ori, m2, 'Parent',h_pole_figure_2)

    title(char(m2))
    axis(h_pole_figure_1,'tight');
    axis( h_pole_figure_1, 'tight');

end
    function [center,unique_comps_new] = find_centerOfComponent(ebsd)
    
    components = ebsd.prop.component;
    unique_comps = unique(components);
    center = zeros(size(unique_comps,1),2);
    unique_comps_new=zeros(size(unique_comps,1),2);
    unique_comps_new(:,1)=unique_comps;

    x = ebsd.prop.x;
    y = ebsd.prop.y;

    for i=1:size(unique_comps,1)
        center(i, 1) = mean( x( components == unique_comps(i) ) );
        center(i, 2) = mean( y( components == unique_comps(i) ) );
        unique_comps_new(i,2) = sum( components == unique_comps(i) );
    end

end
function Update_coordinates(~,eventdata)
    x=h_x_coord_real.Value;
    y=h_y_coord_real.Value;
%     y=eventdata.Position(2)
    disp(floor(x./MicroscopeData.XSTEP));
    h_x_coord.Value=floor(x./MicroscopeData.XSTEP+MicroscopeData.XSTEP/2);
    h_y_coord.Value=floor(y./MicroscopeData.YSTEP+MicroscopeData.XSTEP/2);
    h_x_coord.String=num2str(h_x_coord.Value);
    h_y_coord.String=num2str(h_y_coord.Value);
% <<<<<<< HEAD
%     h_x_coord_real.String=num2str(y);
%     h_y_coord_real.String=num2str(x);
% =======
    h_y_coord_real.String=num2str(y);
    h_x_coord_real.String=num2str(x);
% >>>>>>> eda632a50ef20e5012d186ec265bd0f5a25b88bf
    guidata(h_x_coord);
    guidata(h_y_coord);
    guidata(h_x_coord_real);
    guidata(h_y_coord_real);
        
end
end