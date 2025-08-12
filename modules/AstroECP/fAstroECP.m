function [Output_Data]=fAstroECP(Input_Data,ECP_Pat)
%This tool enables users to load an ECP and visualize this against a simulation
%
% Input_Data is a structure, with the following:
%
% optional - defaults will be used if these are not set
% PC_in = [PCx, PCy, PCz] in fractions of the width of the (square) ECP
%
% S_Start = [Sx, Sy, Sz] in degrees for the initial stage configuration
%
% eangs = [phi1, PHI, phi2] in degrees for the crystal orientation
%
% DeltaRs = [Drx, DrY, Drz] in degrees - the default step for the buttons
% to move the pattern around
%
% ECP_Loc = location of the experimental ECP on the disk, 
% this should be a TIF file

Output_Data=1;

%default the band coloring
if isfield(Input_Data,'cset')
else
    cset=[1,0,0; 0,1,0; 0,0,1; 0,1,1; 1,0,1;1,1,0];
    cset=[cset;cset*0.5];
    Input_Data.cset=cset;
end

%default the delta buttons
if isfield(Input_Data,'DeltaRs')
else
    Input_Data.DeltaRs=[0.5 0.5 2]; % delta values for controlling tilts etc, in degrees
end

% set some fields
if isfield(Input_Data,'PC_in')
    PC_start=Input_Data.PC_in;
else
    PC_start=[0.5,0.5,3.8];
end

if isfield(Input_Data,'Stage_in')
    S_start=Input_Data.Stage_in;
else
    S_start=[0 0 0];
end

if isfield(Input_Data,'eangs')
    eangs=Input_Data.eangs;
else
    eangs=[0 0 0];
end

if isfield(Input_Data,'DeltaRs')
    DeltaRs=Input_Data.DeltaRs;
else
    DeltaRs=[0.1 0.1 0.1];
end

if isfield(Input_Data,'ECP_Loc')
    ECP_pattern=imread(Input_Data.ECP_Loc);
    if size(ECP_pattern,3) == 3
        ECP_pattern=rgb2gray(ECP_pattern);
    end

    ECP_pattern=double(flipud(ECP_pattern));
    ECP_loaded=1;
else
    ECP_loaded=0;
end

%code in the rotation set if they are not specified
if isfield(Input_Data,'Rz')
else
Input_Data.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
end

if isfield(Input_Data,'Ry')
else
Input_Data.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
end

if isfield(Input_Data,'Rx')
else
Input_Data.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
end

%Detector tilt - should be an identity matrix for ECP, but could be
%specified if the user wants
if isfield(Input_Data,'Detector_tilt')
else
    Input_Data.Detector_tilt=eye(3);
end

if isfield(ECP_Pat,'pattern')
    ECP_pattern=ECP_Pat.ECP_Pat_BG;
    ECP_loaded=1;
end

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','on','Color','w');
f.OuterPosition = [200, 200, 800, 600];
f.Name = 'AstroECP v1';
% f.Position = [150,150,800,600];
f.Units = 'normalized';

%% figure box positioning
xsep=80; xstart=50; xwid=70; ystart=100; ysep=-18; yhig=15;

h_close=uicontrol('style','pushbutton','string','Close',...
    'position',[xstart+xsep*8 ystart+5*ysep xwid yhig]);
h_close.Units='normalized';
h_close.Callback=@closefun;



%% Miller

if isfield(Input_Data,'miller1')
    miller1_text=[num2str(Input_Data.miller1(1)) ',' num2str(Input_Data.miller1(2)) ',' num2str(Input_Data.miller1(3)) ];
    miller1_textui=uicontrol('style','text','string',miller1_text,'Units','normalized','Position',[0.2 0.85 0.1 0.02],'BackgroundColor','w'); %#ok<*NASGU>
end

if isfield(Input_Data,'miller2')
    miller2_text=[num2str(Input_Data.miller2(1)) ',' num2str(Input_Data.miller2(2)) ',' num2str(Input_Data.miller2(3)) ];
    miller2_textui=uicontrol('style','text','string',miller2_text,'Units','normalized','Position',[0.3 0.85 0.1 0.02],'BackgroundColor','w');

end
%% Phase list

phases_list=phasefolder_read(Input_Data);
h_phases_text=uicontrol('style','text','String','Phase file selection' ,'BackgroundColor','w',...
    'position',[xstart+3*xsep ystart+4.4*ysep xwid*1.5 yhig],'Callback',@Phase_Update);
h_phases_text.Units='normalized';

h_phases=uicontrol('style','popupmenu','String',phases_list,...
    'position',[xstart+3*xsep ystart+5*ysep xwid*1.5 yhig],'Callback',@Phase_Update);
h_phases.Units='normalized';

try
    phase_val=find(logical(1-cellfun('isempty',strfind(phases_list,Input_Data.Phase_Input{1}))) == true);
    if isempty(phase_val)
        phase_val=1;
    end
catch
    warning('Phase not found, trying to load silicon 20 kV patterns')
    phase_val=find(logical(1-cellfun('isempty',strfind(phases_list,'Si_20kV'))) == true); %#ok<STRCL1>
    % phase_val=1;
end

h_phases.Value=phase_val(1);

[Crystal_UCell,Crystal_Family,screen_int,Family_List,RTM_info]=phase_data(Input_Data); %#ok<*SETNU>

%% stage rotations (useful for ECP analysis)

h_s_x=uicontrol('style','text','string', ...
    'Stage X-tilt',...
    'position',[xstart+xsep*6 ystart-1.6*ysep xwid yhig],'BackgroundColor','w');
h_s_x.Units = 'normalized';

h_s_y=uicontrol('style','text','string', ...
    'Stage Y-tilt',...
    'position',[xstart+xsep*7 ystart-1.6*ysep xwid yhig],'BackgroundColor','w');
h_s_y.Units = 'normalized';

h_s_z=uicontrol('style','text','string', ...
    'Stage Z-rotation',...
    'position',[xstart+xsep*8 ystart-1.6*ysep xwid yhig],'BackgroundColor','w');
h_s_z.Units = 'normalized';


h_s_xe=uicontrol('style','edit','string', ...
    num2str(round(S_start(1),3)),...
    'position',[xstart+xsep*6 ystart-ysep xwid yhig],'Callback',@Update_PC);
h_s_xe.Units = 'normalized';
h_s_xe.UserData=PC_start;

h_s_ye=uicontrol('style','edit','string', ...
    num2str(round(S_start(2),3)),...
    'position',[xstart+xsep*7 ystart-ysep xwid yhig],'Callback',@Update_PC);
h_s_ye.Units = 'normalized';

h_s_ze=uicontrol('style','edit','string', ...
    num2str(round(S_start(3),3)),...
    'position',[xstart+xsep*8 ystart-ysep xwid yhig],'Callback',@Update_PC);
h_s_ze.Units = 'normalized';

%% Pattern centre boxes
h_pc_x=uicontrol('style','text','string', ...
    'PCx',...
    'position',[xstart+xsep*6 ystart+0.4*ysep xwid yhig],'BackgroundColor','w');
h_pc_x.Units = 'normalized';

h_pc_y=uicontrol('style','text','string', ...
    'PCy',...
    'position',[xstart+xsep*7 ystart+0.4*ysep xwid yhig],'BackgroundColor','w');
h_pc_y.Units = 'normalized';

h_pc_z=uicontrol('style','text','string', ...
    'DD',...
    'position',[xstart+xsep*8 ystart+0.4*ysep xwid yhig],'BackgroundColor','w');
h_pc_z.Units = 'normalized';

h_pc_xe=uicontrol('style','edit','string', ...
    num2str(round(PC_start(1),3)),...
    'position',[xstart+xsep*6 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_xe.Units = 'normalized';
h_pc_xe.UserData=PC_start;

h_pc_ye=uicontrol('style','edit','string', ...
    num2str(round(PC_start(2),3)),...
    'position',[xstart+xsep*7 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_ye.Units = 'normalized';

h_pc_ze=uicontrol('style','edit','string', ...
    num2str(round(PC_start(3),3)),...
    'position',[xstart+xsep*8 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_ze.Units = 'normalized';

h_dpc_x=uicontrol('style','text','string', ...
    'phi_1',...
    'position',[xstart+xsep*6 ystart+2.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_dpc_x.Units = 'normalized';

h_dpc_y=uicontrol('style','text','string', ...
    'Phi',...
    'position',[xstart+xsep*7 ystart+2.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_dpc_y.Units = 'normalized';

h_dpc_z=uicontrol('style','text','string', ...
    'phi_2',...
    'position',[xstart+xsep*8 ystart+2.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_dpc_z.Units = 'normalized';

h_dpc_xe=uicontrol('style','edit','string', ...
    num2str(round(eangs(1),3)),...
    'position',[xstart+xsep*6 ystart+3*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_dpc_xe.Units = 'normalized';

h_dpc_ye=uicontrol('style','edit','string', ...
    num2str(round(eangs(2),3)),...
    'position',[xstart+xsep*7 ystart+3*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_ye.Units = 'normalized';

h_dpc_ze=uicontrol('style','edit','string', ...
    num2str(round(eangs(3),3)),...
    'position',[xstart+xsep*8 ystart+3*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_ze.Units = 'normalized';

e_ang_equiv_phi1=uicontrol('style','text','string', ...
    'NaN',...
    'position',[xstart+xsep*6 ystart+4*ysep xwid yhig],'BackgroundColor','w');
e_ang_equiv_phi1.Units = 'normalized';

e_ang_equiv_Phi=uicontrol('style','text','string', ...
    'NaN',...
    'position',[xstart+xsep*7 ystart+4*ysep xwid yhig],'BackgroundColor','w');
e_ang_equiv_Phi.Units = 'normalized';

e_ang_equiv_phi2=uicontrol('style','text','string', ...
    'NaN',...
    'position',[xstart+xsep*8 ystart+4*ysep xwid yhig],'BackgroundColor','w');
e_ang_equiv_phi2.Units = 'normalized';

h_update=uicontrol('style','pushbutton','string','Update',...
    'position',[xstart+xsep*6 ystart+5*ysep xwid yhig]);
h_update.Units='normalized';
h_update.Callback=@draw_pattern;
h_update.UserData=eangs;


h_exp_bands=uicontrol('style','checkbox','string','Exp Bands',...
    'position',[xstart+xsep*5.5 ystart+1*ysep xwid/2 yhig],'BackgroundColor','w');
h_exp_bands.Units='normalized';
h_exp_bands.Callback=@bands_check;
h_exp_bands.Value=1;

h_ecp_bands=uicontrol('style','checkbox','string','Sim Bands',...
    'position',[xstart+xsep*5 ystart+1*ysep xwid/2 yhig],'BackgroundColor','w');
h_ecp_bands.Units='normalized';
h_ecp_bands.Callback=@bands_check;
h_ecp_bands.Value=1;

h_exp_PC =uicontrol('style','checkbox','string','Exp Center',...
    'position',[xstart+xsep*5.5 ystart+2*ysep xwid/2 yhig],'BackgroundColor','w');
h_exp_PC.Units='normalized';
h_exp_PC.Callback=@bands_check;
h_exp_PC.Value=1;

h_ecp_PC=uicontrol('style','checkbox','string','Sim Center',...
    'position',[xstart+xsep*5 ystart+2*ysep xwid/2 yhig],'BackgroundColor','w');
h_ecp_PC.Units='normalized';
h_ecp_PC.Callback=@bands_check;
h_ecp_PC.Value=1;





h_pc_EAngEq=uicontrol('style','pushbutton','string','Update EA',...
    'position',[xstart+xsep*7 ystart+5*ysep xwid yhig]);
h_pc_EAngEq.Units='normalized';
h_pc_EAngEq.Callback=@Update_EA;

%% Controls for rotating the frame

xstart=50; xsep=80;  xwid=70;  yhig=15; ystart=140;  ysep=-20;

%labels
h_drx=uicontrol('style','text','string', ...
    'dRx',...
    'position',[xstart+xsep*0 ystart+0.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_drx.Units = 'normalized';

h_dry=uicontrol('style','text','string', ...
    'dRy',...
    'position',[xstart+xsep*1 ystart+0.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_dry.Units = 'normalized';

h_drz=uicontrol('style','text','string', ...
    'dRz',...
    'position',[xstart+xsep*2 ystart+0.4*ysep xwid yhig],'Callback',@Update_PC,'BackgroundColor','w');
h_drz.Units = 'normalized';

% text - abs value
h_drabs=uicontrol('style','text','string','Value',...
    'position',[xstart-xsep*0.7 ystart+1*ysep xwid yhig],'BackgroundColor','w');
h_drabs.Units='normalized';

%delta value
h_drdel=uicontrol('style','text','string','Delta',...
    'position',[xstart-xsep*0.7 ystart+2*ysep xwid yhig],'BackgroundColor','w');
h_drdel.Units='normalized';

% Rx

% text - abs value
h_drx_text=uicontrol('style','edit','string','0',...
    'position',[xstart+xsep*0 ystart+1*ysep xwid yhig]);
h_drx_text.Units='normalized';

%delta value
h_drx_dtext=uicontrol('style','edit','string',num2str(DeltaRs(1)),...
    'position',[xstart+xsep*0 ystart+2*ysep xwid yhig]);
h_drx_dtext.Units='normalized';

%buttons to modify
h_drx_minus=uicontrol('style','pushbutton','string','-    ↑','FontSize', 12,...
    'position',[xstart+xsep*0 ystart+3*ysep xwid yhig],'Callback',@rXm);
h_drx_minus.Units='normalized';

h_drx_plus=uicontrol('style','pushbutton','string','+    ↓','FontSize', 12,...
    'position',[xstart+xsep*0 ystart+4*ysep xwid yhig],'Callback',@rXp);
h_drx_plus.Units='normalized';

% Ry

% text - abs value
h_dry_text=uicontrol('style','edit','string','0',...
    'position',[xstart+xsep*1 ystart+1*ysep xwid yhig]);
h_dry_text.Units='normalized';

%delta value
h_dry_dtext=uicontrol('style','edit','string',num2str(DeltaRs(2)),...
    'position',[xstart+xsep*1 ystart+2*ysep xwid yhig]);
h_dry_dtext.Units='normalized';

%buttons to modify
h_dry_minus=uicontrol('style','pushbutton','string','-    →','FontSize', 12,...
    'position',[xstart+xsep*1 ystart+3*ysep xwid yhig],'Callback',@rYm);
h_dry_minus.Units='normalized';

h_dry_plus=uicontrol('style','pushbutton','string','+    ←','FontSize', 12,...
    'position',[xstart+xsep*1 ystart+4*ysep xwid yhig],'Callback',@rYp);
h_dry_plus.Units='normalized';

% Rz

% text - abs value
h_drz_text=uicontrol('style','edit','string','0',...
    'position',[xstart+xsep*2 ystart+1*ysep xwid yhig]);
h_drz_text.Units='normalized';

%delta value
h_drz_dtext=uicontrol('style','edit','string',num2str(DeltaRs(3)),...
    'position',[xstart+xsep*2 ystart+2*ysep xwid yhig]);
h_drz_dtext.Units='normalized';

%buttons to modify
h_drz_minus=uicontrol('style','pushbutton','string','-     ↻','FontSize', 12,...
    'position',[xstart+xsep*2 ystart+3*ysep xwid yhig], 'Callback',@rZm);
h_drz_minus.Units='normalized';

h_drz_plus=uicontrol('style','pushbutton','string','+     ↺','FontSize', 12,...
    'position',[xstart+xsep*2 ystart+4*ysep xwid yhig],'Callback',@rZp);
h_drz_plus.Units='normalized';


%% Band Labels

h_band_labels = gobjects(10,1);  % or in your handles struct, etc.

for ii = 1:10

    h_band_labels(ii) = uicontrol('Style', 'text', ...
        'Position', [xstart-(xsep/4), ystart + yhig*(12 - ii), xwid/2, yhig], ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'String', '');  % empty initially
    h_band_labels(ii).Units = 'normalized';
end

%% plot the Experimental EBSP/ECP
xsep=80; xstart=50; xwid=70; ystart=100; ysep=-18; yhig=15;

h_pc_ye.UserData=Input_Data;
h_pc_ze.UserData=screen_int;
h_drz_dtext.UserData=ECP_loaded;

if ECP_loaded ==1
    h_drz_text.UserData=ECP_pattern;

    h_clim=uicontrol('style','checkbox','string', ...
        'CLim',...
        'position',[xstart+xsep*8.3 ystart-26.5*ysep xwid yhig],'BackgroundColor','w');
    h_clim.Units = 'normalized';
    h_clim.Value=0;

    h_climp=uicontrol('style','edit','string', ...
        '',...
        'position',[xstart+xsep*8.3 ystart-25.5*ysep xwid/2 yhig]);
    h_climp.Units = 'normalized';

    h_climm=uicontrol('style','edit','string', ...
        '',...
        'position',[xstart+xsep*8.3 ystart-24.5*ysep xwid/2 yhig]);
    h_climm.Units = 'normalized';
end


%use these UserData fields to store plot info
h_dry_plus.UserData=[];
h_drz_minus.UserData=[];
h_drz_plus.UserData=[];
h_drx_plus.UserData=[]; %crystal plot
h_s_xe.UserData=[]; %pf plot1
h_s_ye.UserData=[]; %pf plot2

if isfield(Input_Data,'ECP_Pat_clim')
    % clim(Input_Data.ECP_Pat_clim)
    h_clim.Value=1;
    h_climm.String=num2str(Input_Data.ECP_Pat_clim(1));
    h_climp.String=num2str(Input_Data.ECP_Pat_clim(2));
end

%angle subtended
h_alpha = uicontrol('style', 'text', ...
    'String', '', ...
    'Position', [xstart+xsep*2.5 ystart-ysep*4.5 xwid yhig], ...
    'FontSize', 10, ...
    'BackgroundColor', 'w');
h_alpha.Units = 'normalized';

%scale bar toggle switch
h_toggle_scalebar = uicontrol('Style', 'togglebutton', ...
    'String', 'Scalebar', ...
    'Position', [xstart+xsep*5 ystart-ysep xwid yhig], ...
    'Callback', @draw_pattern);
h_toggle_scalebar.Units = 'normalized';

%draw the simulation
draw_pattern;

% label exp and simulated patterns
h_label_expECP = uicontrol('style', 'text', 'string', 'Experimental SA-ECP', ...
    'position', [xstart+xsep*5.5 ystart-ysep*4.5 xwid*1.5 yhig], ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
h_label_expECP.Units = 'normalized';


h_label_simECP = uicontrol('style', 'text', 'string', 'Simulated SA-ECP', ...
    'position', [xstart+xsep*0.5 ystart-ysep*4.5 xwid*1.5 yhig], ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
h_label_simECP.Units = 'normalized';

%fix some label positions
h_exp_old=get(h_drz_minus,'UserData');
h_label_expECP.Position=[h_exp_old.Position(1) h_exp_old.Position(2)-0.05 h_label_expECP.Position(3) h_label_expECP.Position(4)]; 

h_ecp_old=get(h_dry_plus,'UserData');
h_label_simECP.Position=[h_ecp_old.Position(1) h_ecp_old.Position(2)-0.05 h_label_simECP.Position(3) h_label_simECP.Position(4)]; 

h_alpha.Position=[h_alpha.Position(1) h_ecp_old.Position(2)-0.05 h_alpha.Position(3) h_alpha.Position(4)]; 

%% Make the GUI visible.
movegui(f,'onscreen')
drawnow;
f.Visible = 'on';

% pause(1); % delay time to have things in place when maximizes

f.WindowState = 'maximized';
f.Units = 'normalized';
% drawnow;
BandLabelsLegend();



%%

%%% START %%%%% PCA by LUKAS June2024 %%%%%%
% modified TBB 2025-08-01

h_refine=uicontrol('style','pushbutton','string','Refine',...
    'position',[xstart+xsep*5 ystart+5*ysep xwid yhig]);

h_refine.Units='normalized';
h_refine.Callback=@refine;

    function refine(~,eventdata) %#ok<*INUSD>
        psize=size(ECP_Pat.pattern);
        PatternInfo.ScreenWidth=psize(1);
        PatternInfo.ScreenHeight=psize(2);
        [PatternIn,SettingsCor]=EBSP_BGCor(ECP_Pat.pattern,ECP_Pat.Settings_Cor_refine);
        eangs_in=[str2double(e_ang_equiv_phi1.String),str2double(e_ang_equiv_Phi.String),str2double(e_ang_equiv_phi2.String)];
        %%
        pc_in=[str2double(h_pc_xe.String),str2double(h_pc_xe.String),str2double(h_pc_ye.String)];
        
        %% Start the refinement method
        old_color=h_refine.BackgroundColor;
        h_refine.BackgroundColor=[1,0,0];
        h_refine.Enable='off';
        drawnow();
        
        disp('Refinement process started - please wait');
        fprintf('\n');
        if isfield(Input_Data,'Mode')
            [PC_refined,eangs_refined,peakhight] =optimize_pc_ori(PC_start,eangs_in*degree,PatternIn,SettingsCor,screen_int,'Mode',Input_Data.Mode);
        else
            [PC_refined,eangs_refined,peakhight] =optimize_pc_ori(PC_start,eangs_in*degree,PatternIn,SettingsCor,screen_int);
        end
        disp('Refinement process ended');
        fprintf('\n');
        h_refine.BackgroundColor=old_color;
        h_refine.Enable='on';
        plot_refinementHQ(PatternIn,SettingsCor,PatternInfo,screen_int,PC_start,PC_refined,eangs_in*degree,eangs_refined, Input_Data)

        %% update pc
        h_pc_xe.String=PC_refined(1);
        h_pc_ye.String=PC_refined(2);
        h_pc_ze.String=PC_refined(3);
        Update_PC();
        %% update eangs
        e_ang_equiv_phi1.String=eangs_refined(1)/degree;
        e_ang_equiv_Phi.String=eangs_refined(2)/degree;
        e_ang_equiv_phi2.String=eangs_refined(3)/degree;
        % h_update.String=eangs_refined;
        % Update_eangs();
        % Update_EA();
        fprintf('Your refinement had a normalised cross correlation effect of: %f', peakhight);
        fprintf('\n');
    end


%%% END %%%%% PCA by LUKAS %%%%%% lukas



%% sub functions
    function bands_check(~,eventdata)
        draw_pattern;
        % h_exp_bands.Callback=@bands_check;
        %         set(h_exp_bands,'value',1-get(h_exp_bands,'value'))
    end


    function Update_EA(~,eventdata)
        % e_ang_equiv_phi1.String=num2str(e_ang_equiv_deg(1));
        % e_ang_equiv_Phi.String=num2str(e_ang_equiv_deg(2));
        % e_ang_equiv_phi2.String=num2str(e_ang_equiv_deg(3));

        eangs(1)=round(str2double(get(e_ang_equiv_phi1,'String')),3);
        eangs(2)=round(str2double(get(e_ang_equiv_Phi,'String')),3);
        eangs(3)=round(str2double(get(e_ang_equiv_phi2,'String')),3);

        h_dpc_xe.String=e_ang_equiv_phi1.String;
        h_dpc_ye.String=e_ang_equiv_Phi.String;
        h_dpc_ze.String=e_ang_equiv_phi2.String;

        h_s_xe.String='0';
        h_s_ye.String='0';
        h_s_ze.String='0';

        h_drx_text.String='0';
        h_dry_text.String='0';
        h_drz_text.String='0';

        h_drx_text.String='0';
        h_dry_text.String='0';
        h_drz_text.String='0';

        % eangs(1)=round(str2double(get(h_dpc_xe,'string')),3);
        % eangs(2)=round(str2double(get(h_dpc_ye,'string')),3);
        % eangs(3)=round(str2double(get(h_dpc_ze,'string')),3);

        h_update.UserData=eangs;
        Update_eangs

    end

    function draw_pattern(~,eventdata)
        %fetch the old pattern frame, so it can be deleted later


        deactivate_buttons;

        h_ecp_old=get(h_dry_plus,'UserData');

        % h_raw = axes('Units','Pixels','Position',[50,400,200,185],'Parent',f);
        h_ecp = axes('Units','normalized','Position',[0.0612 0.3317 0.3750 0.5],'Parent',f);
        h_dry_plus.UserData=h_ecp;

        InputData=get(h_pc_ye,'UserData');

        %nudge rotations
        drX=str2num(get(h_drx_text,'String'));
        drY=str2num(get(h_dry_text,'String'));
        drZ=str2num(get(h_drz_text,'String'));

        % use the stage convention - Rz*Rx*Ry - BUT NOTE THAT THE TESCAN STAGE USES -Ry
        g_nudge=InputData.Rz(drZ*pi/180)*InputData.Rx(drX*pi/180)*InputData.Ry(drY*pi/180);

        % changing the convention back to Rz*Rx*Rz
        % g_nudge=Input_Data.Rz(drZ*pi/180)*Input_Data.Rx(drY*pi/180)*Input_Data.Rz(drX*pi/180);


        % stage rotations
        sX=str2num(get(h_s_xe,'String')); %#ok<*ST2NM>
        sY=str2num(get(h_s_ye,'String'));
        sZ=str2num(get(h_s_ze,'String'));

        % use the stage convention - Rz*Rx*Ry - BUT NOTE THAT THE STAG USES -Sy
        %working before
        g_stage=Input_Data.Rz(sZ*pi/180)*Input_Data.Rx(sY*pi/180)*Input_Data.Rz(-sX*pi/180);

        % g_stage=InputData.Rz(sZ*pi/180)*InputData.Rx(sY*pi/180)*InputData.Rz(-sX*pi/180);


        eangs_n=get(h_update,'UserData');
        eangs_n=eangs_n*pi/180;

        gmatrix=Input_Data.Rz(eangs_n(3))*Input_Data.Rx(eangs_n(2))*Input_Data.Rz(eangs_n(1));


        % RotMat
        RotMat=gmatrix*Input_Data.Detector_tilt*g_stage*g_nudge;


        %calculate the equivalent g_matrix
        e_ang_equiv=conv_G_to_EA(RotMat);
        e_ang_equiv_deg=round(e_ang_equiv,3)*180/pi;

        e_ang_equiv_phi1.String=num2str(e_ang_equiv_deg(1));
        e_ang_equiv_Phi.String=num2str(e_ang_equiv_deg(2));
        e_ang_equiv_phi2.String=num2str(e_ang_equiv_deg(3));

        % gmatrix2=InputData.Rz(e_ang_equiv(3))*InputData.Rx(e_ang_equiv(2))*InputData.Rz(e_ang_equiv(1));
        PC_start=get(h_pc_xe,'UserData');

        [EBSD_Geometry ] = EBSP_Gnom( ECP_Pat,PC_start); %you can change PC_in if you want

        sf=sf_calc(EBSD_Geometry);

        if PC_start(3) > 1 %conic sampling values
            cone_params=720;
        else
            cone_params=360;
        end

        % [ Pattern_Sim ] = EBSP_gen( EBSD_Geometry,RotMat,screen_int,screen_int.isHex ); %generate the EBSP for this iteration

        % Correction for matching patterns with Bruker DynamicS convention (March 11, 2024), 
        % this fix is required to deal with Pbnm vs Pnma symmetry issue in DynamicS
        % the issue was observed when using olivine in DynamicS
        % 
        % description tidied up 2025-08-01, TBB

        g_dynamics=Input_Data.Rx(pi/2)*Input_Data.Rz(pi/2);
        [Pattern_Sim]=EBSP_gen( EBSD_Geometry,g_dynamics*RotMat,screen_int);

        i_ecp=imagesc(EBSD_Geometry.x_screen,EBSD_Geometry.y_screen,Pattern_Sim,'Parent',h_ecp);
        colormap('gray');
        % axis tight; axis xy; axis equal;
        axis(h_ecp,'equal','xy','tight');
        axis off;
        hold(h_ecp,'on');


        %plot the band conic sections
        cset=Input_Data.cset;
        if h_ecp_bands.Value == 1

            %reusing the previous method
            % Plot_EBSPAnnotated_TZ( Pattern_Sim,EBSD_Geometry,[],RotMat,Crystal_UCell{1},Crystal_Family{1},h_ecp, Input_Data.V_in);

            if ~isfield(Crystal_UCell{1},'lambda_p')
                Crystal_UCell{1}.lambda_p=0.0859;
            end

            if ~isempty(Family_List)
                for n=1:size(Family_List,1)
                    HKLs=Family_List{n}(:,2:4);
                    
                
                    for p=1:size(HKLs,1)
                        %calculate the bands
                        [bands] = Cone_Build(HKLs(p,:),Crystal_UCell{1}.Astar,Crystal_UCell{1}.lambda_p,RotMat,cone_params);
                        %plot them on the EBSP/ECP
                        pBand(bands,cset(n,:),h_ecp,sf)

                    end
                end
            end
        end
        h_ecp.XLim=([EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end)]);
        h_ecp.YLim=([EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end)]);
        h_ecp.FontSize=10;

        xl = h_ecp.XLim;        %saved for angle scale bar
        yl = h_ecp.YLim;








        %plot the exp ECP if needed
        if h_drz_dtext.UserData == 1

            %plot the ecp
            h_exp_old=get(h_drz_minus,'UserData');
            h_exp = axes('Units','normalized','Position',[0.5613 0.3317 0.3750 0.5],'Parent',f);

            h_drz_minus.UserData=h_exp;

            ECP_data=h_drz_text.UserData;
            InputExp.size=size(ECP_data);
            [EBSD_Geometry_exp ] = EBSP_Gnom( InputExp,PC_start); %you can change PC_in if you want
            i_ecp=imagesc(EBSD_Geometry_exp.x_screen,EBSD_Geometry_exp.y_screen,ECP_data,'Parent',h_exp);

            colormap('gray');
            % axis tight; axis xy; axis equal;

            axis(h_exp,'equal','xy','tight', 'off');
            hold(h_exp,'on');

            if h_clim.Value == 1

                % h_climm.String=get(Input_Data.ECP_Pat_clim(1));
                % h_climp.String=num2str(Input_Data.ECP_Pat_clim(2));

                clims=sort([str2double(h_climm.String),str2double(h_climp.String)]);
                try
                    h_exp.CLim=clims;
                catch
                    disp('The clim values you have included do not work - try again')
                end

            end
            %plot the band conic sections
            cset=Input_Data.cset;

            theta_family=zeros(size(Family_List,1),1);

            if h_exp_bands.Value == 1


                % Plot_EBSPAnnotated_TZ( ECP_data,EBSD_Geometry,[],RotMat,Crystal_UCell{1},Crystal_Family{1},h_exp, Input_Data.V_in);



                if ~isempty(Family_List)
                    for n=1:size(Family_List,1)
                        HKLs=Family_List{n}(:,2:4);

                        for p=1:size(HKLs,1)
                            % %calculate the bands
                            [bands] = Cone_Build(HKLs(p,:),Crystal_UCell{1}.Astar,Crystal_UCell{1}.lambda_p,RotMat,cone_params);
                            % %plot them on the EBSP/ECP
                            pBand(bands,cset(n,:),h_exp,sf);

                        end

                        theta_family(n)=bands.theta;
                    end
                end
            end

            h_exp.XLim=([EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end)]);
            h_exp.YLim=([EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end)]);



            % display angle subtended
            DD = str2double(get(h_pc_ze, 'String'));
            alpha_rad = 2 * atan(1 / (2 * DD));     % radians
            alpha_deg = rad2deg(alpha_rad);
            alpha_str = sprintf('2α = %.2f°', alpha_deg);
            set(h_alpha, 'String', alpha_str, 'BackgroundColor', 'w');

            % display angular scale bar
            if isgraphics(h_ecp)
                xl = get(h_ecp, 'XLim');
                yl = get(h_ecp, 'YLim');
            end
            h_alpha_step = floor(alpha_deg / 5);    %recompute bar length in data units so it exactly spans 'h_alpha_step'
            fullWidth = xl(2) - xl(1);
            barLen_data = (h_alpha_step / alpha_deg) * fullWidth;
            yOffset = 0.05 * (yl(2) - yl(1));

            dataMarginX = 0.02 * fullWidth;
            x_start = xl(2) - dataMarginX - barLen_data;
            x_end   = xl(2) - dataMarginX;
            y_line  = yl(1) + yOffset;

            % line([x_start, x_end], [y_line, y_line], ...
            %     'Color','k', 'LineWidth', 2.5, 'Parent', h_ecp);
            %
            % text((x_start + x_end)/2, y_line - 0.01*(yl(2) - yl(1)), ...
            %      sprintf('%d°', h_alpha_step), ...
            %      'Parent', h_ecp, ...
            %      'HorizontalAlignment', 'center', ...
            %      'VerticalAlignment', 'top', ...
            %      'BackgroundColor', 'none', 'FontSize', 13, 'FontWeight', 'bold');


            if get(h_toggle_scalebar, 'Value') == 1
                set(h_toggle_scalebar, 'String', 'Show scalebar on Sim pattern');
                line([x_start, x_end], [y_line, y_line], ...
                    'Color','w', 'LineWidth', 3, 'Parent', h_exp);


                text((x_start + x_end)/2, y_line - 0.005*(yl(2) - yl(1)), ...
                    sprintf('%d°', h_alpha_step), ...
                    'Parent', h_exp, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'top', ...
                    'BackgroundColor', 'none', 'FontSize', 14, 'FontWeight', 'bold', 'Color','w');
            else

                set(h_toggle_scalebar, 'String', 'Show scalebar on Exp pattern');
                line([x_start, x_end], [y_line, y_line], ...
                    'Color','k', 'LineWidth', 3, 'Parent', h_ecp);


                text((x_start + x_end)/2, y_line - 0.005*(yl(2) - yl(1)), ...
                    sprintf('%d°', h_alpha_step), ...
                    'Parent', h_ecp, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'top', ...
                    'BackgroundColor', 'none', 'FontSize', 14, 'FontWeight', 'bold', 'Color','k');


            end





            if h_exp_PC.Value == 1
                %pattern center EXP
                f1(2) = scatter(0, 0, 150, 'wo', 'LineWidth', 1.5, 'MarkerEdgeColor', 'k','Parent', h_exp);
                f1(3) = scatter(0, 0, 150, 'kx', 'Parent', h_exp);
            end

            if h_ecp_PC.Value == 1
                %pattern center EXP
                f1(2) = scatter(0, 0, 150, 'wo', 'LineWidth', 1.5,'MarkerEdgeColor', 'k','Parent', h_ecp);
                f1(3) = scatter(0, 0, 150, 'kx', 'Parent', h_ecp);
            end


            %store the theta list somewhere
            h_drabs.UserData=theta_family;

            %plot the histogram
            h_exp2_old=get(h_drz_plus,'UserData');
            %longer
            % h_exp2 = axes('Units','normalized','Position',[0.5613 0.8717 0.3750 0.1],'Parent',f);
            %shorter
            h_exp2 = axes('Units','normalized','Position',[0.62 0.9 0.25 0.07],'Parent',f);
            h_drz_plus.UserData=h_exp2;
            histogram(ECP_data(:),100,'Parent',h_exp2);
            xlabel('Intensity, a.u', FontSize=6);
            ylabel('Frequency, a.u.', FontSize=6);

            xlim(h_exp2, [min(ECP_data(:)), max(ECP_data(:))]);  % Set x-axis limits

            ylim(h_exp2, [0 50000]);
            % ylim auto     % Automatically adjust y-axis limits
            h_exp2.FontSize=9;
            h_exp.FontSize=9;

            delete(h_exp_old);
            delete(h_exp2_old);
            % delete(h_hist_old);
        end
        delete(h_ecp_old);


        % plot the crystal shape unit cell

        if isfield(Input_Data,'crystal_shape')

            h_cry_old=get(h_drx_plus,'UserData');
            % bigger shape size
            % h_cry = axes('Units','normalized','Position',[0.37 0.73 0.25 0.25],'Parent',f);
            % normal shape size



            h_pf1_old=get(h_s_xe,'UserData');
            h_pf2_old=get(h_s_ye,'UserData');

            h_pf1 = axes('Units','normalized','Position',[0.2 0.8817 0.1 0.1],'Parent',f,'Visible','off');
            hold on;
            h_pf2 = axes('Units','normalized','Position',[0.3 0.8817 0.1 0.1],'Parent',f,'Visible','off');
            hold on;

            % axis off; axis equal;
            h_s_xe.UserData=h_pf1;
            h_s_xe.UserData=h_pf2;

            try
                cs_phase=loadCIF(RTM_info.cif_file);
                assignin('base','cs_phase',cs_phase)%% better in the future
                prop.x=0;
                prop.y=0;
                ori_single = rotation('Matrix',RotMat');

                if strcmpi(Input_Data.crystal_shape,'cube')
                    cS=crystalShape.cube(cs_phase);
                elseif  strcmpi(Input_Data.crystal_shape,'hex')
                    cS=crystalShape.hex(cs_phase);

                elseif  strcmpi(Input_Data.crystal_shape,'orthorhombic')
                    cs_phase = crystalSymmetry.load(RTM_info.cif_file); %load cif file


                    mm = Miller({1,0,0},cs_phase);
                    rr = Miller({0,1,0},cs_phase);
                    zz = Miller({0,0,1},cs_phase);
                    NN = [rr,zz,mm];
                    cS = crystalShape(NN);
                    % cS=crystalShape.forsterite;
                end
            catch
                warning('Failed to load the CIF properly for the crystal shape information')
            end

            %plot the unit crystal
            try
                h_cry = axes('Units','normalized','Position',[0.1 0.8817 0.1 0.1],'Parent',f);
                h_drx_plus.UserData=h_cry;
                axis on; axis equal; hold on;

                plot(ori_single * cS,'LineStyle', '-', 'FaceColor','none','Parent',h_cry);
                plot(ori_single * cS(Miller(0,0,1,cs_phase)),'FaceColor','LightSkyBlue','FaceAlpha',0.8,'lineWidth', 1 ,'Parent',h_cry, axis, 'off')
                plot(ori_single * cS(Miller(0,0,-1,cs_phase)),'FaceColor','LightSkyBlue','FaceAlpha',0.8,'lineWidth', 1 ,'Parent',h_cry)
                plot(ori_single * cS(Miller(1,0,0,cs_phase)),'FaceColor','DarkSalmon','FaceAlpha',0.8,'lineWidth', 1  ,'Parent',h_cry)
                plot(ori_single * cS(Miller(-1,0,0,cs_phase)),'FaceColor','DarkSalmon','FaceAlpha',0.8,'lineWidth', 1  ,'Parent',h_cry)
                plot(ori_single * cS(Miller(0,1,0,cs_phase)),'FaceColor','LightSeaGreen','FaceAlpha',0.8,'lineWidth', 1  ,'Parent',h_cry)
                plot(ori_single * cS(Miller(0,-1,0,cs_phase)),'FaceColor','LightSeaGreen','FaceAlpha',0.8,'lineWidth', 1  ,'Parent',h_cry)

            catch
                warning('Unit Cell Plot failed - you can copy "AstroEBSD\modules\AstroECP\plotCS_Update.m" to the folder "\mtex-5.11.1\geometry\@crystalShape" and replace plot.m with this alternative version')
                
                fileLocation=which('plotCS_Update.m');
                disp('Copying file automatically')
                copyfile([Input_Data.mtex_location 'geometry\@crystalShape\plot.m'],[Input_Data.mtex_location 'geometry\@crystalShape\plot_old.m']);
                copyfile(fileLocation,[Input_Data.mtex_location 'geometry\@crystalShape\plot.m']);
                

            end

            try
                if isfield(Input_Data,'miller1')

                    ebsd_single = EBSD(ori_single, 1,{'notIndexed',cs_phase},prop);
                    plotPDF(ebsd_single.orientations,[Miller(Input_Data.miller1(1),Input_Data.miller1(2),Input_Data.miller1(3),cs_phase)],'NoTitle','MarkerSize',3,'parent',h_pf1,'Color','w');
                    axis off
                else
                    h_pf1.Visible='off';
                end

                if isfield(Input_Data,'miller2')

                    ebsd_single = EBSD(ori_single, 1,{'notIndexed',cs_phase},prop);
                    plotPDF(ebsd_single.orientations,[Miller(Input_Data.miller2(1),Input_Data.miller2(2),Input_Data.miller2(3),cs_phase)],'NoTitle','MarkerSize',3,'parent',h_pf2,'Color','w');
                    axis off
                else
                    h_pf2.Visible='off';
                end

            catch
                warning('MTEX ODF plots not working for this')
            end
            delete(h_pf1_old);
            delete(h_pf2_old);
            delete(h_cry_old);
        end




        drawnow();
        guidata(f);
        activate_buttons;




        % % display angle subtended
        %  DD = str2num(get(h_pc_ze, 'String'));
        %     alpha_rad = 2*atan(1 / (2 * DD));     % radians
        %     alpha_deg = rad2deg(alpha_rad);
        %      alpha_str = sprintf('2α = %.2f°', alpha_deg);
        % set(h_alpha, 'String', alpha_str, 'BackgroundColor', 'w');

        % % Draw angular scale bar
        %
        % % --- pick a “nice” step (e.g. 1/5 of total, rounded down to whole degrees) ---
        % h_alpha_step = floor(alpha_deg / 5);     % e.g. 12.63 → 12°
        %
        % % --- recompute bar length in data units so it exactly spans 'h_alpha_step' ---
        % % xl = h_ecp.XLim;  %defined earlier
        % fullWidth = xl(2) - xl(1);
        % barLen_data = (h_alpha_step / alpha_deg) * fullWidth;
        %
        % % --- vertical offset (5% up from bottom) ---
        % % yl = h_ecp.YLim;  %defined earlier
        % yOffset = 0.05 * (yl(2) - yl(1));
        %
        % % --- define start/end X in data coords, with small data‑margin ---
        % dataMarginX = 0.02 * fullWidth;
        % x_start = xl(2) - dataMarginX - barLen_data;
        % x_end   = xl(2) - dataMarginX;
        % y_line  = yl(1) + yOffset;
        %
        % % --- draw the bar ---
        %  line([x_start, x_end], [y_line, y_line], ...
        %      'Color','k', 'LineWidth', 2.5, 'Parent', h_ecp_old);
        %
        % % --- label it with the rounded-down value --\
        % text((x_start + x_end)/2, y_line - 0.01*(yl(2) - yl(1)), ...
        %      sprintf('%d°', h_alpha_step), ...
        %      'Parent', h_ecp_old, ...
        %      'HorizontalAlignment', 'center', ...
        %      'VerticalAlignment', 'top', ...
        %      'BackgroundColor', 'none', 'FontSize', 12, 'FontWeight', 'bold');







    end

    function deactivate_buttons
        h_ring.Enable='off'; %#ok<STRNU>
        h_update.Enable='off';
        h_drx_minus.Enable='off';
        h_dry_minus.Enable='off';
        h_drz_minus.Enable='off';

        h_drx_plus.Enable='off';
        h_dry_plus.Enable='off';
        h_drz_plus.Enable='off';

    end

    function activate_buttons
        h_ring.Enable='on'; %#ok<STRNU>
        h_update.Enable='on';
        h_drx_minus.Enable='on';
        h_dry_minus.Enable='on';
        h_drz_minus.Enable='on';

        h_drx_plus.Enable='on';
        h_dry_plus.Enable='on';
        h_drz_plus.Enable='on';

    end

    function pBand(bands,csetn,h_ecp,sf)
        % bands_clean_centre=cBand(bands.centre_zp,sf);
        % if ~isempty(bands_clean_centre.x)
        %     plot(bands_clean_centre.x,bands_clean_centre.y,'color',csetn,'LineWidth',1,'LineStyle','--','Parent',h_ecp)
        % end


        bands_clean_upper=cBand(bands.upper_zp,sf);
        if ~isempty(bands_clean_upper.x) && ~isempty(bands_clean_upper.y)
            plot(bands_clean_upper.x,bands_clean_upper.y,'color',csetn,'LineWidth',1,'LineStyle','-','Parent',h_ecp)  % line width Line size
            % plot(bands_clean_upper.x,'color',csetn,'LineWidth',0.5,'LineStyle','-','Parent',h_ecp)  % line width Line size
        end


        bands_clean_lower=cBand(bands.lower_zp,sf);
        if ~isempty(bands_clean_lower.x) && ~isempty(bands_clean_lower.y)
            plot(bands_clean_lower.x,bands_clean_lower.y,'color',csetn,'LineWidth',1,'LineStyle','-','Parent',h_ecp)
        end



    end

    function bands_clean=cBand(bzp,sf)



        % % Problem of duplicate bands on one side
        % bands_clean.x=bzp(1,:)./bzp(3,:);
        % bands_clean.y=bzp(2,:)./bzp(3,:);
        % bands_clean.y=bands_clean.y(bands_clean.x>=sf.xmin & bands_clean.x<=sf.xmax);
        % bands_clean.x=bands_clean.x(bands_clean.x>=sf.xmin & bands_clean.x<=sf.xmax);
        %
        % bands_clean.x=bands_clean.x(bands_clean.y>=sf.ymin & bands_clean.y<=sf.ymax);
        % bands_clean.y=bands_clean.y(bands_clean.y>=sf.ymin & bands_clean.y<=sf.ymax);


        x = bzp(1,:) ./ bzp(3,:);
        y = bzp(2,:) ./ bzp(3,:);

        % Create joint mask to filter only valid (x,y) pairs
        mask = x >= sf.xmin & x <= sf.xmax & ...
            y >= sf.ymin & y <= sf.ymax;

        % Apply the mask once, preserving pairing
        bands_clean.x = x(mask);
        bands_clean.y = y(mask);



    end




    function Update_PC(source,eventdata)
        Update_eangs

        PC_start(1)=round(str2double(get(h_pc_xe,'string')),3);
        PC_start(2)=round(str2double(get(h_pc_ye,'string')),3);
        PC_start(3)=round(str2double(get(h_pc_ze,'string')),3);
        h_pc_xe.UserData=PC_start;

        guidata(f);
    end

    function rXp(source,eventdata)
        rvals(1)=round(str2double(get(h_drx_text,'string')),3);
        drvals(1)=round(str2double(get(h_drx_dtext,'string')),3);
        rvals(1)=rvals(1)+drvals(1); %x+
        h_drx_text.String=num2str(rvals(1));
        guidata(f);
        draw_pattern;
    end

    function rXm(source,eventdata)
        rvals(1)=round(str2double(get(h_drx_text,'string')),3);
        drvals(1)=round(str2double(get(h_drx_dtext,'string')),3);
        rvals(1)=rvals(1)-drvals(1); %x+
        h_drx_text.String=num2str(rvals(1));
        guidata(f);
        draw_pattern;
    end


    function rYp(source,eventdata)
        rvals(1)=round(str2double(get(h_dry_text,'string')),3);
        drvals(1)=round(str2double(get(h_dry_dtext,'string')),3);
        rvals(1)=rvals(1)+drvals(1); %x+
        h_dry_text.String=num2str(rvals(1));
        draw_pattern;
        guidata(f);

    end

    function rYm(source,eventdata)
        rvals(1)=round(str2double(get(h_dry_text,'string')),3);
        drvals(1)=round(str2double(get(h_dry_dtext,'string')),3);
        rvals(1)=rvals(1)-drvals(1); %x+
        h_dry_text.String=num2str(rvals(1));
        draw_pattern;
        guidata(f);

    end


    function rZp(source,eventdata)
        rvals(1)=round(str2double(get(h_drz_text,'string')),3);
        drvals(1)=round(str2double(get(h_drz_dtext,'string')),3);
        rvals(1)=rvals(1)+drvals(1); %x+
        h_drz_text.String=num2str(rvals(1));
        draw_pattern;
        guidata(f);
    end

    function rZm(source,eventdata)
        rvals(1)=round(str2double(get(h_drz_text,'string')),3);
        drvals(1)=round(str2double(get(h_drz_dtext,'string')),3);
        rvals(1)=rvals(1)-drvals(1); %x+
        h_drz_text.String=num2str(rvals(1));
        draw_pattern;
        guidata(f);

    end
    function Update_eangs(source,eventdata)
        eangs(1)=round(str2double(get(h_dpc_xe,'string')),3);
        eangs(2)=round(str2double(get(h_dpc_ye,'string')),3);
        eangs(3)=round(str2double(get(h_dpc_ze,'string')),3);
        h_update.UserData=eangs;
        draw_pattern;
        guidata(f);

    end



    function pha_names=phasefolder_read(InputUser)
        if exist(InputUser.Phase_Folder,'dir') ~=7
            error(['The folder ' InputUser.Phase_Folder ' does not exist']);
        end
        phasfolder=fullfile(InputUser.Phase_Folder,'phasefiles'); %bug fix due to phase folder restructure
        folder_dir = dir(phasfolder);
        [~,idx] = sort([folder_dir.datenum]);
        folder_dir=folder_dir(idx);
        folder_files={folder_dir.name};

        %find the pha names
        pha_names=folder_files(logical(1-cellfun('isempty',strfind(folder_files,'pha')))); %#ok<STRCL1>
        for n=1:size(pha_names,2)
            pha_names{n}=pha_names{n}(1:end-4);
        end

    end


    function sf=sf_calc(EBSD_Geometry)
        if EBSD_Geometry.x_gn_min > 0
            sf.xmin=0.9*EBSD_Geometry.x_gn_min;
        else
            sf.xmin=1.1*EBSD_Geometry.x_gn_min;
        end

        if EBSD_Geometry.y_gn_min > 0
            sf.ymin=0.9*EBSD_Geometry.y_gn_min;
        else
            sf.ymin=1.1*EBSD_Geometry.y_gn_min;
        end

        if EBSD_Geometry.x_gn_max > 0
            sf.xmax=1.1*EBSD_Geometry.x_gn_max;
        else
            sf.xmax=0.9*EBSD_Geometry.x_gn_max;
        end

        if EBSD_Geometry.y_gn_max > 0
            sf.ymax=1.1*EBSD_Geometry.y_gn_max;
        else
            sf.ymax=0.9*EBSD_Geometry.y_gn_max;
        end
    end

    function closefun(source,eventdata)
        Update_eangs;

        Output_Data=struct;
        eangs(1)=round(str2double(get(e_ang_equiv_phi1,'String')),3);
        eangs(2)=round(str2double(get(e_ang_equiv_Phi,'String')),3);
        eangs(3)=round(str2double(get(e_ang_equiv_phi2,'String')),3);

        Output_Data.eangs=eangs;
        Output_Data.drx=str2num(get(h_drx_text,'String'));
        Output_Data.dry=str2num(get(h_dry_text,'String'));
        Output_Data.drz=str2num(get(h_drz_text,'String'));

        Output_Data.sX=str2num(get(h_s_xe,'String'));
        Output_Data.sY=str2num(get(h_s_ye,'String'));
        Output_Data.sZ=str2num(get(h_s_ze,'String'));
        Output_Data.dd=str2num(get(h_pc_ze,'String'));

        %store the theta list somewhere
        Output_Data.theta_list=get(h_drabs,'UserData');
        assignin('base', 'Output_Data', Output_Data);

        close(gcf)
    end

    function [Crystal_UCell,Crystal_Family,screen_int,Family_List,RTM_info]=phase_data(Input_Data)
        Phase_Names=get(h_phases,'String');
        Phase_Name=Phase_Names(get(h_phases,'Value'));

        [ Crystal_UCell,Crystal_Family,~,~,~, RTM_info ] = Phase_Builder_RTM(Phase_Name,Input_Data.Phase_Folder);

        [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

        %extract the crystal families
        Family_List=Crystal_Family{1};

    end

    function Phase_Update(source,eventdata)
        [Crystal_UCell,Crystal_Family,screen_int,Family_List,RTM_info]=phase_data(Input_Data);
        draw_pattern;
        guidata(f);
        BandLabelsLegend();

    end

% % % % % % % % % % % % % % % % % % %  Adding new features  % % % % % % % % % % % % % % % % % % % %

%% Create a new button for marking reference
h_mark_reference = uicontrol('style', 'pushbutton', 'string', 'Make as Ref. Position',...
    'Position', [xstart+xsep*0, ystart+5*ysep, 70, yhig], 'Callback', @ref_position);
h_mark_reference.Units='normalized';

    function ref_position(~,~)
        % store values to use later
        x_ref = num2str (get(h_drx_text,'String'));
        y_ref = num2str (get(h_dry_text,'String'));
        z_ref = num2str (get(h_drz_text,'String'));
        setappdata(gcf, 'x_ref', x_ref);
        setappdata(gcf, 'y_ref', y_ref);
        setappdata(gcf, 'z_ref', z_ref);
        h_mark_position.Enable='on';
    end

%% Create a new button for marking new positions
h_mark_position = uicontrol('style', 'pushbutton', 'string', 'Mark New Point',...
    'Position', [xstart+xsep*1, ystart+5*ysep, 70, yhig], 'Callback', @new_position);
h_mark_position.Units='normalized';
h_mark_position.Enable='off';

    function new_position(~,~)
        Update_EA;
        [x,y]=ginput(1);                              % picks XY from user click
        scatter(x , y, 'filled','o','LineWidth',2);   % draws a point where user clicked
        alpha = asind(y/sqrt(x^2 + y^2));             % measures the angle
        x=rad2deg(x);                                 % converts the XY to degrees
        y=rad2deg(y);
        %Feeds the XY coordinates to drx and dry in degrees
        h_drx_text.String = sprintf('%.3f', y);
        h_dry_text.String = sprintf('%.3f', x);

        % store values to use later
        x_nav = num2str (get(h_drx_text,'String'));
        y_nav = num2str (get(h_dry_text,'String'));
        z_nav = num2str (get(h_drz_text,'String'));

        setappdata(gcf, 'x_nav', x_nav);
        setappdata(gcf, 'y_nav', y_nav);
        setappdata(gcf, 'z_nav', z_nav);

        h_navigate.Enable='on';

    end

%% Create button for navigating to new position and display Stage Parameters
h_navigate = uicontrol('style', 'pushbutton', 'string', 'Navigate to New Point',...
    'Position', [xstart+xsep*2, ystart+5*ysep, 70, yhig], 'Callback', @navigate);
h_navigate.Units='normalized';
h_navigate.Enable = 'off';

    function navigate(~,~)
        draw_pattern;
        Update_EA;

        %calling back the ref and new positions
        x_ref = getappdata(gcf, 'x_ref');
        y_ref = getappdata(gcf, 'y_ref');
        z_ref = getappdata(gcf, 'z_ref');
        x_nav = getappdata(gcf, 'x_nav');
        y_nav = getappdata(gcf, 'y_nav');
        z_nav = getappdata(gcf, 'z_nav');

        %Calculations for navigation
        stage_nav_x = str2double(x_nav) + str2double(x_ref);
        stage_nav_y = - (str2double(y_nav) + str2double(y_ref));
        stage_nav_z = str2double(z_nav) + str2double(z_ref);

        %Feeds the navigations to stage controls sx, sy, sz in degrees
        h_s_xe.String = sprintf('%.3f', stage_nav_x);
        h_s_ye.String = sprintf('%.3f', stage_nav_y);
        h_s_ze.String = sprintf('%.3f', stage_nav_z);
        eangs(1) = '0';
        eangs(2) = '0';
        eangs(3) = '0';

    end

%% Create button for loading ECP

% %corner position
% h_load_ecp = uicontrol('style', 'pushbutton', 'string', 'Load New ECP',...
%     'Position', [xstart+xsep*5 ystart-ysep xwid yhig], 'Callback', @load_ecp);
%
%     h_load_ecp.Units='normalized';
%
%     function load_ecp(~,~)
% %% Use uiopen to let the user choose a new file
% [filename, pathname] = uigetfile('*.tif', 'Select an image file');
%
% % Check if the user clicked 'Cancel'
% if isequal(filename,0) || isequal(pathname,0)
%     disp('User canceled the file selection.')
%     return;  % Exit the script if the user cancels the file selection
% end
% close all;
% % Update Input_Data.image_folder and Input_Data.image_name with the new values
% Input_Data.image_folder = pathname;
% Input_Data.image_name = filename;
% Input_Data.ExpImage_Image = fullfile(Input_Data.astro_location, Input_Data.image_folder,Input_Data.image_name);
% Input_Data.ECP_Pat_clim=[2 5]; % default settings of histogram
% Input_Data.PC_in=[0.5 0.5 0.8]; % starting PC - AstroEBSD convention [PCx, PCy, DD]
% Input_Data.Stage_in=[0 0 0]; % stage rotations, in degrees [Rx, Ry, Rz]
% Input_Data.eangs=[0 0 0];
%
% %% Experimental pattern load
% if isfield(Input_Data,'ExpImage_Image')
%     if isfield(Input_Data,'image_frame') %a tescan frame
%         ExpImage_Filename=fullfile(Input_Data.image_folder,Input_Data.image_name);
%         [Input_Data.ExpImage_Image,data1,td] = TescanFrame_Load(ExpImage_Filename,Input_Data.image_frame);
%         Input_Data.size = size(Input_Data.ExpImage_Image); %size of the library patterns and the resize of the raw EBSP
%         Input_Data.ECP_Pat=flipud(double(Input_Data.ExpImage_Image));
%     else
%         %load this frame - normal image loader
%         Input_Data.ECP_Pat=imread(Input_Data.ExpImage_Image);
%         if size(Input_Data.ECP_Pat,3) == 3
%             Input_Data.ECP_Pat=rgb2gray(Input_Data.ECP_Pat);
%         end
%         Input_Data.ECP_Pat=double(flipud(Input_Data.ECP_Pat));
%         Input_Data.size =size(Input_Data.ECP_Pat);
%    end
% end
% %%
% Input_Data.PC_in=[0.5 0.5 3.88]; % starting PC - AstroEBSD convention
% Input_Data.eangs=[0  0  0]; % [phi1, Phi, phi2]
% Input_Data.DeltaRs=[0.1 0.1 2]; % delta values for controlling tilts etc, in degrees
%
% %% run the GUI
% [Output_Data]=f_AstroECP(Input_Data);
%
% end

%% Display info from the loaded ECP hdr file

try
    ftif = Input_Data.image_name;
variables_to_extract = ["OrigFileName", "ScanMode", "AcceleratorVoltage", "WD", "DwellTime", "PixelSizeX", "PredictedBeamCurrent", "StageRotation", "StageTilt", "StageTiltY'" ];

% Call the function header_read_tescan
[data1, table_data] = header_read_tescan(ftif, variables_to_extract); %#ok<*ASGLU>
% [text_retrieved] = fReadTFSHeaderPair(text_to_find,tif_info_string);
if data1 == 0
    try %lets try to load a TFS data set
        info1 = imfinfo(fullfile(Input_Data.image_folder,Input_Data.image_name));
        i1_text=info1.UnknownTags.Value;
        i1_text_pairs=regexp(i1_text, '[\f\n\r]', 'split');
        nonEmptyCells = ~cellfun('isempty', i1_text_pairs);
        i1_text_pairs_full=i1_text_pairs(nonEmptyCells);
        t=1;
        data1=cell(1,10);
        data1{1}=Input_Data.image_name;
        data1{2}=fReadTFSHeaderPair('UseCase',i1_text_pairs_full);
        data1{3}=str2double(fReadTFSHeaderPair('HV',i1_text_pairs_full))/1000; % in kV
        data1{4}=str2double(fReadTFSHeaderPair('WD',i1_text_pairs_full))*1000; % in mm;
        data1{5}=str2double(fReadTFSHeaderPair('BeamCurrent',i1_text_pairs_full))*1E9; % in nA;
        data1{6}=str2double(fReadTFSHeaderPair('FrameTime',i1_text_pairs_full)); %in s
        data1{7}=str2double(fReadTFSHeaderPair('AngularFieldWidth',i1_text_pairs_full))*180/pi;
        data1{8}=str2double(fReadTFSHeaderPair('StageT',i1_text_pairs_full))*180/pi;
        data1{9}=str2double(fReadTFSHeaderPair('StageR',i1_text_pairs_full))*180/pi;
        data1{10}=str2double(fReadTFSHeaderPair('ScanRotation',i1_text_pairs_full))*180/pi;
        
        %do some rounding for the display
        data1{4}=round(data1{4}*10)/10;
        data1{7}=round(data1{7}*100)/100;
        data1{8}=round(data1{8}*10)/10;
        data1{9}=round(data1{9}*10)/10;
        data1{10}=round(data1{10}*10)/10;

        % Display data with proper formatting
        data2 = sprintf('File Name: %s\n\nUse Case: %s\nBeam Energy: %s keV\nWorking Distance: %s mm\nBeam Current: %s nA\nFrameTime: %s S\nAngularField: %s deg\n\nMICROSCOPE SETUP \nStage Tilt: %s°\nStage Rot: %s°\nScan Rot: %s°', ...
        data1{1}, data1{2}, data1{3}, data1{4}, data1{5}, data1{6}, data1{7}, data1{8}, data1{9}, data1{10});

    catch %it messed up, so just populate with unknown
        data2=['Loaded pattern, but no added info'];
    end
else %it was tescan data
    data1{3} = num2str(str2double(data1{3}) / 1000);
    data1{4} = sprintf('%.2f', str2double(data1{4}) * 1e3);
    data1{5} = num2str(str2double(data1{5}) * 100000);
    data1{6} = sprintf('%.1f', str2double(data1{6}) * 1e9);
    data1{7} = sprintf('%.1f', str2double(data1{7}) * 1e9);
    data1{8} = sprintf('%.3f', data1{8});

    % Display data with proper formatting
    data2 = sprintf('File Name: %s\n\nScan Mode: %s\nBeam Energy: %s keV\nWorking Distance: %s mm\nDwell Time: %s µs\nPixel Size: %s nm\nBeam Current: %s nA\n\nMICROSCOPE SETUP \nStage Rotation: %s°\nStage Tilt X: %s°\nStage Tilt Y: %s°', ...
        data1{1}, data1{2}, data1{3}, data1{4}, data1{5}, data1{6}, data1{7}, data1{8}, data1{9}, data1{10});

end


% change the units


catch

end


data2 = strtrim(data2); % Remove newline characters

h_ecp_info=uicontrol('style','text','string', data2,'position',[xstart-(xsep/2) ystart+yhig*15 xwid yhig*8],'BackgroundColor', 'white', 'FontSize', 8);
h_ecp_info.Units = 'normalized';

%% Add Band Labels Legend

    function BandLabelsLegend(~,~)

        % Get phase selection
        selected_index = get(h_phases, 'Value');
        phases_list = get(h_phases, 'String');
        selected_phase = phases_list{selected_index};
        PhaseFilePath = fullfile(Input_Data.Phase_Folder, 'phasefiles', [selected_phase '.pha']);

        if ~exist(PhaseFilePath, 'file')
            warning('Phase file not found: %s', PhaseFilePath);
            return;
        end

        % Read and extract $PlaneReflector section
        PhaseContent = fileread(PhaseFilePath);
        token = regexp(PhaseContent, '\$PlaneReflector\s*(.*?)\$cif', 'tokens', 'once');
        if isempty(token)
            warning('No $PlaneReflector section found.');
            return;
        end
        planeReflectorSection = token{1};

        lines = regexp(planeReflectorSection, ';', 'split');
        lines = strtrim(lines);
        lines = lines(~cellfun('isempty', lines));

        % Miller indices
        millerIndices = cell(size(lines));
        lineText = cell(size(lines));
        for iis = 1:numel(lines)
            v = sscanf(lines{iis}, '%d,%d,%d');
            if numel(v) == 3
                millerIndices{iis} = v(:)';
                lineText{iis} = sprintf('{%d %d %d}', v);

            end
        end

        % Update label text & color
        maxLabels = min(numel(lineText), 10);
        for iis = 1:10
            if iis <= maxLabels
                set(h_band_labels(iis), ...
                    'String', lineText{iis}, ...
                    'ForegroundColor', Input_Data.cset(iis,:), 'BackgroundColor', '#B4B4B4');
            else
                set(h_band_labels(iis), ...
                    'String', '', ...
                    'ForegroundColor', 'w', 'BackgroundColor', 'w');  % or gray it out
            end
        end

    end


end