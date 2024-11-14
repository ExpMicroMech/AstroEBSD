%% Reconstruction_from_pattern
% Created by T. Ben Britton in June 2023
% Edited and annotated by Tianbi Zhang in November 2024

% This script takes a Kikuchi pattern, index it, reproject it onto the
% diffraction sphere and create a stereogram.

% This script accompanies the following article(s):
% "Multi-exposure diffraction pattern fusion applied to enable wider-angle 
% transmission Kikuchi diffraction with direct electron detectors"
% Tianbi Zhang, T. Ben Britton 
% "xxx"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% Requirements: 
% (1) MATLAB toolboxes: image processsing, statistics and machine learning,
% parallel computing
% (2) AstroEBSD package - this script is available as a part of the latest
% AstroEBSD distribution.
% (3) MTEX toolbox (https://mtex-toolbox.github.io/)

% If you wish to sum up and average multiple reconstructions, simply run on
% different patterns and save "pat_maps_norm", then do post-processing.
% For reconstruction then pband profiling, see Reconstruction_Band_profiling.m.

%% Start of the script
close all;
clear;
home;

%% start astro & mtex - please edit the paths
location_astro='D:\AstroEBSD\';
location_mtex='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';

%start mtex if needed
try EBSD;
catch
    run(fullfile(location_mtex,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(location_astro,"start_AstroEBSD.m"));
end

%% load the pattern 
% change accordingly

path= 'C:\Users\billy\OneDrive - UBC\PhD\TKD\Kikuchi_Geometry_Comparison\Al_offaxis_tkd\';
patternName = 's15_100f_002s_bgcor.tif';
expPatPath = fullfile(path, patternName);
expPat=double(flipud(imread(expPatPath))); % this is the demo TKP

%% Enter pattern information and set up indexing
% Euler angles (Bunge) and PC - they must be very accurate for the
% reprojection.
% change accordingly based on the information provided
eangs_pat1=[ 139.7766 53.1139 227.3212 ]; %demo, degs
pc_pat1=[0.45425 -0.15789 0.66814]; %demo

MicroscopeData.TotalTilt=0;

% Normalise intensities
[EBSP_One.PatternIn,Settings_Cor ] = EBSP_BGCor( expPat,[]);
InputUser.Phase_Folder = fullfile(location_astro,'phases');
InputUser.Phase_Input  = {'Al_bwkd'};
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
cs_al=loadCIF(RTM_info.cif_file); 

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% Figure 1: plot the original pattern
figure;
imagesc(expPat); axis image; axis xy; colormap('gray');
[patsize,~] = size(expPat);

%% Simulate the pattern

eangs=eangs_pat1*pi/180;
PC_simulation=pc_pat1;
pattern_info=struct;
pattern_info.size=size(expPat);
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

% % InputUser.Phase_Folder
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},[location_astro 'phases']);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
% 
[EBSD_simulation ] = EBSP_Gnom( pattern_info,PC_simulation); %you can change PC_in if you want

gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));

% If you want to adjust where the reprojection is to match the location of
% the fundamental zone, use the following (change n to see if it works for you):
% gmatrix = cs_al.rot(n).matrix*RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));

%% Figure 2: plot alongside a simulated pattern (refinement skipped)

rotmat=gmatrix*Detector_tilt;
[Pat_sim_ref,r_exp]= EBSP_gen( EBSD_simulation,rotmat,screen_int); %generate the EBSP for this iteration

figure;
subplot(1,2,1);
pPattern(Pat_sim_ref,EBSD_simulation);
subplot(1,2,2);
pPattern(expPat,EBSD_simulation);

%% Now try to plot a sterogram of the dynamically simulated pattern + 
% Figure 3: the full simulated stereogram
num_stereo=2001;
stereo_range=1;
stereo_xl=linspace(-stereo_range,stereo_range,num_stereo);
stereo_yl=linspace(-stereo_range,stereo_range,num_stereo);

[stereo_x,stereo_y]=meshgrid(stereo_xl,stereo_yl);

stereo_r=stereo_x.^2+stereo_y.^2;

stereo_nu=2./(stereo_r+1);

stereo_vx=stereo_nu.*stereo_x;
stereo_vy=stereo_nu.*stereo_y;
stereo_vz=-1+stereo_nu;

stereo_vl=sqrt(stereo_vx.^2+stereo_vy.^2+stereo_vz.^2);
stereo_vx=stereo_vx./stereo_vl;
stereo_vy=stereo_vy./stereo_vl;
stereo_vz=stereo_vz./stereo_vl;


%sample from the sphere
[stereo_pat] = Cube_Sample(stereo_vx(:),stereo_vy(:),stereo_vz(:),screen_int,0);
figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat,num_stereo,num_stereo));
axis image; axis xy; colormap('gray');

%% Now take the experimental pattern and generate a look up, so we can do the spherical sampling...
% Figure 5: experimental pattern overlaid on the stereogram
% rotmat=EBSP_Fe_One.rotdata{1}.detector;

r_lambda=1./(r_exp(:,3)+1);

%convert from r_exp to the stereogram
P_rx=r_exp(:,1).*r_lambda;
P_ry=r_exp(:,2).*r_lambda;

% switch to simulation
pat1_plot=Pat_sim_ref;
pat1_plot=pat1_plot/std(expPat(:));
pat1_plot=pat1_plot-mean(pat1_plot(:));

% create a intensity normalized stereogram
stereo_pat_plot=stereo_pat;
stereo_pat_plot=stereo_pat/std(stereo_pat_plot(:));
stereo_pat_plot=stereo_pat_plot-mean(stereo_pat_plot(:));

% use the interpolant to generate the pattern
% create the interpolant
stero_interp_exp=scatteredInterpolant(P_rx(:),P_ry(:),pat1_plot(:),'natural','none');

%read the intepolant
stereo_remap=stero_interp_exp(stereo_x(:),stereo_y(:));

stereo_remap_2d=reshape(stereo_remap,num_stereo,num_stereo);
stereo_mask_2d=isnan(stereo_remap_2d);

% plot the outline of the detector (experimental pattern)
line_top=[EBSD_simulation.xpts_screen(:,1),EBSD_simulation.ypts_screen(:,1)];
line_bottom=flipud([EBSD_simulation.xpts_screen(:,end),EBSD_simulation.ypts_screen(:,end)]);
line_right=[EBSD_simulation.xpts_screen(end,:);EBSD_simulation.ypts_screen(end,:)]';
line_left=flipud([EBSD_simulation.xpts_screen(1,:);EBSD_simulation.ypts_screen(1,:)]');

line_allx=[line_top(:,1);line_right(:,1);line_bottom(:,1);line_left(:,1)];
line_ally=[line_top(:,2);line_right(:,2);line_bottom(:,2);line_left(:,2)];
line_allz=line_ally*0+1;

%rotate these points around the sphere
line_r = [line_allx(:), line_ally(:), line_allz(:)].*1./sqrt((line_allx(:).^2+line_ally(:).^2+1));
line_r2 = line_r*rotmat';

%convert these points to the stereogram
line_r_lambda=1./(line_r2(:,3)+1);

%convert from r_exp to the stereogram
line_P_rx=line_r2(:,1).*line_r_lambda;
line_P_ry=line_r2(:,2).*line_r_lambda;

figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray'); axis off;
hold on;

%plot the blue box
plot(line_P_rx,line_P_ry,'b','LineWidth',2);

%
%do the fundamental zone
%write three axes, and construct the lines that go between one and the next
ax1=[0 0 1];
ax2=[1 1 1];
ax3=[1 0 1];

%normalize length
ax1=ax1./norm(ax1);
ax2=ax2./norm(ax2);
ax3=ax3./norm(ax3);
nl=100;

%inscribe along the great circles and equispace
ang_12=acos(dot(ax1,ax2));
ax_12=cross(ax1,ax2);
ax_12=ax_12/norm(ax_12);
ax_2p=cross(ax_12,ax1);
ax_2p=ax_2p/norm(ax_2p);
sample_spacing12=linspace(0,ang_12,nl);
line_12=ax1'*cos(sample_spacing12)+ax_2p'*sin(sample_spacing12);

ang_23=acos(dot(ax2,ax3));
ax_23=cross(ax2,ax3);
ax_23=ax_23/norm(ax_23);
ax_3p=cross(ax_23,ax2);
ax_3p=ax_3p/norm(ax_3p);
sample_spacing23=linspace(0,ang_23,nl);
line_23=ax2'*cos(sample_spacing23)+ax_3p'*sin(sample_spacing23);

ang_31=acos(dot(ax3,ax1));
ax_31=cross(ax3,ax1);
ax_31=ax_31/norm(ax_31);
ax_1p=cross(ax_31,ax3);
ax_1p=ax_1p/norm(ax_1p);
sample_spacing31=linspace(0,ang_31,nl);
line_31=ax3'*cos(sample_spacing31)+ax_1p'*sin(sample_spacing31);

%projected on stereogram
line_set=[line_12';line_23';line_31'];
line_set_lambda=1./(line_set(:,3)+1);
line_set_rx=line_set(:,1).*line_set_lambda;
line_set_ry=line_set(:,2).*line_set_lambda;

% uncomment the following to plot the fundamental zone
% you may need to edit the fundamental zone for different crystal
% structures (default is cubic)
plot(line_set_rx(1:end),line_set_ry(1:end),'y', 'LineWidth',2);

figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray'); axis off;
hold on;
i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2d);
axis image; axis xy; colormap('gray'); axis off;
i1.AlphaData=1-stereo_mask_2d;
%plot the blue box
plot(line_P_rx,line_P_ry,'b','LineWidth',2);
plot(line_set_rx(1:end),line_set_ry(1:end),'y', 'LineWidth',2);

%% lets inflate this up to the full pattern, if we can remember how

% this is actually redundant now that astro also reads directly from cif, can consolidate...

EBSP_av=EBSD_simulation;
[Pat_sim_astroind,r_exp]= EBSP_gen( EBSD_simulation,rotmat,screen_int); %generate the EBSP for this iteration

r_n=zeros(size(r_exp,1),size(r_exp,2),24);

[n_sym,~] =size(cs_al.rot);
for n=1:n_sym
    sym_matrix=cs_al.rot(n).matrix;
    r_n(:,:,n) = r_exp*sym_matrix;
end

nline=20;
line_P_rxn=zeros(4*nline,n_sym);
line_P_ryn=zeros(4*nline,n_sym);
stereo_remap_maskp=zeros(num_stereo,num_stereo,n_sym);
stereo_remap_2dnp=stereo_remap_maskp;

stereo_x_ok=stereo_x(:);
stereo_y_ok=stereo_y(:);

% calculate and plot each symmetrically equivalent pattern in the
% stereogram
parfor p=1:n_sym

    x_lookup=EBSP_av.xpts_screen;
    y_lookup=EBSP_av.ypts_screen;

    x_lookup=x_lookup(:,2:end-1);
    y_lookup=y_lookup(:,2:end-1);
    pat1_plot_n=pat1_plot(:,2:end-1);

    x_lookup=x_lookup(2:end-1,1);
    y_lookup=y_lookup(2:end-1,1);
    pat1_plot_n=pat1_plot_n(2:end-1,:);

    r = [x_lookup(:), y_lookup(:), y_lookup(:)*0+1].*1./sqrt((x_lookup(:).^2+y_lookup(:).^2+1));
    r_exp2 = r*rotmat';

    r_3=r_n(:,3,p);
    r_3(r_3<0)=-r_3(r_3<0);
    r_lambda=1./(r_3+1);

    %convert from r_exp to the stereogram
    P_rxn=r_n(:,1,p).*r_lambda;
    P_ryn=r_n(:,2,p).*r_lambda;

    % interpolate
    stero_interp_exp_n=scatteredInterpolant(P_rxn(:),P_ryn(:),pat1_plot(:),'natural','none');

    %plot the outline of the detector
    line_top=[EBSD_simulation.xpts_screen(:,1),EBSD_simulation.ypts_screen(:,1)];
    line_bottom=flipud([EBSD_simulation.xpts_screen(:,end),EBSD_simulation.ypts_screen(:,end)]);
    line_right=[EBSD_simulation.xpts_screen(end,:);EBSD_simulation.ypts_screen(end,:)]';
    line_left=flipud([EBSD_simulation.xpts_screen(1,:);EBSD_simulation.ypts_screen(1,:)]');
    
    s2=linspace(1,patsize,nline);
    s2(2:end-1)=round(s2(2:end-1));

    line_allx=[line_top(s2,1);line_right(s2,1);line_bottom(s2,1);line_left(s2,1)];
    line_ally=[line_top(s2,2);line_right(s2,2);line_bottom(s2,2);line_left(s2,2)];
    line_allz=line_ally*0+1;

    %rotate these points around the sphere
    line_r = [line_allx(:), line_ally(:), line_allz(:)].*1./sqrt((line_allx(:).^2+line_ally(:).^2+1));
    line_r2 = line_r*rotmat'*cs_al.rot(p).matrix;
    line_r2_3=line_r2(:,3);
    line_r2_3(line_r2_3<0)=-line_r2_3(line_r2_3<0);
    %convert these points to the stereogram
    line_r_lambda=1./(line_r2_3+1);

    %convert from r_exp to the stereogram
    line_P_rx=line_r2(:,1).*line_r_lambda;
    line_P_ry=line_r2(:,2).*line_r_lambda;
    
    %create the blank frame
    stereo_remap_2dn=zeros(num_stereo,num_stereo);
    stereo_remap_mask=ones(num_stereo,num_stereo);
    num_el=numel(stereo_remap_mask);
    % %find the points in the detector frame
    % [in,on] = inpolygon(stereo_x_ok,stereo_y_ok,line_P_rx,line_P_ry);
    in=true(num_el,1);
    stereo_remap_mask_0=stereo_remap_mask*0;
    %read the intepolant
    stereo_remap_n=stero_interp_exp_n(stereo_x_ok(in),stereo_y_ok(in));
    stereo_remap_2dn(in)=stereo_remap_n;
    stereo_remap_mask(in)=stereo_remap_mask_0(in);

    % stereo_remap_2dn=reshape(stereo_remap_n,num_stereo,num_stereo);
    % stereo_mask_2dn=isnan(stereo_remap_2dn);
    
    % Figure 6,7
    %plot the base dynamical pattern
    % Comment out to recude the number of windows 
    % figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
    % axis image; axis xy; colormap('gray');
    % hold on;
    % 
    % %plot the remapped experiment
    % i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2dn);
    % axis image; axis xy; colormap('gray');
    % i1.AlphaData=1-stereo_remap_mask;
    % 
    % %plot the blue box
    % plot(line_P_rx,line_P_ry,'b');
    % xlim([-1 1]);
    % ylim([-1 1]);

    title(int2str(p));

    line_P_rxn(:,p)=line_P_rx;
    line_P_ryn(:,p)=line_P_ry;
    stereo_remap_maskp(:,:,p)=stereo_remap_mask;
    stereo_remap_2dnp(:,:,p)=stereo_remap_2dn;
end

%% Plot the sterogram and try to work out the points for this projection

%read the experimental pattern
pat1_plot2=pat1_plot;
pat_c=pat1_plot2(:);

pat_maps=zeros(2001,2001,n_sym);
pat_mapn=pat_maps;

%pick the symemtry
parfor p=1:n_sym
    %extract the points of this pattern
    %rotate these points of the detector around the sphere
    box_x=EBSD_simulation.xpts_screen(:);
    box_y=EBSD_simulation.ypts_screen(:);
    box_z=box_y*0+1;

    % box_r=[box_x,box_y,box_z];
    box_r = [box_x(:), box_y(:), box_z(:)].*1./sqrt((box_x(:).^2+box_y(:).^2+1));

    box_r2=box_r*rotmat'*cs_al.rot(p).matrix;
    box_r2x=box_r2(:,1);
    box_r2y=box_r2(:,2);
    box_r2z=box_r2(:,3);

    %now find those which are in the positive and negative hemisphere ('p'
    %and 'n' points) based upon their z coords
    box_r2z_pz=find(box_r2z>0);
    box_r2z_nz=find(box_r2z<0);

    %subselect the XYZ positives
    box_r2xp=box_r2x(box_r2z_pz);
    box_r2yp=box_r2y(box_r2z_pz);
    box_r2zp=box_r2z(box_r2z_pz);

    %subselect the XYZ negatives
    box_r2xn=-box_r2x(box_r2z_nz);
    box_r2yn=-box_r2y(box_r2z_nz);
    box_r2zn=-box_r2z(box_r2z_nz);

    %extract the pattern data (we are going to use scatter for a lazy plot)
    pat_cp=pat_c(box_r2z_pz);
    pat_cn=pat_c(box_r2z_nz);

    %sterographic projection
    box_r2_lambdap=1./(box_r2zp+1);
    box_r2_xp=box_r2xp.*box_r2_lambdap;
    box_r2_yp=box_r2yp.*box_r2_lambdap;

    %find the boundary of the point cloud
    boundary_r2_p=boundary(box_r2_xp,box_r2_yp);
    boundary_r2_px=box_r2_xp(boundary_r2_p);
    boundary_r2_py=box_r2_yp(boundary_r2_p);

    %sterographic projection
    box_r2_lambdan=1./(box_r2zn+1);
    box_r2_xn=box_r2xn.*box_r2_lambdan;
    box_r2_yn=box_r2yn.*box_r2_lambdan;

    %find the boundary of the point cloud
    boundary_r2_n=boundary(box_r2_xn,box_r2_yn);
    boundary_r2_nx=box_r2_xn(boundary_r2_n);
    boundary_r2_ny=box_r2_yn(boundary_r2_n);
  
    %find the pounts within each boundary
    sx=stereo_x;
    sy=stereo_y;

    [in_p]=inpolygon(sx,sy,boundary_r2_px(:),boundary_r2_py(:));
    [in_n]=inpolygon(sx,sy,boundary_r2_nx(:),boundary_r2_ny(:));

    x_lookup=EBSD_simulation.xpts_screen;
    y_lookup=EBSD_simulation.ypts_screen;

    r = [x_lookup(:), y_lookup(:), y_lookup(:)*0+1].*1./sqrt((x_lookup(:).^2+y_lookup(:).^2+1));
    r_exp2 = r*rotmat'*cs_al.rot(p).matrix;
    %create the interpolant
    r_3=r_exp2(:,3);
    r_lambda=1./(r_3+1);

    %convert from r_exp to the stereogram
    P_rxn=r_exp2(:,1).*r_lambda;
    P_ryn=r_exp2(:,2).*r_lambda;

    % generate interpolate
    stero_interp_exp_s=scatteredInterpolant(P_rxn(:),P_ryn(:),pat1_plot(:),'natural','none');

    stero_interp_exp_sp=scatteredInterpolant(box_r2_xp(:),box_r2_yp(:),pat_cp(:),'natural','none');
    stero_interp_exp_sn=scatteredInterpolant(box_r2_xn(:),box_r2_yn(:),pat_cn(:),'natural','none');

    P_stero_xp=sx(in_p);
    P_stero_yp=sy(in_p);
    Pat_interp_p=stero_interp_exp_sp(P_stero_xp,P_stero_yp);

    P_stero_xn=sx(in_n);
    P_stero_yn=sy(in_n);
    Pat_interp_n=stero_interp_exp_sn(P_stero_xn,P_stero_yn);

    pat_map=0*sx;
    pat_num=0*sx;
    pat_num(in_p)=pat_num(in_p)+1;
    pat_num(in_n)=pat_num(in_n)+1;

    pat_map(in_p)=Pat_interp_p;
    pat_map(in_n)=Pat_interp_n;

    pat_maps(:,:,p)=pat_map;
    pat_mapn(:,:,p)=pat_num;
end

%% Normalize the overlaid stereogram and plot
pat_maps_num=sum(pat_mapn,3);
pat_maps_sum=sum(pat_maps,3);
pat_maps_norm=pat_maps_sum./pat_maps_num;

%% final stereogram
figure; 
subplot(1,3,1);
imagesc(stereo_xl,stereo_yl,pat_maps_sum); % total stereogram
axis image; axis xy; colormap('gray'); axis off;
subplot(1,3,2);
imagesc(stereo_xl,stereo_yl,pat_maps_norm); % stereogram normalized by weight
axis image; axis xy; colormap('gray'); axis off;
% sum/weight hemisphere
subplot(1,3,3); imagesc(stereo_xl,stereo_yl,pat_maps_num); % overlapped patterns
axis image; axis xy; colormap('gray'); axis off; 

% Difference between simulated and reprojected patterns
% normalize the reprojected dynamics pattern
pat_maps_dyn=reshape(stereo_pat_plot,num_stereo,num_stereo);
pat_maps_dynn=(pat_maps_dyn-mean(pat_maps_dyn(:)))./std(pat_maps_dyn(:));
pat_maps_normn=(pat_maps_norm-mean(pat_maps_norm(:)))./std(pat_maps_norm(:));

figure; imagesc(stereo_xl,stereo_yl,pat_maps_norm-pat_maps_dyn);
axis image; axis xy; colormap('gray');
hold on;
clim([-6 6]);
