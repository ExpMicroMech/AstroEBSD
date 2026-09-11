%% This deck is supposed, to bundle the various methods, that can be used 
% to get from the gnomonoic projection, i.e. the screen of the EBSD
% Detector to the sphere or to the stereographic projection
% I will try, to make this more understandable and outsurce a lot of stuff
% into funcitons, such that this script will work more modularly
%       Needed Inputs:  1 or more EBSD patterns
%                       Orientation of the EBSD pattern (with frame of
%                       reference convention known)
%                       Geometry -> Projection center coordinates in Astro
%                       EBSD convention (see "Which way is up", doi ...)
%                       Phase information, to get the fundamental zone
%                       (should be provided in the AstroEBSD phase folder
%                       as a .pha file). A Dynamic Template is helpful for
%                       comparison
%
%       Outputs:        Depends on what you want:
%                       We aim to have the sphere and the stereogram, both
%                       in fourier space (spherical harmonics) and
%                       discretely
%
%
%
%           Author of this script is Lukas Berners, with some
%           functionallity used from the initial script used in (Tianbis
%           Paper), by Ben Britton and Tianbi Zhang
%
%
%
%
%





%% clear and reset
% clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;
%Load Astro and MTEX if not loaded yet % I tried to set up a different script, so that one does not have to change the file name in every file
% InputUser.MTEX_loc='C:\Users\berners\mtex-5.7.0';
% InputUser.MTEX_loc='C:\Users\berners\MTEX_dev\mtex';
% try EBSD;
% catch
%     run(fullfile(MTEX_loc,"startup.m"));
% end
% 
% try astro_loadcheck
% catch
%     run(fullfile(Astro_loc,"start_AstroEBSD.m"));
% end

%% Set Up our inputs

InputUser.Phase_Folder = Phase_Folder;
% run([InputUser.Astro_loc,'\start_AstroEBSD']);
% error('end')
InputUser.Phase = {'Aluminum'};
InputUser.PatternFile = {'stitched_flat_gaus.tif'}; % expand this to point to a folder in the future


InputUser.Pattern = double( flipud( imread( InputUser.PatternFile{1} ) ) );
InputUser.Pattern=InputUser.Pattern-mean(InputUser.Pattern(:));
InputUser.Pattern=InputUser.Pattern/std(InputUser.Pattern(:));
InputUser.eangs = [140.4739,23.8782,246.8339];
% InputUser.eangs = [254.5510 33.2576 115.8512];
% eangs_pat1=[140.4739,23.8782,246.8339];%demo
% pc_pat1=[0.5454366,0.4350097,0.7578597];
InputUser.eangs = InputUser.eangs*pi/180;
InputUser.PC_in = [0.5454366,0.4350097,0.7578597];
% InputUser.PC_in =[0.4839393 0.4376833 0.5816051];
InputUser.pattern_info = struct();
InputUser.pattern_info.size = size( InputUser.Pattern );
MicroscopeData.TotalTilt = 0;


% eangs_pat1=; pc_pat1=[0.4839393 0.4376833 0.5816051];
RTM.screensize=size(InputUser.Pattern,1);
RTM.SamplingFreq=8;
RTM.iterations=10;
RTM.LPTsize=RTM.screensize;
% error('end')
%%
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num,RTM.RTM_info] = Phase_Builder_RTM( InputUser.Phase,InputUser.Phase_Folder );
[InputUser.screen_int] = Cube_Generate(RTM.RTM_info.bin_file, RTM.RTM_info.isHex);
%%

InputUser.geometry = EBSP_Gnom(InputUser.pattern_info, InputUser.PC_in)
InputUser.cs=loadCIF(RTM.RTM_info.cif_file);
% InputUser=refineInputPattern(InputUser,RTM);
InputUser.gmatrix_refined = conv_EA_to_G(InputUser.eangs)     
ang_coverage=atan(1/double(InputUser.PC_in(3)));
grid_res=(ang_coverage/double(InputUser.pattern_info.size(1)));
%% Build an gridded Interpolant from the pattern, that can be used with the discrete method.
InputUser.geometry = EBSP_Gnom(InputUser.pattern_info, InputUser.PC_in)
InputUser.gaussianRange=5;
InputUser.gaussianRangeLarge =5 ;
InputUser.disc_sizeLarge = 0.95;
InputUser.disc_size=0.3;
%%

Settings_Cor.MeanCentre=1;
Settings_Cor.gfilt = 1;
% Settings_Cor.gfiltMasked=1;
Settings_Cor.gfilt_s=4;
Settings_Cor.hotpixel =1;
Settings_Cor.hot_thresh = 0.5;
% InputUser.PatternMask = rescale(-gauss2_disc(InputUser.geometry,InputUser.gaussianRange,InputUser.disc_size)+...
% gauss2_square(InputUser.pattern_info.size,InputUser.gaussianRangeLarge,InputUser.disc_sizeLarge));

% InputUser.Pattern = reshape(normalizeddata(:,pat_index), [254 254]) ./ normfactor;

%%(InputUser.Pattern-mean(InputUser.Pattern(:)))./std(InputUser.Pattern(:))
InputUser.pattern_info = struct();
InputUser.pattern_info.size = size( InputUser.Pattern );
InputUser.PatternMask = ones(size(pat));
InputUser.PatternMask(satdomains{pat_index})=0;
se = offsetstrel('ball',5,5);
nhood =  [0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1 ; 0 1 1 1 0]
% dialte the mask a pit (actually dilation
figure
imagesc(InputUser.PatternMask)
erodedI  =InputUser.PatternMask~=0
for i = 1:10
erodedI = imerode(erodedI, nhood);
end
InputUser.PatternMask = erodedI*1
figure
imagesc(InputUser.PatternMask)
colorbar()
%%

InputUser.Pattern(satdomains{pat_index})= mean(InputUser.Pattern(InputUser.PatternMask ~= 0));
% set black bits to mean -> they should not influence the other bits to
% much then.
InputUser.Pattern = cor_hotpix(InputUser.Pattern,0.5);
[InputUser.PatternMask,num_hot]=cor_hotpix(InputUser.PatternMask,0.5);
InputUser.PatternMask = imgaussfilt(InputUser.PatternMask,InputUser.gaussianRange);

InputUser.PatternMask = rescale(InputUser.PatternMask.*gauss2_square(InputUser.pattern_info.size,InputUser.gaussianRangeLarge,InputUser.disc_sizeLarge))%.*(InputUser.Pattern~=0))
%%
% gmatrix = InputUser.gmatrix_refined
% rotation.byMatrix(gmatrix')
[gridded_interpolant, InputUser] = PatternToInterpolant(InputUser,'ApplyMask',InputUser.PatternMask,'bgCor'); %% try gaussian, imgaussfilt, rectgauss
% InputUser.Pattern_masked = 
%
% figure
% InputUser.Pattern_masked = EBSP_BGCor(InputUser.Pattern_masked,Settings_Cor)
% imagesc(InputUser.Pattern_masked)
% InputUser.Pattern_masked(satdomains{pat_index})= mean(InputUser.Pattern_masked(InputUser.PatternMask ~= 0));
figure
imagesc(InputUser.Pattern_masked)
colorbar
%%  Build the Harmonic Represnentation
pHarm1=reprojections_harm(InputUser)
%%
InputData.V_in=30E3; %voltage in eV
InputData.lattice_param_scale=1E-10; % Angstroms
[InputData.lambda] = calc_lambda(InputData.V_in); %Wavelenght
InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale; %Wavelength in angstroms
% InputUser.Phase_Folder = Phase_Folder;
InputData.maxBNum = 6; % setup for structure factor calculation
maxHKL = 10;
InputUser.cutoff = 0.2;
InputUser.maxInd = 15;

[m_unique,f_hkl] = calcStructureFactor(InputUser.Phase  , InputData.lambda_p ,maxHKL,InputUser);
I = f_hkl.*conj(f_hkl);
I_norm = I./max(I(:));
m_Plot = m_unique(I_norm>0.2);
%%
Settings_Cor.MeanCentre=1;
Settings_Cor.gfilt=1;
Settings_Cor.gfilt_s=4;
label_bandsStructureFactor(InputUser.Pattern_masked, InputUser.geometry, m_Plot, rotation.byMatrix(InputUser.gmatrix_refined'), InputData)
%% Plot the complete upper hemisphere of the spherical harmonics representation
figure
plot(pHarm1,'complete','upper')
colormap 'gray'
%% If this is not used, with the fundamental zone, it is very likely to cause memory issues
% profile -memory on
setpref('profiler','showJitLines',1);
img_out=reprojections_discrete(InputUser,gridded_interpolant,'img_size',1000,'fundamental');
% profile viewer
% 
% profile off
imwrite(img_out,'resampled_image.tif')
%%
figure
imagesc(EBSP_BGCor(img_out,Settings_Cor))
clim([-3.5 3.5])
axis('equal','xy')
colormap 'gray'

%%
r= rotation.byEuler(0,0,0)
InputUser.cs=crystalSymmetry('m-3m')
proj_ster=eangleProjection;
bandwidth=512;

[plan.S2G] = quadratureS2Grid(2*get_option('bandwidth',128));
plan.W = plan.S2G.weights(:)
InputUser.Stereo=double(imread('resampled_image.tif'));
% do some normalisation
InputUser.Stereo=InputUser.Stereo-mean(InputUser.Stereo(:));
InputUser.Stereo=InputUser.Stereo/std(InputUser.Stereo(:));
plan.S2G_rot=inv(r)*plan.S2G;

resolution=2/(size(InputUser.Stereo,1)-1);
[gr_x,gr_y]=ndgrid(-1:resolution:1,-1:resolution:1);

stereo_pattern=griddedInterpolant(gr_x,gr_y,InputUser.Stereo,'cubic','none');

%%
[stereo_x,stereo_y]=proj_ster.project(plan.S2G_rot);
values=stereo_pattern(stereo_x,stereo_y);
% values=
%%
pHarm_fromStereo=S2FunHarmonicSym.quadrature(plan.S2G_rot,values,'weights',plan.W(:),'bandwidth',bandwidth,'quadrature',InputUser.cs,'delta',0.2);
%%
% EBSP_BGCor

figure
plot(pHarm_fromStereo,'complete','upper')
colormap 'gray'

%% Test if we can extract the bands from the stereographic projection and plot it
h=Miller({1,0,0},{1,1,0},{1,1,1},InputUser.cs)
i=1
m=h(i)
omega=[pi/2-90*degree:0.05*degree:pi/2+90*degree]
omega_size=size(omega,2);
theta=(0:0.1:360)*degree;
sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
[x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
vals=reshape(vals,[],omega_size);
figure
imagesc(omega/degree,theta/degree,vals)
axis equal;
colormap 'gray'
%% compare 100 010 and 001 to see if it causes issues
figure
hold on
for i=1:length(h)
h=Miller({1,0,0},{0,1,0},{0,0,1},InputUser.cs);
%     i=1
    disp(i)
    m=h(i)
    
    sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
    [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
    vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
    vals=reshape(vals,[],omega_size);
    sum_vals=sum(vals,1)/length(theta);
    
    hold on
    plot(omega/degree,sum_vals+i,'DisplayName',char(h(i)))
    
%     title(char(m))
end
legend()

%%
h=Miller({1,0,0},{1,1,0},{1,1,1},InputUser.cs);
for i = 1:length(h)
    [~,psi_stereo{i}] = pHarm_fromStereo.symmetrise(h(i));
    psi_stereo{i}.A(1) = 0; 
    [~,psi_harm{i}] = pHarm1.symmetrise(h(i));
    psi_harm{i}.A(1) = 0;   
end
%%
figure

omega=[pi/2-20*degree:0.05*degree:pi/2+20*degree];
omega_size=size(omega,2);
theta=(0:0.05:360)*degree;
for i = 1:length(h)
    
%   plot the harmonics
    subplot(1,3,i)
    profile_stereo=psi_stereo{i}.eval(cos(omega)) 
    profile_stereo=profile_stereo-mean(profile_stereo);
    profile_stereo=profile_stereo/std(profile_stereo);
    profile_harm=psi_harm{i}.eval(cos(omega))
    profile_harm=profile_harm-mean(profile_harm);
    profile_harm=profile_harm/std(profile_harm);
%     plot(omega/degree,profile_stereo)
    hold on
    plot(omega/degree,profile_harm,'DisplayName','original','LineWidth',2)
    plot(omega/degree,profile_stereo,'DisplayName','stereo','LineWidth',2)
% and then the discrete profiles
    m=h(i)
    sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
    [x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
    vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
    vals=reshape(vals,[],omega_size);
    sum_vals=sum(vals,1)/length(theta);
    sum_vals= sum_vals-mean(sum_vals);
     sum_vals=sum_vals/std(sum_vals);
    
    hold on
    plot(omega/degree,sum_vals,'DisplayName','Discrete','LineWidth',2)
 % psi_master_nb6co7_bruker{i}.eval(cos(omega)) psi_master_nb6co7_em{i}.eval(cos(omega))
 legend()
 title(char(h(i)))
% exportgraphics(gcf,sprintf('bandprofiles_division%s_sim.png',char(h(i))),'Resolution',300)
% close all
end
%%
for i = 1:length(h)
% m=Miller(1,0,0,InputUser.cs)
m=h(i)
omega=[pi/2-20*degree:0.005*degree:pi/2+20*degree]
omega_size=size(omega,2);
theta=(0:0.05:360)*degree
sampling=circles(vector3d(m.x,m.y,m.z),theta,omega);
[x_sample_discrete,y_sample_discrete]=proj_ster.project(sampling);
vals=stereo_pattern(x_sample_discrete,y_sample_discrete);
vals=reshape(vals,[],omega_size);
% figure
% imagesc(omega/degree,theta/degree,vals)
% axis equal;
% colormap 'gray'

figure
plot(omega/degree,sum(vals,1))
title(char(m))
end
% fv=QuadCubeMesh; % also try this with `QuadRhombDodecMesh.m`
% %%
% % Subdivide base mesh and visualize the results%%
% pattern_info = struct();
% pattern_info.size = [200 200];
% pc_in=[0.5 0.5 0.5]
% [EBSD_geometry ] = EBSP_Gnom(pattern_info,pc_in);
% z_comp=ones(size(EBSD_geometry.xpts_screen(:)));
% %%
% faces=cat(3,[EBSD_geometry.xpts_screen(:),EBSD_geometry.ypts_screen(:),z_comp],...
%         [EBSD_geometry.xpts_screen(:),EBSD_geometry.ypts_screen(:),-z_comp],...
%         [EBSD_geometry.xpts_screen(:),z_comp,EBSD_geometry.ypts_screen(:)],...
%         [EBSD_geometry.xpts_screen(:),-z_comp,EBSD_geometry.ypts_screen(:)],...
%         [z_comp,EBSD_geometry.xpts_screen(:),EBSD_geometry.ypts_screen(:)],...
%         [-z_comp,EBSD_geometry.xpts_screen(:),EBSD_geometry.ypts_screen(:)]);
%     %%
% 
% Cube_Sample()
% for i=1:6
%     scatter3(faces(
function c= circles(normals,theta,omega)
    for i = 1:length(normals)

      n = normals.subSet(i);

      % generate circles
      c = axis2quat(n,theta) * axis2quat(orth(n),omega)*n;

      % plot circles
%       h = [h,line(c,varargin{:},'hold','doNotDraw')]; %#ok<AGROW>

    end
end
%%
function G=gauss2_disc(EBSDgeometry,sigma,r_thresh)
    % r_thresh in gnomonic coordinates
      % = EBSDgeometry.size;
     G = zeros(size(EBSDgeometry.xpts_screen));
     r_grid=sqrt(double(EBSDgeometry.xpts_screen).^2+double(EBSDgeometry.ypts_screen).^2);
     G(r_grid<=r_thresh) = 1;
     G=imgaussfilt(G,sigma);
     %%
end

function G=gauss2_square(size2,sigma,disc_size)
     G = zeros(size2(1),size2(2));
     G( uint16(size2(1)*(1-disc_size)/2) : uint16(size2(1)*(1-(1-disc_size)/2)) ,uint16( size2(2)*(1-disc_size)/2) : uint16( size2(2)*(1-(1-disc_size)/2)) ) = 1;
     
%%
     G=imgaussfilt(G,sigma);
%      Y(size(2)*(disc_size/2):size(2)*(1-disc_size/2)=1
end

function [EBSP2,num_hot]=cor_hotpix(EBSP,hthresh)

%Image correction - modified from JH code
%median filter

EBSP_med = medfilt2(EBSP);

%hot pixel correct
h_pix=find(abs(EBSP-EBSP_med)>hthresh); %1000 is chosen as an arbitary number
EBSP2=EBSP;
EBSP2(h_pix)=EBSP_med(h_pix);
num_hot=size(h_pix);

end