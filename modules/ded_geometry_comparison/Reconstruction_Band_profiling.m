%% Reconstruction_Band_Profiling.m
% Created by Lukas Berners, 2024
% Edited and annotated by Lukas Berners and Tianbi Zhang

% This script takes a Kikuchi pattern, index it, reproject it onto the
% diffraction sphere and create a stereogram. This stereogram is
% approximated by a spherical harmonic series for band profiling.

% This script accompanies the following article(s):
% "Comparison of Kikuchi Diffraction Geometries in Scanning Electron Microscope"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% Requirements: 
% (1) MATLAB toolboxes: image processsing, statistics and machine learning,
% parallel computing
% (2) AstroEBSD package - this script is available as a part of the latest
% AstroEBSD distribution, specifically the TKD-DED branch.
% (3) MTEX toolbox (https://mtex-toolbox.github.io/)
% (4) The custom class detector (written by Lukas Berners) included in this
% release.

%%
home;
clear; 
close all;

%% 
InputUser.MTEX_loc='C:\Users\billy\Documents\MATLAB\mtex-5.11.1\';

try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end

InputUser.Phase_Folder="D:\AstroEBSD\phases\";
InputUser.Astro_loc='D:\AstroEBSD\';

% start astro if needed
try astro_loadcheck
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end

%% load input pats
% Euler angles (Bunge) and PC - they must be very accurate for the
% reprojection.
% change accordingly
ebsd_pat=double(imread('C:\Users\billy\OneDrive - UBC\PhD\TKD\Kikuchi_Geometry_Comparison\Al_offaxis_tkd\s15_100f_002s_bgcor.tif'));
eangs_pat1=[139.7766 53.1139 227.3212] * degree; %demo
pc_pat1=[0.45425 -0.15789 0.66814];

%
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% setup phase and crystal structure
InputUser.Phase_Input={'Aluminium'};
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
cs=loadCIF(RTM_info.cif_file); 
% [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM({InputUser.Phase_Input{1}},InputUser.Phase_Folder);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Bands do you want to analyse
h=Miller({0,0,2},{2,2,0},{1,1,1},cs,'hkl');
cset=[1,0,0;
    0,1,0;
    0,0,1]; % colors for band overlay

gmatrix=RTM.Rz(eangs_pat1(3))*RTM.Rx(eangs_pat1(2))*RTM.Rz(eangs_pat1(1));
pc_patin=pc_pat1;
PatternIn=ebsd_pat;
ori=orientation.byMatrix(gmatrix,cs);
pattern_info=struct;
pattern_info.size=size(PatternIn);
pattern_info.size=uint16(pattern_info.size);

bandwidth = 384;
%% corrections
Settings_Cor.resize=0; %resize correction (1)
Settings_Cor.size=[pattern_info.size(1) pattern_info.size(2)];
Settings_Cor.SquareCrop=0;
Settings_Cor.Square=0;

% Normalise intensities
[PatternIn,Settings_Cor ] = EBSP_BGCor( PatternIn,Settings_Cor);

[EBSD_geometry ] = EBSP_Gnom( pattern_info,pc_patin);

%% get the detector -> this is a class
det = detector(pattern_info.size(2),pattern_info.size(1),pc_patin(3),[pc_patin(1),pc_patin(2)]);

%% loop over all symmetric equivalent patterns and get the grids and intensities
all_harmonics=struct();
count=1;
% tic

for i=1:size(cs.rot(:),1)
    [plan,v_det]= det.pattern2Fun_multi_pat(flipud(PatternIn)','bandwidth',bandwidth,'quadrature','delta',0.2,'ori',cs.rot(i)*ori);
    all_harmonics(count).plan=plan;
    all_harmonics(count).v_det=v_det;
    count=count+1;
end

%% Sum the intensities
multiples=zeros(size(v_det,1),size(all_harmonics,2));
intensities=zeros(size(multiples));

for i=1:size(all_harmonics,2)
    multiples(:,i)= all_harmonics(i).plan.mult;
    intensities(:,i)=all_harmonics(i).v_det;
end

multiples=sum(multiples,2);
zero_mults=multiples==0;
multiples=1./multiples;
multiples(zero_mults)=0;

intensities=sum(intensities,2);
inp=multiples.*intensities;

% save(inp, "thisstereogram.mat");
% If the pattern doesn't cover the entire fundamental zone and you wish to
% do some averaging first, stop here and export "inp", do this for multiple 
% patterns, perform the corresponding normalization steps as needed, 
% then come back.

% std(inp)
% inp=inp-mean(inp(:));
% inp=inp/std(inp(:));
% input=inp;

%% Calculat the harmonic pattern

pHarm = S2FunHarmonicSym.quadrature(plan.S2G,inp,'weights',plan.W,'bandwidth',bandwidth ,'quadrature',cs,'delta',0.2,'complete');
figure
plot(pHarm,'complete')
hold on
plot(Miller(1,1,1,cs).symmetrise(),'plane')
hold on
plot(Miller(1,1,0,cs).symmetrise(),'plane')
hold on
plot(Miller(1,0,0,cs).symmetrise(),'plane')
colormap gray;

%% build the base pattern based on simulation
isHex=screen_int.isHex;
InputUser.Phase_Input={'Aluminium'};

[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
cs=loadCIF(RTM_info.cif_file); 
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
template_al = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));
templateHarm_al = S2FunHarmonicSym.quadrature(template_al,cs,'bandwidth', bandwidth);
template_direct_al=Cube_Sample(plan.S2G.x(:),plan.S2G.y(:),plan.S2G.z(:),screen_int,isHex);

%% Plot the spheres
figure
f = subplot(1,1,1);
plot(templateHarm_al,'complete','upper','colorrange',[-3 3])
% f.clim(-3,3)
title('pattern from simulation')
colormap gray;
f=mtexColorbar();
% subplot(1,2,2)
figure;
subplot(1,1,1);
plot(pHarm,'complete','upper','colorrange',[-3 3])
title('exp pat')
colormap gray;
% f.clim(-3,3)
mtexColorbar();

%% Get band profiles
for i = 1:length(h)
    [~,profile_simulation{i}] = templateHarm_al.symmetrise(h(i));
    profile_simulation{i}.A(1) = 0;
    [~,profile_pattern{i}] = pHarm.symmetrise(h(i));
    profile_pattern{i}.A(1) = 0;
end

%% plot bands on pattern

% As you can see, the approximated pattern is blurry, and this figure is
% just for reference. For overlaying on a high resolution stereogram,
% consider using the plotting block in "Reconstruction_from_pattern".

figure; 
  plot(pHarm,'DisplayName','simulation','linewidth',2,'upper','complete')
  colormap gray
for i = 1:length(h)
    hold on
    mils=h(i).symmetrise;
    for j=1:size(mils)
        if j==1
            plot(mils(j),'DisplayName',char(h(i)),'linewidth',2,'Plane','color',cset(i,1:3))
        
        else
            plot(mils(j),'HandleVisibility','off','linewidth',2,'Plane','color',cset(i,1:3))
        end
    end
end
legend();

%% plot band profiles 
% Note that this is for a single pattern
% For plotting Figure 5, see Plot_Band_Profiles.m in this folder.
figure
for i = 1:length(h)
  nexttile;
  
  plot(profile_simulation{i},'DisplayName','simulation','linewidth',3)
  hold on
  plot(profile_pattern{i},'DisplayName','experiment','linewidth',3)
  title(char(h(i)))
  xlim([80,100])
end
legend();






