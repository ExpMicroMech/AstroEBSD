%% Plot_band_profiles
% Created by Tianbi Zhang in November 2024

% This script plots band profiles obtained from multiple patterns,
% corresponding to different diffraction geometries. These profiles are in
% the form of spherical harmonics.

% This script only plots the profile, it does not calculate the profiles.
% For that, use "Reconstruction_Band_profiling.m".

% This script accompanies the following article(s):
% "Comparison of Kikuchi Diffraction Geometries in Scanning Electron Microscope"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% Requirements: 
% (1) AstroEBSD package - this script is available as a part of the latest
% AstroEBSSD distribution.
% (2) MTEX toolbox (https://mtex-toolbox.github.io/)
% N.B. You may still need other dependencies to run the otehr scripts in
% AstroEBSD.

%% Load AstroEBSD and MTEX - edit the directories as necessary
InputUser.Astro_loc='C:\Users\billy\Documents\GitHub\AstroEBSD_v2_standalone\';
InputUser.mtex_location='C:\Users\billy\Documents\MATLAB\mtex-5.11.1';

%start mtex if needed
try EBSD;
catch
    run(fullfile(InputUser.mtex_location,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end

% which bands to plot?
cs=loadCIF('Al-Aluminum'); 
h=Miller({0,0,2},{0,2,2},{1,1,1},cs,'hkl');

%% Bragg Angles


% Band profiling by the spherical harmonics approach doesn't give the band
% opening directly, so we need to centre the plot at 90 degree

centre_angle = 90;

% We calculated the angles elsewhere and just put the value down here.
braggs = [1.2147, 2.4351, 3.6497; %002, 004, 006
    1.7189, 3.4435, 5.1681; % 022, 044, 066
    1.2147, 2.4351, 3.6497; %111, 222, 333
    ];

% Color set for the shading
cset =  parula(256);

% Plot limits
ylims = [-2 3;-1.5 1.5; -2 4];

%% Plot
fig = figure;
for i = 1:length(h)
    
  ax1 = subplot(3,1,i, 'Parent',fig);
  
  hold on
  plot(ebsd_prof{i},'DisplayName','EBSD','linewidth',3)
  plot(tkdon_prof{i},'DisplayName','RKD','linewidth',3)
  plot(tkdon_prof{i},'DisplayName','On-axis TKD','linewidth',3)
  plot(tkdoff_prof{i},'DisplayName','Off-axis TKD','linewidth',3)
  plot(sim_prof{i},'DisplayName','Simulation','linewidth',2, 'color','black')
  title(char(h(i)))
  
  xlim([80,100]);
  set(ax1,'XTickLabel', {}, 'YTickLabel', {});

    % axis off;
    
    % uncomment to get the legend - but it will be messed up due to the shades
    % legend();
    
    % plot the shades
  for j=1:3
      
      x_fill = [centre_angle-braggs(i,j) centre_angle-braggs(i,j) centre_angle+braggs(i,j) centre_angle+braggs(i,j)];
      y_fill = [-5 5 5 -5];
      fill(x_fill, y_fill, cset(j,1:3), 'FaceAlpha',0.08, 'EdgeColor','none');

  end
  ylim([ylims(i,1),ylims(i,2)]);
  
end
