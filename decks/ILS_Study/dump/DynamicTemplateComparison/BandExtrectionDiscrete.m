% Bandprofile extraction
% based on a gridded interpolant
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
InputUser.Phase_Folder=Phase_Folder
InputUser.Phase_Input={'Nb6Ni7_em'}
root_path = 'E:\sciebo\IMM\BandProfiles'
%% set Colormaps
InputData.cset=[1,0,0; % Colors of the Bands in rgb (1 =255)
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1];
InputData.lineset =["-",":","-.","--" ]
% InputData.cset=[InputData.cset;InputData.cset*0.5];
InputData.cset = tab10;
%% 
1
%% get Wavelength
InputData.V_in=20E3; %voltage in eV
InputData.lattice_param_scale=1E-10; %lattice paramters are provided in this unit, often Angstrom (1E-10);
i=1
[InputData.lambda] = calc_lambda(InputData.V_in);

%scale the wavelength to the same units as the dhkl
InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale;

%%

%% Build Crystal structure from Astro Phase file
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( {InputUser.Phase_Input{1}},InputUser.Phase_Folder );
Family_List=Crystal_Family{1};
[cs,~]=loadCIF(RTM_info.cif_file); 
cs = crystalAlignEdax(cs);
InputUser.cs = cs
%%
InputUser.maxInd = 30;
% f_j * exp( - 2*pis*i * (h*x(j) + k*y(j)+ l*z(j) ))
% InputUser
InputUser.cutoff = 0.05
InputUser.minD = 0.2
InputUser.AlignEdax =1 
[m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{1}}, InputData.lambda_p ,InputUser);

%%
cutoff = 0.2
[m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,cutoff)
%% Build the interpolant, to use for pattern projection
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
export_path=fullfile(root_path,InputUser.Phase_Input{i})
mkdir(export_path)
phaseName = InputUser.Phase_Input{i};
phaseNameClean = cleanNameLatex(phaseName);
%% Build the dynamical template interpolant
screen_int=Cube_Generate(RTM_info.bin_file,0); 
%% Setup some Projection center parameters (ecp will likely be pcx=0.5 pcy=0.5 pcz=3.5-5)
PatternInfo.screensize=1024;
InputUser.pattern_info.size = [PatternInfo.screensize, PatternInfo.screensize] %% size of plotted pattern
PC= [0.5 0.5 0.6]
InputUser.PC_in = PC%% pcx pcy pcz (or pcz also known as Detector Distance)
EBSD_geometry=EBSP_Gnom(PatternInfo,PC); %Build the geometry
%% Plot the Spherical pattern (on a spher, will be able to rotate)   
[x,y,z]=sphere(1000);
i_data=Cube_Sample(x,y,z,screen_int,screen_int.isHex);
  % Plot the sphere for rotation
figure
surf(x,y,z,i_data,'EdgeColor','none')
colormap gray;
axis equal;
title(phaseNameClean)
exportgraphics(gcf,fullfile(export_path,sprintf('%s_sphere.png',InputUser.Phase_Input{i}))) 
% end
% close all
%%
%% Specify zone axes that should be in the centre
h=Miller({0,0,0,1},cs)%,{0,0,1},{1,1,1},cs)
h = Miller({-2,0,2,1,'UVTW'},cs)%% h is for the zone axis along which we want to observe the ecp
zoneAxis = Miller({0,0,0,1},{1,0,-1,0},{1,1,-2,0},cs,'uvtw')
m = Miller({0,0,0,12},{2,2,-4,0},{6,-6,0,0},cs)
m = Miller({6,-6,0,0},cs)
ori=[]

for m_num=1:size(h,2)
    rot_mat=get_h_rot(h(m_num));
    % rot_mat=get_h_rot
    ori=[ori rot_mat]
end
% ori = 
%% {2,2,0},{0,0,4}, % Bands we want to plot
InputData.IPlot = I
figure
for j=1:size(ori,2)
        r=inv(ori(j).matrix);
        % euler=ori(j).Euler();
        pattern=EBSP_gen(EBSD_geometry,r,screen_int,0);
        fname=fullfile(export_path, sprintf( '%s_gnomonic_%s.png' ,InputUser.Phase_Input{i}, char(h(j))));
        imwrite(rescale(flipud(pattern)),fname)
        % title(sprintf('Rz = %.0f °, Rx =%.0f°, Rz2 = %.0f°',angles(j),angles(i),0))
        label_bands(pattern,EBSD_geometry,m,inv(rotation.byMatrix(r)),InputData)
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_with_bands.png' ,InputUser.Phase_Input{i}, char(h(j))));
        exportgraphics(gcf,fname2)
        % close all
        % error('end')
        %
end
%%
distances = [-6 0]
dSize = size(distances,2);
rectifiedBands = extractBandsDiscrete(pattern,EBSD_geometry,m,r,InputData,'resolution',0.1,'kikuchiWidth',5);
rectifiedBandsExtracted = getBandProfileZoneDist(zoneAxis,distances*degree,rectifiedBands);
%%

figure
t1 = tiledlayout(1,2,'TileSpacing','Compact');
t2 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3 = tiledlayout(t1,'flow','TileSpacing','Compact');
t4 = tiledlayout(t3,'flow','TileSpacing','Compact');
nexttile(t2)

% subplot(1,2,1)
imagesc( EBSD_geometry.x_screen,  EBSD_geometry.y_screen,pattern)
axis xy;
axis equal;
hold on
colormap gray;
for iii = 1:size(rectifiedBandsExtracted,2)
for jjj = 1:size(rectifiedBandsExtracted(iii).bandProfile,2)
bandDist = rectifiedBandsExtracted(iii).bandProfile(jjj).zoneAxisPositonDist;
gnCoords = rectifiedBandsExtracted(iii).gn_coords(bandDist,:,:);
% gnSize = size(gnCoords,1)
for kkk = 1: size(gnCoords,1)
    plot(gnCoords(kkk,:,1),gnCoords(kkk,:,2),'Color',InputData.cset(iii,:)*(1-0.3*(dSize-kkk+1)+0.3),'linewidth',3)
%     disp('plotting')
end
end
end
t3.Layout.Tile = 2;
% nexttile(t3)
hold on

for iii = 1:size(rectifiedBandsExtracted,2)
for jjj = 1:size(rectifiedBandsExtracted(iii).bandProfile,2)
    band = rectifiedBandsExtracted(iii).bandProfile(jjj).band;
    if ~isempty(band)
        nexttile(t4)
        hold on
        bandDist = rectifiedBandsExtracted(iii).bandProfile(jjj).zoneAxisPositonDist;
        gnCoords = rectifiedBandsExtracted(iii).gn_coords(bandDist,:,:);
        
        for kkk = 1: size(gnCoords,1)
%             band = rectifiedBandsExtracted(iii).bandProfile(jjj).band
            plot(rectifiedBandsExtracted(iii).omega/degree,band(kkk,:),'Color',InputData.cset(iii,:)*(1-0.3*(dSize-kkk+1)+0.3),'linewidth',3,'DisplayName', sprintf("%.1f °", distances(kkk)))
        
%         disp('plotting')
        end
        pN = rectifiedBandsExtracted(iii).planeNormal;
        zA = rectifiedBandsExtracted(iii).bandProfile(jjj).zAxis;
        disp(zA)
        title(sprintf('%s @ %s',pN, zA))
        t4.Padding = 'tight';
        t4.TileSpacing = 'tight';
        legend()
    end
end
end
%%
% figure
% plot(rectifiedBands(