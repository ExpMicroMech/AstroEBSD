clear; close all;
InputUser.MTEX_loc='C:\Users\Lukas\mtex-5.10.2\mtex-5.10.2';
try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end

InputUser.Phase_Folder="C:\Users\Lukas\AstroEBSD_v2\phases";
InputUser.Astro_loc='C:\Users\Lukas\AstroEBSD_v2';
run([InputUser.Astro_loc,'\start_AstroEBSD']);
%% load input analysis
% input_folder='D:\Data_Old\Lukas\NbNi\Nb-49Ni Specimen 1 Site 3 Map Data 4ph5o_EBSD\vt0_2'
input_folder='D:\Data_Old\Lukas\Nb-44Ni\h5oina\Nb-44Ni Specimen 1 Site 6 Map Data 5ph5o_EBSD\vt0,2'
% input_folder='D:\Data_Old\Lukas\NbCo\Nb52Co\20240227'
input_file='Results.mat'
input_pats='extracted_patterns.mat'
load(fullfile(input_folder,input_file));
load(fullfile(input_folder,input_pats));
%
[ MapInfo.MapData,MicroscopeData,~,~ ] = breadh5oina(InputUser); %change input user, if the file moved
% MapInfo.MapData.g_num=MapInfo.g_num;
%%
if isfield(MapInfo,'g_num')
    MapInfo.MapData.g_num=MapInfo.g_num
end
extractions=weight_oris_pca_from_map(tile,MapInfo.MapData,extractions);
%% load extracted patterns, in case analysis was not run at full resolution
% load('D:\Data_Old\Lukas\NbCo\Nb52Co\20240220\Nb52Co_1\Nb52Co_1 Specimen 1 Site 9 Map Data 9ph5o_EBSD\\vt1_extracted_256\extracted_patterns.mat');
%%

names = fieldnames(tile);
match = {};
for i = 1:length(names)
  if ~isempty(regexp(names{i},'Tile\d*'))
      match{end+1} = names{i};
  end
end
t_size=size(match,2)
%%
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
RTM_info=RTM;
%run the analysis code
p_cur=1
EBSP_One.PatternIn=extractions.(match{1}).patterns(:,:,p_cur);


%% start the hough based refinement Input parameters in analogy to input deck single, but we can take some rom the loaded data

%%0.54,0.105,0.608]
%PC search ranges

InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path
InputUser.FolderOut=InputUser.SavingFolder
InputUser.FileOut='temp'
InputUser.Phase_Input={'Nb6Ni7_bruker_more_reflex'}
InputUser.Phases={'Nb6Ni7_bruker_more_reflex'}
InputUser.Mode = 'Isolated'
InputUser.PCSearch=0;
InputUser.PatternFlip=0;

pcinx=mean(MapInfo.MapData.PCX)
pciny=mean(MapInfo.MapData.PCY)
pcindd=mean(MapInfo.MapData.DD)
Settings_PCin.start=[pcinx pciny pcindd]; 
PC_start=Settings_PCin.start
% Settings_PCin.start=[0.5261,0.2990,0.7089]%% would be better to read it from the ebsd file directly. but somehow this inormation got lost in th epcae. check to export this

%% correction settings
Settings_Cor_rtm=Settings_Cor;
Settings_Cor.gfilt=1;
Settings_Cor.gfilt_s=4;
Settings_Cor.radius=0; %use a radius mask
Settings_Cor.radius_frac=0.65;
Settings_Cor.resize=1;
Settings_Cor.Square=0; %square crop (1)
Settings_Cor.SquareCrop=1;
Settings_Cor.MeanCentre=0;

pattern_info.size=size(EBSP_One.PatternIn)/4;
pattern_info.size=uint32(pattern_info.size);
Settings_Cor.size=pattern_info.size;
% Settings_Cor.Square
% S
Refine.ss=0.02; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=15; %number of interations
Refine.debug=1;



%% colors
 [cmap_div]=cbrewer('div','RdBu',30);
cmap_div(cmap_div<0)=0;
%%
InputUser.DeckPath=mfilename('fullpath');
%obtain the deck path

tile_name='Tile1'
p_cur=1
% EBSP_single.Pattern
Geometry.PC=Settings_PCin.start;
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );

%% try to find the bit from toms refienment using rtm. overwrite some settings, as they are useful for ebsd but not for rtm
pc_rtm_in=Settings_PCin.start;
PC_best=pc_rtm_in;
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
RTM_info.screensize =double(Settings_Cor.size(1)); %size of the library patterns and the resize of the raw EBSPs
RTM_info.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM_info.iterations = 6;%Set the number of iterations to do in the refinement step
RTM_info.LPTsize = 500;
RTM_info.isHex=1;
RTM_info.Rx=RTM.Rx;
RTM_info.Rz=RTM.Rz;
RTM_info.Ry=RTM.Ry;
RTM_setup=RTM_info %% at some point i got lost which info needs to be in which one
RTM_setup.InputMode=InputUser.Mode;
% RTM_setup.issorted=InputUser.issorted;
RTM_setup.PatternFlip=0;
RTM_setup.apply_gauss=0;
% RTM_setup.screensize=RTM_info.screensize;

Detector_tilt = RTM.Rx(MicroscopeData.SampleTilt);
detector_Euler=MicroscopeData.CameraEuler;
g_camera=RTM.Rz(detector_Euler(3))*RTM.Rx(detector_Euler(2))*RTM.Rz(detector_Euler(1));

Settings_Cor_rtm.resize=1  %% we want the full resolution
Settings_Cor.InvSquareCrop=0 %% make a rectangle on top  and bottom
Settings_Cor.SquareCrop=0
Settings_Cor.max_var_pc=0.2;
Settings_Cor.max_var_ori=10*degree;

%%
% detector_Euler(2)/degree-90
%% 
% EBSP_One.PatternIn=flipud(EBSP_One.PatternIn)

for tnum=1:t_size
    %%
    % tnum=1
    %%
    % tnum=1
    tile_name=match{tnum}
    num_comp=size(extractions.(tile_name).patterns,3)
    %%
    for p_cur=1:num_comp
        %%
        p_cur=6
        %% Correction on exp. pattern
        PatternIn=extractions.(tile_name).patterns(:,:,p_cur);
        [ EBSP_single.ebsp_cor,EBSP_single.PatternInfo ] = EBSP_BGCor(PatternIn,Settings_Cor );
        EBSP_single.PatternInfo.size=double(EBSP_single.PatternInfo.size)
        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,pc_rtm_in );
        [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
        [EBSD_simulation ] = EBSP_Gnom( RTM_info,pc_rtm_in); %you can change PC_in if you want
        % if rotdata{1}.error>1.6*degree
        % disp('selecting ori from map')
        %% Setup input orientation matrix
        eangs=extractions.(tile_name).oris(p_cur,:);
        gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
        Detector_tilt = RTM.Rx(MicroscopeData.SampleTilt);
        gmatrix=gmatrix*Detector_tilt*g_camera';
        % [ Pat_sim_eang ] = EBSP_gen( EBSD_simulation,gmatrix,screen_int);
        %%  Check for pseudo symmetry
        Settings_Cor.resize=1;
        Settings_Cor.SquareCrop=0;
        pseudo_rot=RTM.Rz(0)*RTM.Rx(0)*RTM.Rz(60*degree);
        %%
        [gmatrix_pseudo,~]=ftestPseudo(pc_rtm_in,EBSP_single.ebsp_cor,gmatrix,EBSP_single.PatternInfo,Refine,Settings_Cor,screen_int,RTM_info,RTM_setup,pseudo_rot);
        % gmatrix_pseudo=gmatrix;
        %% Refine patterncenter
        Settings_Cor.resize=1;
        Settings_Cor.SquareCrop=0;
        Settings_Cor.gfilt=0;
        % [gmatrix_ref,PC_best1,EBSP_ref,phbest] = fPCRefine_b(pc_rtm_in,EBSP_single.ebsp_cor,gmatrix_pseudo,EBSP_single.PatternInfo,Refine,Settings_Cor,screen_int,RTM_info,RTM_setup,Settings_Cor);
        eangs=conv_G_to_EA(gmatrix_pseudo);
%%
        if ~any(isnan(eangs))
            %%
            [pc_best_lpr,eangs_out_lpr,phbest_lpr,LPR_lpr]=optimize_pc_ori(pc_rtm_in,eangs,PatternIn,Settings_Cor,screen_int,'Mode','fminsearch','LPRatio');
            [pc_best,eangs_out,phbest,LPR]=optimize_pc_ori(pc_rtm_in,eangs,PatternIn,Settings_Cor,screen_int,'Mode','fminsearch');
            %%
        else
            extractions.(tile_name).rotation_data_refined(:,:,p_cur)=NaN;
            extractions.(tile_name).refined_pc(:,p_cur)=NaN;
            extractions.(tile_name).peak_hight(:,p_cur)=NaN;
            extractions.(tile_name).LPR(:,p_cur)=NaN;

            continue;
        end
        % close;
        %%
        gmatrix_ref=conv_EA_to_G(eangs_out);
        gmatrix_ref_lpr=conv_EA_to_G(eangs_out_lpr);
        %% Save the information into a struct for later analysis and generate images for plotting
        %
        EBSP_single.PatternInfo.resize=1;
        EBSP_single.PatternInfo.size=size(extractions.(tile_name).patterns(:,:,p_cur));
        
        Settings_Cor2=Settings_Cor   
        Settings_Cor2.size=EBSP_single.PatternInfo.size;
        Settings_Cor2.radius=0;
        Settings_Cor2.SquareCrop=0;
        Settings_Cor2.resize=0;
        Settigns_Cor2.gfilt=0;%use a radius mask
        % Settings_Cor.radius_frac=0.75;
        
        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,Settings_PCin.start);
        [ Pat_sim_eang] = EBSP_gen( EBSD_geometry,gmatrix,screen_int);
        [ EBSD_geometry_ref ] = EBSP_Gnom( EBSP_single.PatternInfo,pc_best);
        [ Pat_sim_eang_ref ] = EBSP_gen( EBSD_geometry_ref,gmatrix_ref,screen_int,'lattice_scaling',LPR);
        [Pat_sim_eang_ref]=EBSP_BGCor(Pat_sim_eang_ref,Settings_Cor2);

        [ Pat_sim_eang_ref_lpr] = EBSP_gen( EBSD_geometry_ref,gmatrix_ref,screen_int,'lattice_scaling',LPR_lpr);
        [Pat_sim_eang_ref_lpr]=EBSP_BGCor(Pat_sim_eang_ref,Settings_Cor2);
        %%
        PatternIn=EBSP_BGCor(extractions.(tile_name).patterns(:,:,p_cur),Settings_Cor2);
        extractions.(tile_name).rotation_data_refined(:,:,p_cur)=gmatrix_ref;
        extractions.(tile_name).refined_pc(:,p_cur)=pc_best;
        extractions.(tile_name).peak_hight(:,p_cur)=phbest;
        extractions.(tile_name).LPR(:,p_cur)=LPR;
        eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
        eangs_ref_i_lpr=conv_G_to_EA(gmatrix_ref_lpr)/degree;
        eangs_hough_i=conv_G_to_EA(gmatrix)/degree;
        %% Plot and export, so we can have a look at it later.
        f=figure

        subplot(2,3,1)
        pPattern(PatternIn,EBSD_geometry)
        axis off;
        axis equal;
        title('exp')
        % text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_hough_i(1),eangs_hough_i(2),eangs_hough_i(3)),'units','normalized')

        
        subplot(2,3,2)
        pPattern(Pat_sim_eang_ref,EBSD_geometry)
        axis off;
        axis equal;
        title(sprintf('Sim no LPR, %.3f',phbest))
        text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_ref_i(1),eangs_ref_i(2),eangs_ref_i(3)),'units','normalized')

        ax_par0=subplot(2,3,3)
        % pPattern_alt(normalize(PatternIn),EBSD_geometry_ref)
        % axis off;
        % axis equal;
        % title('original,normalization')

        imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,normalize(PatternIn)-normalize(Pat_sim_eang_ref))
        axis xy;
        axis off;
        axis equal;
        
        % 0
        colormap(ax_par0,cmap_div)
        clim([-3,3])
        % colormap(ax_par0,cmap_div)
        title('Sim, no LPR - Exp')


        subplot(2,3,5)

        imagesc(Pat_sim_eang_ref_lpr)
        axis xy%,EBSD_geometry_ref)
        axis off;
        axis equal;
        title(sprintf('Sim with LPR %.3f',phbest_lpr))
        % title()
        text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_ref_i_lpr(1),eangs_ref_i_lpr(2),eangs_ref_i_lpr(3)),'units','normalized')
        text(0, -0.3, sprintf('LPR= %.4f , %.4f, %.4f', LPR_lpr(1),LPR_lpr(2),LPR_lpr(3)),'units','normalized')
        ax_par=subplot(2,3,6)


        imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,normalize(PatternIn)-normalize(Pat_sim_eang_ref_lpr))
        axis xy;
        axis off;
        axis equal;
        title('Ref - Exp, norm.')
        
        colormap(ax_par,cmap_div)
        clim([-3,3])
        colormap(ax_par,cmap_div)
        title('Sim, with LPR - Exp')
        % text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_ref_i_lpr(1),eangs_ref_i_lpr(2),eangs_ref_i_lpr(3)),'units','normalized')
        % % % % ax_par2=subplot(2,3,6)
        % % % % imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,normalize(Pat_sim_eang)-normalize(Pat_sim_eang_ref))
        % % % % axis xy;
        % % % % axis off;
        % % % % axis equal;
        % % % % title('Sim_{hough} - Sim_{ref}')
        % % % 
        % % colormap(ax_par2,cmap_div)
        % clim([-3,3])
        colorbar('position',[0.92,0.16,0.02,0.22])
        % % % colormap(ax_par0,cmap_div)
        exportgraphics(f,fullfile(input_folder,sprintf('Refinement_pca_tile_%i_comp%i_LPR_comparison.png',tnum,p_cur)),'Resolution','600')
        % close all;

        
        % error('end')
    end
end
    %%
output_pats='output_pats_LPR_comparison.mat'
%%
save(fullfile(input_folder,output_pats),"extractions_no_LPR");
%% Space for testing stuff
% function