%% setup mtex and Astro
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
input_folder='D:\Data_Old\Lukas\Haroon_Orientation_refinement\ECP'
% input_folder='D:\Data_Old\Lukas\NbCo\Nb52Co\20240227'
input_files={'2_cropped.tif','5_cropped.tif','4_cropped.tif'}
input_img=flipud(double(imread(fullfile(input_folder,input_files{1}))));
%% specify phases
InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path
% InputUser.FolderOut=InputUser.SavingFolder
InputUser.FileOut='temp'
InputUser.Phase_Input={'Si'}
InputUser.Phases={'Si'}

%% plot the ECP, to make sure we load it corretly
figure
imagesc(input_img)
colormap gray
%% define other parameters (pc, euler, constraints...)
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation



pc_in=[0.5, 0.5, 4];
% euler_in=[];    %% if you define a list of images , define a list of eulers with 2nd dimension corresponding to number of images
eangs_in=[[108.06, 2.40, -60.56];...
        [42.62, 0.34,-60.10];...
    [127.83, 1.55, -112.41]]*degree;

stage_rot= [0 0 0]*degree; %% stage rotation x y z
sX=stage_rot(1);
sY=stage_rot(2);
sZ=stage_rot(3);

g_stage=RTM.Rz(sZ)*RTM.Rx(sX)*RTM.Ry(-sY);

Detector_tilt = RTM.Rz(0)*RTM.Rx(0)*RTM.Rz(0);


pattern_crop_factor=4;
%% constraints for the refinement
% Settings_Cor.max_var_pc=0.2; % you can use this one, if you only want the
% same constraints on all pc coordinates
Settings_Cor.max_var_pc_x=0;
Settings_Cor.max_var_pc_y=0;
Settings_Cor.max_var_pc_z=1;
Settings_Cor.max_var_ori=10*degree;

%%

%% correction settings
% Settings_Cor_rtm=Settings_Cor;
Settings_Cor.gfilt=0;
Settings_Cor.gfilt_s=4;
Settings_Cor.radius=1; %use a radius mask
Settings_Cor.radius_frac=0.72;
Settings_Cor.resize=1;
Settings_Cor.Square=0; %square crop (1)
Settings_Cor.SquareCrop=0;
Settings_Cor.MeanCentre=0;
Settings_Cor.InvSquareCrop=0 %% make a rectangle on top  and bottom
%%setup_rescaling
pattern_info.size=size(input_img)/pattern_crop_factor;
pattern_info.size=uint32(pattern_info.size);
Settings_Cor.size=pattern_info.size;

%% colors
c_num=255
[cmap_div]=cbrewer('div','RdBu',255);
center_color=ceil(c_num/2);
cmap_div(center_color,:)=[1,1,1];
cmap_div(cmap_div<0)=0;
%%
Geometry.PC=pc_in;
% [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );

%% try to find the bit from toms refienment using rtm. overwrite some settings, as they are useful for ebsd but not for rtm
% pc_rtm_in=pc_in;
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM({InputUser.Phase_Input{1}},InputUser.Phase_Folder);
RTM_info.screensize =double(Settings_Cor.size(1)); %size of the library patterns and the resize of the raw EBSPs
RTM_info.isHex=0;
RTM_info.Rx=RTM.Rx;
RTM_info.Rz=RTM.Rz;
RTM_info.Ry=RTM.Ry;
RTM_setup=RTM_info %% at some point i got lost which info needs to be in which one
% RTM_setup.InputMode=InputUser.Mode;
RTM_setup.PatternFlip=0;
RTM_setup.apply_gauss=0;

for p_num=1:size(input_files,2)
    
    input_img=flipud(double(imread(fullfile(input_folder,input_files{p_num}))));
    eangs=eangs_in(p_num,:);
    gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
    gmatrix=gmatrix*Detector_tilt*g_stage;
    %% Correction on exp. pattern
    [PatternIn,PatternInfo]=EBSP_BGCor(input_img,Settings_Cor);
    [ EBSD_geometry ] = EBSP_Gnom( PatternInfo,pc_in );
    [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

    %% Call the refinement
    % [gmatrix_ref,PC_best1,EBSP_ref,phbest] = fPCRefine_b(pc_rtm_in,EBSP_single.ebsp_cor,gmatrix_pseudo,EBSP_single.PatternInfo,Refine,Settings_Cor,screen_int,RTM_info,RTM_setup,Settings_Cor);
    eangs=conv_G_to_EA(gmatrix);
    [pc_best,eangs_out,phbest]=optimize_pc_ori(pc_in,eangs,PatternIn,Settings_Cor,screen_int,'Mode','fmincon');
    gmatrix_ref=conv_EA_to_G(eangs_out);
    
    %% Generate images for plotting
    [ EBSD_geometry ] = EBSP_Gnom( PatternInfo,pc_in);
    [ Pat_sim_eang ] = EBSP_gen(EBSD_geometry,gmatrix,screen_int);
    [ Pat_sim_eang ]=EBSP_BGCor(Pat_sim_eang,Settings_Cor);
    
    [ EBSD_geometry_ref ] = EBSP_Gnom(PatternInfo,pc_best);
    [ Pat_sim_eang_ref ] = EBSP_gen( EBSD_geometry_ref,gmatrix_ref,screen_int);
    [ Pat_sim_eang_ref ]=EBSP_BGCor(Pat_sim_eang_ref,Settings_Cor);
    
    eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
    eangs_hough_i=conv_G_to_EA(gmatrix)/degree;
    
    %% Plot and export, so we can have a look at it later.
    f=figure;
    
    subplot(2,3,1);
    pPattern(Pat_sim_eang,EBSD_geometry);
    axis off;
    axis equal;
    title('Sim, hough')
    text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f, pc= %.3f %.3f %.3f %s',eangs_hough_i(1),eangs_hough_i(2),eangs_hough_i(3),pc_in(1),pc_in(2),pc_in(3)),'units','normalized')
    hold on
    scatter (0,0,100,[0 1 1],'pentagram','filled');
    
    subplot(2,3,2);
    pPattern(PatternIn,EBSD_geometry)
    axis off;
    axis equal;
    hold on
    scatter (0,0,100,[0 1 1],'pentagram','filled')
    title('original')
    
    ax_par0=subplot(2,3,3);
    % pPattern_alt(normalize(PatternIn),EBSD_geometry_ref)
    % axis off;
    % axis equal;
    % title('original,normalization')
    diff_img=normalize(PatternIn(:))-normalize(Pat_sim_eang(:));
    diff_img=reshape(diff_img,size(PatternIn));
    imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,normalize_radius(PatternIn,Settings_Cor)-normalize_radius(Pat_sim_eang,Settings_Cor))
    axis xy;
    axis off;
    axis equal;
    
    % 0
    colormap(ax_par0,cmap_div);
    clim([-3,3]);
    colormap(ax_par0,cmap_div);
    title('Hough - Exp');
    
    
    subplot(2,3,4)
    
    
    pPattern(Pat_sim_eang_ref,EBSD_geometry_ref)
    axis off;
    axis equal;
    title('Sim, refined')
    text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f, pc %.3f %.3f %.3f',eangs_ref_i(1),eangs_ref_i(2),eangs_ref_i(3),pc_best(1),pc_best(2),pc_best(3)),'units','normalized')
    hold on
    scatter (0,0,100,[0 1 1],'pentagram','filled')
    ax_par=subplot(2,3,5);
    
    diff_img=normalize(PatternIn(:))-normalize(Pat_sim_eang_ref(:));
    diff_img=reshape(diff_img,size(PatternIn));
    imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,normalize_radius(PatternIn,Settings_Cor)-normalize_radius(Pat_sim_eang_ref,Settings_Cor));
    axis xy;
    axis off;
    axis equal;
    title('Ref - Exp, norm.')
    
    colormap(ax_par,cmap_div)
    clim([-3,3])
    colormap(ax_par,cmap_div)
    title('Ref - Exp, norm.')
    
    ax_par2=subplot(2,3,6);
    diff_img=normalize(Pat_sim_eang(:))-normalize(Pat_sim_eang_ref(:));
    diff_img=reshape(diff_img,size(PatternIn));
    imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,normalize_radius(Pat_sim_eang,Settings_Cor)-normalize_radius(Pat_sim_eang_ref,Settings_Cor));
    axis xy;
    axis off;
    axis equal;
    title('Sim_{hough} - Sim_{ref}')
    
    colormap(ax_par2,cmap_div)
    clim([-3,3])
    colorbar('position',[0.92,0.16,0.02,0.22])
    colormap(ax_par0,cmap_div)
    if ~isfield(Settings_Cor,'max_var_pc')
    exportgraphics(f,fullfile(input_folder,sprintf('Refinement_ecp_pcx_con%f.2_pcy_con%f.2_pcz_con%f.2_pattern_%s',Settings_Cor.max_var_pc_x,Settings_Cor.max_var_pc_y,Settings_Cor.max_var_pc_z,input_files{p_num})),'Resolution','600');
    else
     exportgraphics(f,fullfile(input_folder,sprintf('Refinement_ecp_pc_con%f.2_pattern_%s',Settings_Cor.max_var_pc,input_files{pnum})),'Resolution','600');
    end
    close all;
    % error('end')
    %     end
% end
end
    %%
% output_pats='output_pats.mat'
%%
% save(fullfile(input_folder,output_pats),"extractions");
%% Space for testing stuff
% function
function EBSP2=normalize_radius(EBSP2,Settings_Cor)
    %% normalizes a radius cropped image based on the cropping radius

    % EBSP2=fix_mean(EBSP2);
    cs=size(EBSP2);
    r_thresh=Settings_Cor.radius_frac*4/3*cs(1)/2;
    
    [xgrid,ygrid]=meshgrid(1:cs(2),1:cs(1));
    xgrid=double(xgrid);
    ygrid=double(ygrid);
    r_grid=sqrt(double(xgrid-size(EBSP2,2)/2).^2+double(ygrid-size(EBSP2,1)/2).^2);
    EBSP2(r_grid>=r_thresh) = 0;
    EBSP2(r_grid<r_thresh)=  normalize(EBSP2(r_grid<r_thresh));

end