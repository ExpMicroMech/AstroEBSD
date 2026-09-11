clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
close all;
InputUser.MTEX_loc=MTEX_loc
InputUser.MTEX_loc='C:\Users\berners\mtex-5.7.0'
try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end

InputUser.Phase_Folder="C:\Users\berners\AstroEBSD_v2\phases\";
InputUser.Astro_loc='C:\Users\berners\AstroEBSD_v2\';
run([InputUser.Astro_loc,'\start_AstroEBSD']);
%% load input analysis
input_folder='J:\Users\Berners_Lukas\MA_Dinah\ErNi2-NS_2024_08_08_EBSD\DB_Files\lukas\Cole\Mapping\Ebsd\MA Dinah\ErNi2-NS\Area 3'
input_folder='J:\Users\Berners_Lukas\Fib\EBSD\C02a_vi\c02_vi(v_reprep)\DB_Files\lukas\Cole\Mapping\Ebsd\SFB_samples\C02_as_cast_vi\Area 3\_EBSD\vt0,2'% input_folder='D:\Data_Old\Lukas\Nb-44Ni\h5oina\Nb-44Ni Specimen 1 Site 6 Map Data 5ph5o_EBSD\vt0,2'
% input_folder='D:\Data_Old\Lukas\NbCo\Nb52Co\20240227'
input_file='Results.mat'
input_pats='extracted_patterns.mat'
load(fullfile(input_folder,input_file));
load(fullfile(input_folder,input_pats));
%%
%[ MapInfo.MapData,MicroscopeData,~,~ ] = breadh5oina(InputUser); %change input user, if the file moved
InputUser.Phase_Input={'Mg_20kv_em'}%'ErNi2_super_structure','ErNi2'};%,'Nb57_Ni43_joubert_c1-8_c2-16'}
InputUser.Phases={'Mg_20kv_em'}%ErNi2_super_structure','ErNi2'}
[~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );
% [cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 
[cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
InputUser.CS=cs

[ MapInfo.MapData,MicroscopeData,~,MapInfo.EBSPData ] = bReadOSC(InputUser); 
% MapInfo.MapData.g_num=MapInfo.g_num;
%%
if isfield(MapInfo,'g_num')
    MapInfo.MapData.g_num=MapInfo.g_num
end
%%
extractions=weight_oris_pca_from_map(tile,MapInfo.MapData,extractions);
%% load extracted patterns, in case analysis was not run at full resolution
% load('D:\Data_Old\Lukas\NbCo\Nb52Co\20240220\Nb52Co_1\Nb52Co_1 Specimen 1 Site 9 Map Data 9ph5o_EBSD\\vt1_extracted_256\extracted_patterns.mat');
%%
appendix_string='with_LPR_new_phases_bethe_var'
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
% rot_edax_conversion=RTM.Rz(90*degree)*RTM.Rx(180*degree);

%% start the hough based refinement Input parameters in analogy to input deck single, but we can take some rom the loaded data

%%0.54,0.105,0.608]
%PC search ranges

InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path
InputUser.FolderOut=InputUser.SavingFolder
InputUser.FileOut='temp'
% InputUser.Phase_Input={'ErNi2_super_structure'};%,'ErNi2'};%,'Nb57_Ni43_joubert_c1-8_c2-16'}
% InputUser.Phases={'ErNi2_super_structure'};%,'ErNi2'};%{'Nb57_Ni43_joubert','Nb57_Ni43_joubert_c1-2_c2-8','Nb57_Ni43_joubert_c1-8_c2-16'}
InputUser.Mode = 'Isolated'
InputUser.PCSearch=0;
InputUser.PatternFlip=0;

if isfield(MapInfo.MapData,'PCX')
    pcinx=mean(MapInfo.MapData.PCX)
    pciny=mean(MapInfo.MapData.PCY)
    pcindd=mean(MapInfo.MapData.DD)
else
    pcinx=0.49
    pciny=0.43
    pcindd=0.65
end
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

pattern_info.size=size(EBSP_One.PatternIn)/2;
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
InputUser.Phase_Folder="C:\Users\berners\AstroEBSD_v2\phases\";
tile_name='Tile1'
p_cur=1
% EBSP_single.Pattern
Geometry.PC=Settings_PCin.start;
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num,RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,InputUser.Phase_Folder );

%% try to find the bit from toms refienment using rtm. overwrite some settings, as they are useful for ebsd but not for rtm
pc_rtm_in=Settings_PCin.start;
PC_best=pc_rtm_in;

% RTM_info_copy_fields=RTM_info;
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
% RTM_info=struct()
RTM_info.screensize =double(Settings_Cor.size(1)); %size of the library patterns and the resize of the raw EBSPs
RTM_info.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM_info.iterations = 6;%Set the number of iterations to do in the refinement step
RTM_info.LPTsize = 500;
RTM_info.isHex=1;
RTM_info.Rx=RTM.Rx;
RTM_info.Rz=RTM.Rz;
RTM_info.Ry=RTM.Ry;
% for fn = fieldnames(RTM_info_copy_fields)'
%    RTM_info.(fn{1}) = RTM_info_copy_fields.(fn{1});
% end
RTM_setup=RTM_info %% at some point i got lost which info needs to be in which one
RTM_setup.InputMode=InputUser.Mode;
% RTM_setup.issorted=InputUser.issorted;
RTM_setup.PatternFlip=0;
RTM_setup.apply_gauss=0;

% RTM_setup.screensize=RTM_info.screensize;
InputUser.SampleTilt=70
InputUser.CameraTilt=5.7
MicroscopeData.SampleTilt=InputUser.SampleTilt
MicroscopeData.CameraTilt=InputUser.CameraTilt
MicroscopeData.TotalTilt=((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
g_camera=RTM.Rx(0*degree);
% g_camera=RTM.Rz(detector_Euler(3))*RTM.Rx(detector_Euler(2))*RTM.Rz(detector_Euler(1));


Settings_Cor_rtm.resize=1  %% we want the full resolution
Settings_Cor.InvSquareCrop=0 %% make a rectangle on top  and bottom
Settings_Cor.SquareCrop=0
Settings_Cor.max_var_pc=0.2;
Settings_Cor.max_var_ori=10*degree;
%%
extractions.Phase_List=InputUser.Phase_Input;
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
        Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);
        gmatrix=gmatrix*Detector_tilt*RTM.Ry(180*degree);
        % [ Pat_sim_eang ] = EBSP_gen( EBSD_simulation,gmatrix,screen_int);
        %%  Check for pseudo symmetry
        Settings_Cor.resize=1;
        Settings_Cor.SquareCrop=0;
        pseudo_rot=RTM.Rz(0)*RTM.Rx(0)*RTM.Rz(60*degree);
        %%
        % [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
        % for fn = fieldnames(RTM_info_copy_fields)'
        %    RTM_info.(fn{1}) = RTM_info_copy_fields.(fn{1});
        % end
        % [gmatrix_pseudo,~]=ftestPseudo(pc_rtm_in,EBSP_single.ebsp_cor,gmatrix,EBSP_single.PatternInfo,Refine,Settings_Cor,screen_int,RTM_info,RTM_setup,pseudo_rot);
        gmatrix_pseudo=gmatrix;
        %% Refine patterncenter
        Settings_Cor.resize=1;
        Settings_Cor.SquareCrop=0;
        Settings_Cor.gfilt=0;
        Settings_Cor.radius=1;
        % [gmatrix_ref,PC_best1,EBSP_ref,phbest] = fPCRefine_b(pc_rtm_in,EBSP_single.ebsp_cor,gmatrix_pseudo,EBSP_single.PatternInfo,Refine,Settings_Cor,screen_int,RTM_info,RTM_setup,Settings_Cor);
        eangs=conv_G_to_EA(gmatrix_pseudo);
        if ~any(isnan(eangs))
            for j=1:Phase_Num
            %%
                [ ~,~,~,~,~, RTM_info_phase ] = Phase_Builder_RTM(  {InputUser.Phase_Input{j}},InputUser.Phase_Folder);
                [screen_int] = Cube_Generate(RTM_info_phase.bin_file,RTM_info_phase.isHex);       
                [pc_best,eangs_out,phbest,LPR]=optimize_pc_ori(pc_rtm_in,eangs,PatternIn,Settings_Cor,screen_int,'Mode','fminsearch');%,'LPRatio');
                %%
                %%
                gmatrix_ref=conv_EA_to_G(eangs_out);
                extractions.(tile_name).phase_index(p_cur,j)=j;
                extractions.(tile_name).rotation_data_refined(:,:,p_cur,j)=gmatrix_ref;
                extractions.(tile_name).refined_pc(:,p_cur,j)=pc_best;
                extractions.(tile_name).peak_hight(:,p_cur,j)=phbest;
                extractions.(tile_name).LPR(:,p_cur,j)=LPR;
                eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
                eangs_hough_i=conv_G_to_EA(gmatrix)/degree;
            end
            %%
            [~,phase_ind]=max(extractions.(tile_name).peak_hight(:,p_cur,:))
            extractions.(tile_name).Phase(p_cur)=phase_ind;
            %,'LPRatio')
            %%
        else
            extractions.(tile_name).Phase(p_cur)=0;
             for j=1:Phase_Num
    
                
                extractions.(tile_name).rotation_data_refined(:,:,p_cur,j)=NaN;
                extractions.(tile_name).refined_pc(:,p_cur,j)=NaN;
                extractions.(tile_name).peak_hight(:,p_cur,j)=NaN;
                extractions.(tile_name).LPR(:,p_cur,j)=NaN;
             end
            continue;
        end
        %%
        gmatrix_ref=extractions.(tile_name).rotation_data_refined(:,:,p_cur,phase_ind);
        pc_best=extractions.(tile_name).refined_pc(:,p_cur,phase_ind);
        
        %% Save the information into a struct for later analysis and generate images for plotting
        %
        EBSP_single.PatternInfo.resize=1;
        EBSP_single.PatternInfo.size=size(extractions.(tile_name).patterns(:,:,p_cur));
        
        Settings_Cor2=Settings_Cor   
        Settings_Cor2.size=EBSP_single.PatternInfo.size;
        Settings_Cor2.radius=1;
        Settings_Cor2.SquareCrop=0;
        Settings_Cor2.resize=0;
        Settigns_Cor2.gfilt=0;%use a radius mask
        % Settings_Cor.radius_frac=0.75;


        cur_phase=extractions.(tile_name).Phase(p_cur)
        %%
        eangs_best=conv_G_to_EA(extractions.(tile_name).rotation_data_refined(:,:,p_cur,cur_phase))/degree;
        %%
        gmatrix_ref=extractions.(tile_name).rotation_data_refined(:,:,p_cur,cur_phase)
        %%
        eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
        LPR=extractions.(tile_name).LPR(:,p_cur,cur_phase);
        peak_hight=extractions.(tile_name).peak_hight(:,p_cur,cur_phase);


        [ ~,~,~,~,~, RTM_info_phase ] = Phase_Builder_RTM(  {InputUser.Phase_Input{cur_phase}},InputUser.Phase_Folder);
        [screen_int] = Cube_Generate(RTM_info_phase.bin_file,RTM_info_phase.isHex); 

        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,Settings_PCin.start);
        [ Pat_sim_eang] = EBSP_gen( EBSD_geometry,gmatrix,screen_int);
        Pat_sim_eang=EBSP_BGCor(Pat_sim_eang,Settings_Cor2);
        [ EBSD_geometry_ref ] = EBSP_Gnom( EBSP_single.PatternInfo,pc_best);
        [ Pat_sim_eang_ref ] = EBSP_gen( EBSD_geometry_ref,gmatrix_ref,screen_int,'lattice_scaling',LPR);
        [Pat_sim_eang_ref]=EBSP_BGCor(Pat_sim_eang_ref,Settings_Cor2);
        %%
        PatternIn=EBSP_BGCor(extractions.(tile_name).patterns(:,:,p_cur),Settings_Cor2);
         % for j=1:Phase_Num
         %    extractions.(tile_name).Phase(p_cur,j)=InputUser.Phase_Input{j}};
         %    extractions.(tile_name).rotation_data_refined(:,:,p_cur,j)=gmatrix_ref;
         %    extractions.(tile_name).refined_pc(:,p_cur,j)=pc_best;
         %    extractions.(tile_name).peak_hight(:,p_cur,j)=phbest;
         %    extractions.(tile_name).LPR(:,p_cur,j)=LPR;
         %    eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
         %    eangs_hough_i=conv_G_to_EA(gmatrix)/degree;
         % end
        
        %% Plot and export, so we can have a look at it later.
        % set(groot,'defaultAxesFontSize',12)

        f=figure;
        f.Position=[100 100 900 700]
        subplot(2,3,1)
        pPattern(Pat_sim_eang,EBSD_geometry)
        axis off;
        axis equal;
        title('Sim, hough')
        text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_hough_i(1),eangs_hough_i(2),eangs_hough_i(3)),'units','normalized')

        
        subplot(2,3,2)
        pPattern(PatternIn,EBSD_geometry)
        axis off;
        axis equal;
        title('original','fontsize',14)

        ax_par0=subplot(2,3,3)
        % pPattern_alt(normalize(PatternIn),EBSD_geometry_ref)
        % axis off;
        % axis equal;
        % title('original,normalization')

        imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,PatternIn-Pat_sim_eang)
        axis xy;
        axis off;
        axis equal;
        
        % 0
        colormap(ax_par0,cmap_div)
        clim([-3,3])
        colormap(ax_par0,cmap_div)
        title('Hough - Exp')


        subplot(2,3,4)

        pPattern(Pat_sim_eang_ref,EBSD_geometry_ref)
        axis off;
        axis equal;
        title('Sim, refined')
        text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_ref_i(1),eangs_ref_i(2),eangs_ref_i(3)),'units','normalized')
        text(0, -0.3, sprintf('LPR= %.4f , %.4f, %.4f; PH= %.4f', LPR(1),LPR(2),LPR(3),phbest),'units','normalized')
        text(0, -0.3, sprintf('PC= %.4f , %.4f, %.4f', PC(1),PC(2),PC(3)),'units','normalized')
        ax_par=subplot(2,3,5)


        imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,PatternIn-Pat_sim_eang_ref)
        axis xy;
        axis off;
        axis equal;
        title('Ref - Exp, norm.')
        
        colormap(ax_par,cmap_div)
        clim([-3,3])
        colormap(ax_par,cmap_div)
        title('Ref - Exp, norm.')

        ax_par2=subplot(2,3,6)
        imagesc(EBSD_geometry_ref.x_screen,EBSD_geometry_ref.y_screen,Pat_sim_eang-Pat_sim_eang_ref)
        axis xy;
        axis off;
        axis equal;
        title('Sim_{hough} - Sim_{ref}')
        
        colormap(ax_par2,cmap_div)
        clim([-3,3])
        colorbar('position',[0.92,0.16,0.02,0.22])
        colormap(ax_par0,cmap_div)
        exportgraphics(f,fullfile(input_folder,sprintf('Refinement_pca_tile_%i_comp%i_%s.png',tnum,p_cur,appendix_string)),'Resolution','600')
        % close all;
        %%

        f=figure;
        f.Position=[100 100 300*Phase_Num+400 700];

        subplot(2,Phase_Num+1,1)
        pPattern(PatternIn,EBSD_geometry)
        axis off;
        axis equal;
        title('experimental','fontsize',14)
        text(0,-0.1, sprintf('phi1= %.2f, Phi= %.2f, phi2= %.2f',eangs_best(1),eangs_best(2),eangs_best(3)),'units','normalized')
        phase_name=InputUser.Phase_Input{cur_phase};
        startIndex = regexp(phase_name,'_');
        for slashes=1:size(startIndex,2)
            phase_name=insertAfter(phase_name,startIndex(slashes)-2+slashes,'\');
        end 
        text(0,-0.5, sprintf('selected Phase : \n %s',phase_name),'units','normalized')
        for p_num=1:Phase_Num

            phase_name=InputUser.Phase_Input{p_num};
            startIndex = regexp(phase_name,'_');
            for slashes=1:size(startIndex,2)
                phase_name=insertAfter(phase_name,startIndex(slashes)-2+slashes,'\');
            end 
            %%
            ax_par=subplot(2,Phase_Num+1,p_num+1);
            % eangs_ref_i=conv_G_to_EA(extractions.(tile_name).rotation_data_refined(:,:,p_cur,cur_phase))/degree
            gmatrix_ref=extractions.(tile_name).rotation_data_refined(:,:,p_cur,p_num);
            eangs_ref_i=conv_G_to_EA(gmatrix_ref)/degree;
            LPR=extractions.(tile_name).LPR(:,p_cur,p_num);
            PC=extractions.(tile_name).refined_pc(:,p_cur,p_num);
            peak_hight=extractions.(tile_name).peak_hight(:,p_cur,p_num);

            [ ~,~,~,~,~, RTM_info_phase ] = Phase_Builder_RTM(  {InputUser.Phase_Input{p_num}},InputUser.Phase_Folder);
            [screen_int] = Cube_Generate(RTM_info_phase.bin_file,RTM_info_phase.isHex); 

            [ Pat_sim_eang_ref ] = EBSP_gen( EBSD_geometry_ref,gmatrix_ref,screen_int,'lattice_scaling',LPR);
            [Pat_sim_eang_ref]=EBSP_BGCor(Pat_sim_eang_ref,Settings_Cor2);

            imagesc(Pat_sim_eang_ref)%,EBSD_geometry_ref);
            axis xy;
            axis off;
            axis equal;
            title(phase_name);
            colormap(ax_par,'gray')
            text(0,-0.1, sprintf('eangs= %.2f, %.2f, %.2f',eangs_ref_i(1),eangs_ref_i(2),eangs_ref_i(3)),'units','normalized');
            text(0, -0.18, sprintf('LPR= %.4f , %.4f, %.4f', LPR(1),LPR(2),LPR(3)),'units','normalized');
            text(0, -0.26, sprintf('PC= %.4f , %.4f, %.4f', PC(1),PC(2),PC(3)),'units','normalized');
            
            % ax_par=subplot(2,3,5)


            ax_par0=subplot(2,Phase_Num+1,Phase_Num+1 +   p_num+1);
            imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,PatternIn-Pat_sim_eang_ref);
            axis xy;
            axis off;
            axis equal;
            % colormap(ax_par0,cmap_div)
            clim([-3,3])
            colormap(ax_par0,cmap_div)
            title(sprintf('Exp - Sim, PH= %.3f',peak_hight))
            text(0, -0.1, sprintf('PC= %.4f , %.4f, %.4f;',PC(1),PC(2),PC(3)),'units','normalized');
            if p_num==Phase_Num
                colorbar('position',[0.92,0.17,0.01,0.22])
            end


        end
        exportgraphics(f,fullfile(input_folder,sprintf('Refinement_pca_tile_%i_comp%i_%s_phase_comparison.png',tnum,p_cur,appendix_string)),'Resolution','600')
        close all;





        % error('end')
    end
end
    %%
output_pats=sprintf('output_pats_phase_comparison_%s.mat',appendix_string);
%%
save(fullfile(input_folder,output_pats),"extractions");
%% Space for testing stuff
% function

for j=1:Phase_Num
            %%
                [ ~,~,~,~,~, RTM_info_phase ] = Phase_Builder_RTM(  {InputUser.Phase_Input{j}},InputUser.Phase_Folder);
                [screen_int] = Cube_Generate(RTM_info_phase.bin_file,RTM_info_phase.isHex);       
                % [pc_best,eang%%

end