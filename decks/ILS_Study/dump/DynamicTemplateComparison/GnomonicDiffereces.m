clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
%%
% L. Berners, 2025. This deck can project gnomonic reprojections form
% dynamic templates. You can specify the zone axis that is supposed to be
% in center, and the euler angles are automatically calculated from this.
% 
% There is also code to include multiple phase files that match a certain
% string pattern.
% The deck will also plot th difference between the first phase and all
% subsequent phases.
%
%
%
%
%
%
%
%






InputUser.Phase_Folder=Phase_Folder

% InputUser.Phase_Input={'Nb_em','Nb_em_10kV','Nb_em_5kV','Ni_em_10kV','Ni','Ni_em_5kV','NbNi3_em','NbNi3_em_5kV','NbNi3_em_10kV','Nb6Ni7_em','Nb6Ni7_em_5kV','Nb6Ni7_em_10kV'};
% InputUser.Phase_Input={'Cu','Ni_em20','10Cu90Ni','20Cu80Ni','30Cu70Ni','40Cu60Ni','50Cu50Ni','60Cu40Ni','70Cu30Ni','80Cu20Ni','90Cu10Ni'}
% InputUser.Phase_Input={'Nb7Ni6_em','Nb6Ni7_em'}%,Ni_em20','10Cu90Ni','20Cu80Ni','30Cu70Ni','40Cu60Ni','50Cu50Ni','60Cu40Ni','70Cu30Ni','80Cu20Ni','90Cu10Ni'}
% InputUser.Phase_Input={'Au3Ni','75Au25Ni'}
% InputUser.Phase_Input={'Nb7Ni6_em_lattice_fixed','Nb6Ni7_em_lattice_fixed'}
% InputUser.Phase_Input={'Si_HR_lowd'}
% InputUser.Phase_Input = {'NbNi_dt0.400_c26.7_a4.9_3aNb_c1-8_c2-20'}%,'NbNi_dt0.400_c26.7_a4.9_3aNb_c2-12','NbNi_dt0.400_c26.7_a4.9_3aNb_c2-20'}
% InputUser.Phase_Input = {'La7Ni16_30kV','La15Ni32_30kV','LaNi2_30kV'}
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.5_a4.9_20kV_2025_07_22__10_59_34_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% InputUser.Phase_Input = {'NbNi_dt0.400_c26.9_a4.9_20kV_2025_07_21__16_01_10_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.250_c26.9_a4.9_20kV_2025_07_21__16_01_05_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.250_c26.5_a4.9_20kV_2025_07_21__16_01_00_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000',...
%     'NbNi_dt0.400_c26.5_a4.9_20kV_2025_07_21__16_00_52_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}
% InputUser.Phase_Input = {'Nb518Ni482_joubert_20.0kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.000_npx-1000',...
%                          'Nb518Ni482_joubert_20.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.800_npx-1000',...
%                         'Nb518Ni482_joubert_20.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.300_npx-1000',...
%                         'Nb518Ni482_joubert_19.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.800_npx-1000',...
%                         'Nb518Ni482_joubert_19.6kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.600_npx-1000',...
%                         'Nb518Ni482_joubert_19.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...}
%                         'Nb518Ni482_joubert_18.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.800_npx-1000',...}
%                         'Nb518Ni482_joubert_18.3kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-18.300_npx-1000',...
%                         'Nb518Ni482_joubert_17.8kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-17.800_npx-1000'   
% InputUser.Phase_Input ={'Nb518Ni482_joubert_20_2025_07_30__21_40_49_dmin-0.020_c1-8.000_c2-30.000_ehistmin-19.800_npx-1000',...
%    'Nb518Ni482_joubert_20_2025_07_30__19_17_43_dmin-0.020_c1-6.000_c2-30.000_ehistmin-19.800_npx-1000',...
%    'Nb518Ni482_joubert_20.0kV_dmin-0.020_c1-5.000_c2-10.000_ehistmin-20.000_npx-1000'}
% InputUser.Phase_Input={'Nb569Ni431_joubert_20kV_2025_07_24__08_56_35_dmin-0_c1-5_c2-10_ehistmin-19_npx-1000'}%,'Nb569Ni431_joubert',}
subFolder = "3aMixed"
PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/3aMixed/*.pha'))' % make sure to only copy them here.
% PhasesToRead = dir(fullfile(Phase_Folder,'phasefiles/Y*Ni*_dmin*.pha'))'

InputUser.Phase_Input = cell(size(PhasesToRead));
for i=1:numel(PhasesToRead)
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
end
% InputUser.Phase_Input = {'NbNi_dt0.250_c26.7_a4.9_19_2025_08_04__10_23_33_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.275_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.325_c26.7_a4.9_19_2025_08_04__10_23_35_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.300_c26.7_a4.9_19_2025_08_04__10_23_34_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.350_c26.7_a4.9_19_2025_08_04__10_23_36_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.375_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     'NbNi_dt0.400_c26.7_a4.9_19_2025_08_04__10_23_37_dmin-0.020_c1-5.000_c2-10.000_ehistmin-19.300_npx-1000',...
%     }
% cat(PhasesToRead.name)
% InputUser.Phase_Input={''}
% InputUser.Phase_Input={'Au','Ni_em20','AuNi3','10Au90Ni','20Au80Ni','30Au70Ni','40Au60Ni','50Au50Ni','60Au40Ni','70Au30Ni','80Au20Ni','90Au10Ni'}
root_path= fullfile('I:\002Simulation_EBSD\NbNi\EnergyComparison',subFolder)
% root_path = 'J:\Users\Berners_Lukas\clara\Energies_Si_wafer'
% root_path='E:\RWTH_Aachen\Promo\CIF\ecp'
%%
InputData.cset=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0];
InputData.cset=[InputData.cset;InputData.cset*0.5];
% 
% % <<<<<<< HEAD
% InputData.V_in=29E3; %voltage in eV
% =======
InputData.V_in=20E3; %voltage in eV
% >>>>>>> 4fea30eee8f80939ebd6aff0c8a2e124695ade58
InputData.lattice_param_scale=1E-10; %lattice paramters are provided in this unit, often Angstrom (1E-10);
 
%%
[cmap_div]=cbrewer('div','RdBu',30);
cmap_div(cmap_div<0)=0;
%% Deal with the wavelegnth
[InputData.lambda] = calc_lambda(InputData.V_in);

%scale the wavelength to the same units as the dhkl
InputData.lambda_p=InputData.lambda/InputData.lattice_param_scale; %lambda normalised by the scale factor of the unit cell
patterns=zeros(2048,2048,size(InputUser.Phase_Input,2),3); % foruth number for size of millers to plot
%%
for i=1:size(InputUser.Phase_Input,2)
    %% [Crystal_UCell<,Crystal_Family,screen_int,Family_List,RTM_info]=phase_data(Input_Data);
    % try
    % phase_name = 
    [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( {InputUser.Phase_Input{i}},InputUser.Phase_Folder );
    Family_List=Crystal_Family{1};
    [cs,~]=loadCIF(RTM_info.cif_file); 
    %%
    % catch
    %     warning(sprintf('%s could not be loaded',InputUser.Phase_Input{i}));
    %     continue;
    % end
    % [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},InputUser.Phase_Folder);
    %%
    [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
     % catch
        % warning(sprintf('%s could not be loaded',InputUser.Phase_Input{i}));
        % continue;
    % end
    %%
    export_path=fullfile(root_path,InputUser.Phase_Input{i});
    mkdir(export_path)
    
    phase_name=InputUser.Phase_Input{i};
    phase_name = regexprep(phase_name,'202.*_dmin-','')
    phase_name = regexprep(phase_name,'ehistmin.*','')
    %%
    % phase_name = phase_name{1}(1:end-3)
    startIndex = regexp(phase_name,'_');
    for slashes=1:size(startIndex,2)
        phase_name=insertAfter(phase_name,startIndex(slashes)-2+slashes,'\');
    end 
    screen_int=Cube_Generate(RTM_info.bin_file,0);
    % [x,y,z]=sphere(1000);
    % i_data=Cube_Sample(x,y,z,screen_int);
    %% Setup some parameter
    PatternInfo.screensize=2048;
    PC= [0.5 0.25 0.7]
    EBSD_geometry=EBSP_Gnom(PatternInfo,PC);
    %%
    
    
    
    
    [x,y,z]=sphere(2000);
    i_data=Cube_Sample(x,y,z,screen_int,screen_int.isHex);
    i_data = i_data-mean(i_data(:));
    i_data = i_data/std(i_data(:));
    %% Plot the sphere for rotation
    % figure
    % surf(x,y,z,i_data,'EdgeColor','none')
    % colormap gray;
    % axis equal;
    % title(phase_name)
    % clim([-4 4])
    % exportgraphics(gcf,fullfile(export_path,sprintf('%s_sphere.png',InputUser.Phase_Input{i}))) 
% end
% close all
% error('end)')
%%
% ori = equispacedSO3Grid(cs,'resolution',15*degree)
% ori = homochoricSO3Grid(cs,cs,'resolution',15*degree)  
% ori=orientation.byEuler(0,0,0,cs);
    %% Plot some different euler angles
    angles=[0,25,45];
    % figure
    % n=1
    % i=1;
    % j=1;
    % r=ori(i).matrix
    cone_params=360;
    %%
    % cs=cryst
    ori=[]
    % h=Miller({0,0,1},{1,1,1},{0,1,1},cs);
    h=Miller({0,0,0,1},{1,1,-2,1},{1,0,-1,0},{-1,-1,2,0},'UVTW',cs)%% h is for the zone axis along which we want to observe the ecpori=[]
for m_num=1:size(h,2)
    rot_mat=get_h_rot(h(m_num));
    % rot_mat=get_h_rot
    ori=[ori rot_mat]
end
    m=Miller({0,0,2},{0,1,1},{2,2,2},{1,1,2},cs);%{2,2,0},{0,0,4}, % Bands we want to plot
    %%
for j=1:size(ori,2)
        
        r=inv(ori(j).matrix);
        % euler=ori(j).Euler();
        pattern=EBSP_gen(EBSD_geometry,r,screen_int,0);
        patterns(:,:,i,j)=pattern;
        fname=fullfile(export_path, sprintf( '%s_gnomonic_%s.png' ,InputUser.Phase_Input{i}, char(h(j))));
        imwrite(rescale(flipud(pattern)),fname)

        % title(sprintf('Rz = %.0f °, Rx =%.0f°, Rz2 = %.0f°',char(h(j).h),char(h(j).k),char(h(j).l)))
        % label_bands(pattern,EBSD_geometry,m,inv(rotation.byMatrix(r)),InputData)
        % 
        % fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_differnce.png' ,InputUser.Phase_Input{i}, char(h(j))));
        % exportgraphics(gcf,fname2)
        % close all
        % error('end')
        %
end
    % for j=1:length(ori)
    %         r=ori(j).matrix;
    %         euler=ori(j).Euler();
    %         pattern=EBSP_gen(EBSD_geometry,r,screen_int,0);
    %         fname=fullfile(export_path, sprintf( '%s_gnomonic_%.1f_%.1f_%.1f_%i.png' ,InputUser.Phase_Input{i}, euler(1),euler(2),euler(3),j));
    %         imwrite(rescale(flipud(pattern)),fname)
    %         % title(sprintf('Rz = %.0f °, Rx =%.0f°, Rz2 = %.0f°',angles(j),angles(i),0))
    %         %
    % end
 % error('end')
 close all
end
%%
for i=2:size(InputUser.Phase_Input,2)
    for j=1:size(ori,2)
        %%
        % phase_name=InputUser.Phase_Input{i};
        % phase_name=InputUser.Phase_Input{i};
        % phase_name = strsplit(phase_name,'kV')
        %
        % phase_name = phase_name{1}
        phase_name=InputUser.Phase_Input{i};
        phase_name = regexprep(phase_name,'202.*_dmin-','')
        phase_name = regexprep(phase_name,'ehistmin.*','')
        startIndex = regexp(phase_name,'_');
        for slashes=1:size(startIndex,2)
            phase_name=insertAfter(phase_name,startIndex(slashes)-2+slashes,'\');
        end 
        %%
        phase_name2 = InputUser.Phase_Input{1}
        phase_name2 = regexprep(phase_name2,'202.*_dmin-','')
        phase_name2 = regexprep(phase_name2,'ehistmin.*','')
        % phase_name2 = strsplit(phase_name2,'kV')
        % phase_name2 = phase_name2{1}
        startIndex = regexp(phase_name2,'_');
        for slashes=1:size(startIndex,2)
            phase_name2=insertAfter(phase_name2,startIndex(slashes)-2+slashes,'\');
        end 

        f=figure
        f.Position=[100 100 1000 1000]
        diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(diff_pat);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('%s \n- %s',phase_name2,phase_name));
        colormap(cmap_div)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        clim([-2,2])
        colorbar()
        exportgraphics(gcf,fname2)
        close all
         f=figure
        f.Position=[100 100 1000 1000]
        diff_pat=patterns(:,:,i,j);
        imagesc(diff_pat);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_blackwhite.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('%s',phase_name));
        colormap gray
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        clim([-3.3,3.3])
        colorbar()
        exportgraphics(gcf,fname2)
        close all
%%

        [cmap_blues]=cbrewer('seq','Blues',30);
        cmap_blues(cmap_blues<0)=0
        f=figure
        f.Position=[100 100 1000 1000]
        diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(diff_pat);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_minus.png' ,InputUser.Phase_Input{i}, char(h(j))));
        % cbrewer()
        title(char(h(j)),sprintf('%s \n- %s',phase_name2,phase_name));
        colormap(cmap_blues)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        clim([-2,0])
        colorbar()
        exportgraphics(gcf,fname2)
        close all


        [cmap_reds]=cbrewer('seq','Reds',30);
        cmap_reds=flipud(cmap_reds)
        cmap_reds(cmap_reds<0)=0
        f=figure
        f.Position=[100 100 1000 1000]
        diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(diff_pat);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_difference_plus.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('%s \n- %s',phase_name2,phase_name));
        colormap(cmap_reds)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        clim([0,2])
        colorbar()
        exportgraphics(gcf,fname2)
        close all
        %%
        [ssimAll,ssimAllMap] = ssim(patterns(:,:,1,j),patterns(:,:,i,j));
        [luminance,luminanceMap] = ssim(patterns(:,:,1,j),patterns(:,:,i,j),'exponents',[1 0 0]);
        [contrast, contrastMap] = ssim(patterns(:,:,1,j),patterns(:,:,i,j),'exponents',[0 1 0] ) ;
        [structureVal, ssimstructureMap] = ssim(patterns(:,:,1,j),patterns(:,:,i,j),'exponents', [0 0 1] ) ;
        % [structureVal, ssimstructureMap] = ssim( patterns(:,:,1,j),patterns(:,:,i,j) );
        %%
        % [cmap_reds]=cbrewer('seq','Reds',30);
        % cmap_reds=flipud(cmap_reds)
        % cmap_reds(cmap_reds<0)=0
        f=figure
        f.Position=[100 100 1000 1000]
        % diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(ssimAllMap);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimAll.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('SSIM : %.4f \n %s \n %s',ssimAll,phase_name2,phase_name));
        colormap(cmap_reds)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        % clim([0,2])
        colorbar()
        exportgraphics(gcf,fname2)
        % close all
        %%
        % [cmap_reds]=cbrewer('seq','Reds',30);
        % cmap_reds=flipud(cmap_reds)
        % cmap_reds(cmap_reds<0)=0
        f=figure
        f.Position=[100 100 1000 1000]
        % diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(luminanceMap);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimLuminance.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('SSIM Luminance : %.4f \n %s \n %s',luminance,phase_name2,phase_name));
        colormap(cmap_reds)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        % clim([0,2])
        colorbar()
        exportgraphics(gcf,fname2)
        % close all
        
        f=figure
        f.Position=[100 100 1000 1000]
        % diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(contrastMap);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimContrast.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('SSIM Contrast: %.4f \n %s \n %s',contrast,phase_name2,phase_name));
        colormap(cmap_reds)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        % clim([0,2])
        colorbar()
        exportgraphics(gcf,fname2)


        f=figure
        f.Position=[100 100 1000 1000]
        % diff_pat=patterns(:,:,1,j)-patterns(:,:,i,j);
        imagesc(ssimstructureMap);
        fname2=fullfile(export_path, sprintf( '%s_gnomonic_%s_ssimStructure.png' ,InputUser.Phase_Input{i}, char(h(j))));

        title(char(h(j)),sprintf('SSIM structure: %.4f \n %s \n %s',structureVal,phase_name2,phase_name));
        colormap(cmap_reds)
        axis equal
        axis xy
        xlim([0,2048])
        ylim([0,2048])
        % clim([0,2])
        colorbar()
        exportgraphics(gcf,fname2)
        % close all

    end
end
%%
        %     figure
        %     f=figure
        %     h_ecp = axes('Units','normalized','Position',[0 0 1 1],'Parent',f);
        %     % h_ecp=imagesc(pattern)
        %     imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern,'Parent',h_ecp)
        %     axis(h_ecp,'equal','xy','tight');
        %     hold(h_ecp,'on');
        %     colormap gray;
        %      cset=InputData.cset;
        %       sf=sf_calc(EBSD_geometry);
        % 
        %      % cset=InputData.cset;
        % % if h_ecp_bands.Value == 1
        % hold on
        %     if ~isempty(Family_List)
        %         for n=1:size(Family_List,1)
        %             HKLs=Family_List{n}(:,2:4);
        % 
        %             for p=1:size(HKLs,1)
        %                 %calculate the bands
        %                 [bands] = Cone_Build(HKLs(p,:),Crystal_UCell{1}.Astar,InputData.lambda_p,r,cone_params);
        %                 %plot them on the EBSP/ECP
        %                 %%
        %                 hold on
        %                 if p==1
        %                     pBand(bands,cset(n,:),h_ecp,sf,'display',HKLs(p,:))
        %                 else
        %                     pBand(bands,cset(n,:),h_ecp,sf)
        %                 end
        %             end
        %         end
        %     end
        % % % end
        % h_ecp.XLim=([EBSD_geometry.x_screen(1),EBSD_geometry.x_screen(end)]);
        % h_ecp.YLim=([EBSD_geometry.y_screen(1),EBSD_geometry.y_screen(end)]);
        % h_ecp.FontSize=8;
        % legend()
    %%

function label_bands(pattern,EBSD_geometry,m,r,InputData)
    theta=(0:1:360)*degree;
    % m_sym=m.symmetrise('unique','antipodal');
    
    figure
    imagesc(EBSD_geometry.x_screen,EBSD_geometry.y_screen,pattern)
    colormap gray
    axis xy;
    axis off;
    axis equal
    hold on
    % fix that we plot too many bands, and transparncy gets gone
    for i=1:size(m,2)
        %%
        m_sym=m(i).symmetrise('unique','antipodal');
        %%
        col_count=1;
        for j=1:size(m_sym,1)
            col_plot=InputData.cset(i,:);
            col_plot=[col_plot 1];
            % if j==1
            %     plot(-100,-99,'Color',col_plot,'LineWidth',1, 'DisplayName',char(m_sym(j,i)));
            % end

            for k=1:3
                %% depending on case, plot band center or band edges
                if k==1
                omega=pi/2;
                elseif k==2
                    omega = pi/2-asin(InputData.lambda_p/(2*m_sym(j).dspacing));
                else
                    omega = pi/2+asin(InputData.lambda_p/(2*m_sym(j).dspacing));
                end
                plane=circels(r*m_sym(j)./norm(r*m_sym(j)),theta,omega);
            
            hold on
            % m_sym(j)
            gn_coords=[plane.x./plane.z, plane.y./plane.z];
            % only select coordinates that are in front
            gn_coords=gn_coords(plane.z>0.001,:);
            % get rid of everything that is bigger than the pattern
            gn_coords=find_limits(gn_coords,EBSD_geometry);
            %plot the bands
            if ~isempty(gn_coords)
                col_plot=InputData.cset(i,:)*(1-(col_count/size(m_sym,1)));
                lw=col_count*2
                col_plot=[col_plot 0.5];
                if k==1
                   plot(gn_coords(:,1),gn_coords(:,2),':','color',col_plot, 'HandleVisibility','off');%,'DisplayName',char(m_sym(j,i)));
                elseif k==2
                      l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'DisplayName',char(m_sym(j)));
                else
                    l=plot(gn_coords(:,1),gn_coords(:,2),'Color',col_plot,'LineWidth',lw, 'HandleVisibility','off');
                end

                % l=patchline()
                % l.lh=[0,0,0,0.5]
            
            end
               hold on
        % a
            end
            if ~isempty(gn_coords)
                col_count=col_count+1;
            end
    %     disp([i j])
        end
    end
    % legend()
        xlim([EBSD_geometry.x_gn_min, EBSD_geometry.x_gn_max])
        ylim([EBSD_geometry.y_gn_min, EBSD_geometry.y_gn_max])
        legend('FontSize',8)
        % axis off
end

    function pBand(bands,csetn,h_ecp,sf,varargin)
    display=0;
        for i=1:size(varargin,1)
            if varargin{i}=='display'
                display=1;
                display_name=varargin{i+1}
            end
        end
        %%
        bands_clean_centre=cBand(bands.centre_zp,sf);
        if ~isempty(bands_clean_centre.x)
            if display==1
                plot(bands_clean_centre.x,bands_clean_centre.y,'color',csetn,'LineWidth',1,'LineStyle','-','Parent',h_ecp,'DisplayName', sprintf("(%.f %.f %.f)",display_name(1),display_name(2),display_name(3)))
            else
                plot(bands_clean_centre.x,bands_clean_centre.y,'color',csetn,'LineWidth',1,'LineStyle','-','Parent',h_ecp)
            end

            end
        bands_clean_upper=cBand(bands.upper_zp,sf);
        if ~isempty(bands_clean_upper.x)
            plot(bands_clean_upper.x,bands_clean_upper.y,'color',csetn,'LineWidth',0.5,'LineStyle',':','Parent',h_ecp)
        end

        bands_clean_lower=cBand(bands.lower_zp,sf);
        if ~isempty(bands_clean_lower.x)
            plot(bands_clean_lower.x,bands_clean_lower.y,'color',csetn,'LineWidth',0.5,'LineStyle',':','Parent',h_ecp)
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
function gn_coords=find_limits(gn_coords,EBSD_geometry)
%%
    y_width=EBSD_geometry.y_gn_max-EBSD_geometry.y_gn_min;
    x_width=EBSD_geometry.x_gn_max-EBSD_geometry.x_gn_min;
    gn_coords(abs(gn_coords)<1e-4)=0;
    condition=gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.2) & gn_coords(:,2)<(EBSD_geometry.y_gn_max+y_width*0.2) & gn_coords(:,1)>(EBSD_geometry.x_gn_min-x_width*0.2) & gn_coords(:,1)<(EBSD_geometry.x_gn_max+x_width*0.2);
    %%
    gn_coords=gn_coords(condition,:);
%     gn_coords=gn_coords(find(gn_coords(:,2)>(EBSD_geometry.y_gn_min-y_width*0.1) & gn_coords(:,2)<(EBSD_geometry.y_gn_min+y_width*0.1)));
   %% 
end
function [lambda_out] = calc_lambda(V_in)
%CALC_LAMBDA calculate the wavelength of a particle, with relativistic correction
%
% V_in = voltage in eV
%
% lambda_out = output wavelength in m
%
% e.g. 
% [lambda_out] = calc_lambda(V_in)
% where:
% V_in=20E3; %voltage of the microscrope in eV

%constants
h=6.626E-34; %Planck constant in Js
me=9.1E-31; %mass of electron in kg
e_charge=1.602E-19; %charge of an electron in C
c=2.998E8; %speed of light in m/s

%equation with relativistic correction
lambda_out = ( h / sqrt(2 * me * e_charge * V_in) ) / (sqrt(1 + e_charge * V_in /2/me/c^2)); % in m

end

 function bands_clean=cBand(bzp,sf)
        bands_clean.x=bzp(1,:)./bzp(3,:);
        bands_clean.y=bzp(2,:)./bzp(3,:);

        bands_clean.y=bands_clean.y(bands_clean.x>=sf.xmin & bands_clean.x<=sf.xmax);
        bands_clean.x=bands_clean.x(bands_clean.x>=sf.xmin & bands_clean.x<=sf.xmax);

        bands_clean.x=bands_clean.x(bands_clean.y>=sf.ymin & bands_clean.y<=sf.ymax);
        bands_clean.y=bands_clean.y(bands_clean.y>=sf.ymin & bands_clean.y<=sf.ymax);
 end

%  function rot_mat=get_h_rot(h)
% % nach h.j. bunge 1988
%     m=h.xyz
%     h=m(1)
%     k=m(2)
%     l=m(3)
%     Phi= acos(l/sqrt(h^2+k^2+l^2));
%     %%;
%     phi2=acos(k/sqrt(h^2+k^2));
%     %%
% %     rot_x=rotation.byAxisAngle(vector3d.X,angle(vector3d.Z,h(1)));
% %     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% % %     rot_x=rotation.byAxisAngle(vector3d.X,55*degree);
% %     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% %     rot_mat=rot_x*rot_z;
%     if isnan(phi2)
%         rot_mat=rotation.byEuler(0,Phi,0);
%     else
%         rot_mat=rotation.byEuler(0,Phi,phi2);
%     end
% 
%  end

 function c= circels(normals,theta,omega)
%     normals=normals(:);
%     c_all=zeros([length(normals),length(theta),length(omega)])
c_all=[];
    for i = 1:length(normals)
%%
      n = normals.subSet(i);

      % generate circles
      c = axis2quat(n,theta) * axis2quat(orth(n),omega)*n;
%       c=c';
%       c_all=[c_all c];
      % plot circles
%       h = [h,line(c,varargin{:},'hold','doNotDraw')]; %#ok<AGROW>

    end
    %%
end