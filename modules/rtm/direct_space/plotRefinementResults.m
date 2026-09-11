function plotRefinementResults(InputUser,rtmConfigRefine,structConfig,MultiPhaseStruct,tileNum)
    disp('---------------- running Plotting of refinements from H5 files')
    disp('################################################################')
    disp('#                                                              #')
    disp('#                                                              #')
    disp('#                                                              #')
    disp('################################################################')
    cmap_div = returnCmapDiv;
    Patterns =h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
    PC = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
                    tileNum));
    %% might be used later (phaseId) % TODO use phaseId to distinguish pahse fore refinement
    % phaseId = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/Phase", ...
                    % tileNum)); %
    % G_in = h5read(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/gmatrixCommercialDetector", ...
                    % tileNum)); %
    if ~isfield(InputUser,'inputEulerMode')
        InputUser.inputEulerMode = 0
    end
    if ~isfield(InputUser,'PatMatchMode')
        InputUser.PatMatchMode = 0
    end
    if ~isfield(InputUser,'PatMatchMode')
        InputUser.inputEulerMode = ""
    end
    if ~isfield(InputUser,'refineLatticeParams')
        InputUser.refineLatticeParams =0;
    end
    if ~isfield(InputUser,'averagedOrisAfterRefinement')
         InputUser.averagedOrisAfterRefinement=0;
    end
    if ~isfield(InputUser,'plotHoughAfterRefinement')
         InputUser.plotHoughAfterRefinement=0;
    end

    if strcmp(InputUser.inputEulerMode,'commercial')
        Eangs_in = h5read(InputUser.h5name, sprintf("/RTM/RTMPCA/%i/EangsCommercial",tileNum));
    % end
    else
        Eangs_in = zeros(size(PC));
    end
    for phaseNumber =1:numel(InputUser.Phase_Input)%,2)
       % read refinement results


       if InputUser.refineLatticeParams ==1
            G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefinedLattice/%i/G_Out", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefinedLattice/%i/RefinedNCC", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            eangsOut = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefinedLattice/%i/eangs_out",...
                InputUser.Phase_Input{phaseNumber},tileNum));
            LatticeRatio = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefinedLattice/%i/LatticeRatio", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
             PC = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefinedLattice/%i/PC", ...
                    InputUser.Phase_Input{phaseNumber},tileNum));
       elseif strcmp(InputUser.PatMatchMode,'full')
            G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/G_Out", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/Refined_PH", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            eangsOut = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/eangs_out",...
                InputUser.Phase_Input{phaseNumber},tileNum))';
             % end
       elseif InputUser.fitPC==0
                  % else
             G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/G_Out", ...
                 InputUser.Phase_Input{phaseNumber},tileNum));
            PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/RefinedNCC", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            eangsOut = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/eangs_out",...
                InputUser.Phase_Input{phaseNumber},tileNum));
              PC = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/PC", ...
                    InputUser.Phase_Input{phaseNumber},tileNum));
       elseif strcmp(InputUser.inputEulerMode,'commercial') 
            G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/G_Out", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/Refined_PH", ...
                            InputUser.Phase_Input{phaseNumber},tileNum));
            eangsOut = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/eangs_out",...
                InputUser.Phase_Input{phaseNumber},tileNum));

       end
       %%
       %sprintf("RTM/OriCorrectedRefinement/%i/oriOFF%.3fdeg_pcOFF%.3f",tileNum,sprintf("RTM/OriCorrectedRefinement/%i/oriOFF%.3fdeg_pcOFF%.3f",tileNum,
          if InputUser.averagedOrisAfterRefinement ==1
              G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/OriCorrectedRefinement/%i/oriOFF0.000deg_pcOFF0.000/gmatrixFilderedAllCombined", ...
                        tileNum));
            % PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMRefined/%i/oriOFF0.000deg_pcOFF0.000", ...
            %                 InputUser.Phase_Input{phaseNumber},tileNum));
            eangsOut = h5read(InputUser.h5name,sprintf("/RTM/OriCorrectedRefinement/%i/oriOFF0.000deg_pcOFF0.000/eangs_outAllCombined", ...
                            tileNum));
            eangsOut = squeeze(eangsOut)
            usedPatterns =  h5read(InputUser.h5name,sprintf("/RTM/OriCorrectedRefinement/%i/oriOFF0.000deg_pcOFF0.000/patternNums", ...
                           tileNum));
          end
          % if strcmp(InputUser.PatMatchMode,'full')
          %       G_SO3 = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/G_SO3", ...
          %                   InputUser.Phase_Input{phaseNumber},tileNum));
          %   PH = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMCoarse/%i/Library_PH", ...
          %                   InputUser.Phase_Input{phaseNumber},tileNum));
          %   % eangsOut = h5read(InputUser.h5name,sprintf("/RTM/%s/RTMPCA/%i/eangs_out",...
          %       % InputUser.Phase_Input{phaseNumber},tileNum));

        %

        if isfield(InputUser,'numPatterns')
            numPatterns = InputUser.numPatterns ;
            if strcmp(InputUser.numPatterns,"")
                 numPatterns = size(Patterns,3);
            end

        else
              numPatterns = size(Patterns,3);
        end
        if InputUser.averagedOrisAfterRefinement ==1
            if numPatterns>numel(usedPatterns)
                numPatterns = numel(usedPatterns)
            end
        end

        InputUser.imDir = fullfile(InputUser.ExportPath,sprintf("Tile_%03i",tileNum),InputUser.Phase_Input{phaseNumber},'RefinementResults');
        mkdir(InputUser.imDir);
        %%
        %%
        InputUser.cs = InputUser.CSList{phaseNumber};
         %%
       % InputUser.imDir = fullfile(InputUser.ExportPath,sprintf("Tile_%03i",tileNum),InputUser.Phase_Input{phaseNumber});

        for patternNumber_loop = 1:numPatterns
                 % tic
            % try
            
            if InputUser.averagedOrisAfterRefinement ==1
                patternNumber = usedPatterns(patternNumber_loop);

            else
                patternNumber = patternNumber_loop;
            end
            disp(sprintf("Plotting Refinement Results: Phase: %i  \t Pattern: %i",...
                    phaseNumber, patternNumber))
            PC_in = PC(:,patternNumber);
            GSO3_in=       G_SO3(:,:,patternNumber_loop);

            Pat_exp = Patterns(:,:,patternNumber);
            Pat_exp_cor=EBSP_BGCor(Pat_exp, rtmConfigRefine.Settings_Cor);
    
            [ EBSD_geom_plot ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                     PC(:,patternNumber) );
            if InputUser.refineLatticeParams ==1
                [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_SO3(:,:,patternNumber_loop),...
                        MultiPhaseStruct(phaseNumber).screen_int,'lattice_scaling',LatticeRatio);
            else
                [ Pat_sim ] = EBSP_gen(EBSD_geom_plot,G_SO3(:,:,patternNumber_loop),...
                        MultiPhaseStruct(phaseNumber).screen_int);
            end
            Pat_sim=EBSP_BGCor(Pat_sim, rtmConfigRefine.Settings_Cor);
            PH = getMatchMetric(Pat_sim,Pat_exp_cor,"NCC");
            % end
            colormap(cmap_div);
            %%
            %% setup structure factor bits for plotting the bands
            structConfig.Phase_Folder = InputUser.Phase_Folder;
            [m_unique,f_hkl] = calcStructureFactor({InputUser.Phase_Input{phaseNumber}}, structConfig.lambda_p ,structConfig);
            [m_plot,I] = cleanAndSortFHKL(m_unique,f_hkl,structConfig.cutoff);
            [m_plot,I] =cleanFHKLDuplicates(m_plot,I);
            structConfig.IPlot = I;
            %% Plot the refinement results
            %%
            f = figure;
            f.Position = [100 100 1500 1000];
            f.Units = 'pixels';
            % figure
            h = tiledlayout(2,3,'TileSpacing','compact');

            % ------------ simulated Pattern --------------------
            ax1 = nexttile(1);
            % subplot(2,3,1)
            pPattern(Pat_sim, EBSD_geom_plot);
            axis off
            % axis equal;
            % axis xy;
            % axis off
            title('Simulated Pattern')
            clim([-3.5,3.5])
            ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % figure
            % ax5=subplot(2,3,4)
           colormap(ax1,'gray');

            
 

            % ------------ simulated Pattern  with bands labeled--------------------
            ax4 = nexttile(4);

            imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim)
            label_bands(Pat_exp_cor,EBSD_geom_plot,m_plot, ...
                inv(rotation.byMatrix(G_SO3(:,:,patternNumber_loop))),... 
                structConfig)
                           % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
            axis ('equal','xy', 'off');
            colormap(ax4,'gray');
            ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            title('Simulated Pattern','fontsize',14)
            clim([-3.5,3.5])
            




            % ------------ experimental Pattern --------------------           
            ax2=nexttile(2);
            
            imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
            % label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
            axis ('equal','xy', 'off');
            colormap(ax2,'gray');
            ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            title('Exp. Pattern','fontsize',14)
            clim([-3.5,3.5])
            



            % ------------ experimental Pattern with Band overlayed--------------
        
            ax5 = nexttile(5);
            
            imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
            label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(G_SO3(:,:,patternNumber_loop))),structConfig)
                             % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
            axis ('equal','xy', 'off');
            colormap(ax5,'gray');
            ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            axis off
            title('Exp. Pattern','fontsize',14)
            clim([-3.5,3.5])
            


            % ------------ difference Pattern                          --------------
        
            ax3 = nexttile(3);
            % ax=subplot(2,3,3)%subplot(1,3,3)
                                    % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            axis ('equal','xy', 'off');
            colormap gray;
            ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            clim([-3,3])
            xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            % axis off
            title('Experiment - Simulation')
            colormap(ax3,cmap_div)
            colorbar('position',[0.91,0.63,0.01,0.25])        
            
            
            % set(gcf, 'Position' , [100 000 1500 900])
            % % ------------ Quotient Patterns                          --------------
            % ax6 = nexttile(6)%subplot(1,3,3)
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,normalize(Pat_exp_cor)./normalize(Pat_sim))
            % axis equal;
            % axis xy;
            % colormap gray;
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % clim([0.01,10])
            % caxis(log([0,10]));
            % xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            % axis off
            % title('Division')
            % colormap(ax6,cmap_div)
            % 
            % ------------ refinement results (belonging to sim pattern --------------------

            text_offset = 0.1;
            text_start = -0.07;
            text_left_offset = -0.08;

            text(text_left_offset,text_start, sprintf('Eangs, out: phi1= %.2f, Phi= %.2f, phi2= %.2f',...
                    eangsOut(1,patternNumber_loop)/degree,...
                    eangsOut(2,patternNumber_loop)/degree,...
                    eangsOut(3,patternNumber_loop)/degree), ...
                    'units', 'normalized', 'fontsize',14)
            text(text_left_offset,text_start-text_offset, sprintf('Eangs, in : phi1= %.2f, Phi= %.2f, phi2= %.2f', ....
                    Eangs_in(1,patternNumber_loop)/degree,...
                    Eangs_in(2,patternNumber_loop)/degree,...
                    Eangs_in(3,patternNumber_loop)/degree),...
                    'units','normalized','fontsize',14)

            %TODO make optional
            % text(0, -0.14, sprintf('LPR= %.4f , %.4f, %.4f', LPR(1),LPR(2),LPR(3)),'units','normalized','fontsize',14)
            
            text(text_left_offset,text_start-text_offset*2, sprintf('PC= %.4f , %.4f, %.4f', ...
                     PC(1,patternNumber_loop), ...
                     PC(2,patternNumber_loop), ...
                     PC(3,patternNumber_loop) ), ...
                    'units', 'normalized', 'fontsize', 14);
            text(text_left_offset,text_start-text_offset*3, sprintf('NCC= %.4f', ...
                     PH),...
                    'Units', 'normalized', 'fontsize',14);
            if InputUser.refineLatticeParams ==1
                text(text_left_offset,text_start-text_offset*4, sprintf(' a/a0 = %.4f, b/b0= %.4f,c/c0= %.4f', ...
                     LatticeRatio(1,patternNumber_loop),LatticeRatio(2,patternNumber_loop),LatticeRatio(3,patternNumber_loop) ),...
                    'Units', 'normalized', 'fontsize',14);
            end
            sgtitle(cleanNameLatex(InputUser.Phase_Input{phaseNumber}),'fontSize',18);
            exportgraphics(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults.png',patternNumber,InputUser.Phase_Input{phaseNumber})),'resolution','300');
          
            exportgraphics(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults.pdf',patternNumber,InputUser.Phase_Input{phaseNumber})),'resolution','300');
            saveas(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults.svg',patternNumber,InputUser.Phase_Input{phaseNumber})),"svg");
            
            
            

            %%
            close all
            %%  Plot a Hough transform of the refinement results to spot differences more easily
            %%

            theta = -90:0.2:89;

            f = figure;
            f.Position = [100 100 1900 1500];
            f.Units = 'pixels';
            % figure
            h = tiledlayout(2,3,'TileSpacing','compact');

            % ------------ simulated Pattern --------------------
            ax1 = nexttile(1);
            % subplot(2,3,1)
            % pPattern(Pat_sim, EBSD_geom_plot);

            colormap gray
              % hard coded for now
            [R_sim,xp_sim] = radon(Pat_sim,theta);
            R_sim = R_sim-mean(R_sim(:));
            R_sim = R_sim./std(R_sim(:));
            % [R,xp] = radon(Pat_sim,-90:0.2:89);
            imagesc(theta,xp_sim,R_sim)
            xlabel('\theta (degrees)')
            ylabel('x''')
            % axis off

            % axis equal;
            % axis xy;
            % axis off
            title('Simulated Pattern, radon Transform')
            % clim([-3.5,3.5])
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % figure
            % ax5=subplot(2,3,4)
           colormap(ax1,'gray');

            
 

            % ------------ simulated Pattern  with bands labeled--------------------
            % ax4 = nexttile(4);
            % 
            % % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_sim)
            % % label_bands(Pat_exp_cor,EBSD_geom_plot,m_plot, ...
            %     % inv(rotation.byMatrix(G_SO3(:,:,patternNumber_loop))),... 
            %     % structConfig)
            %                % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
            % axis ('equal','xy', 'off');
            % colormap(ax4,'gray');
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            % title('Simulated Pattern','fontsize',14)
            clim([-3.9,3.9])
            % 




            % ------------ experimental Pattern --------------------           
            ax2=nexttile(2);
            
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
            [R_exp,xp_exp] = radon(Pat_exp_cor,theta);
            % [R,xp] = radon(Pat_sim,-90:0.2:89);
            R_exp = R_exp-mean(R_exp(:));
            R_exp = R_exp./std(R_exp(:));
            imagesc(theta,xp_exp,R_exp)

            % label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(g_refined)),InputData)
            % axis ('equal','xy', 'off');
            colormap(ax2,'gray');
            xlabel('\theta (degrees)')
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            title('Exp. Pattern, radon Transform','fontsize',14)
            clim([-3.9,3.9])
            



            % ------------ experimental Pattern with Band overlayed--------------
        
            % ax5 = nexttile(5);
            % 
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor)
            % label_bands(pattern,EBSD_geom_plot,m_plot,inv(rotation.byMatrix(G_SO3(:,:,patternNumber_loop))),structConfig)
            %                  % plotBandsStructureFactor(EBSD_geometry,m_plot,inv(rotation.byMatrix(r)),InputData)
            % axis ('equal','xy', 'off');
            % colormap(ax5,'gray');
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            % axis off
            % title('Exp. Pattern','fontsize',14)
            % clim([-3.5,3.5])
            


            % ------------ difference Pattern                          --------------
        
            ax3 = nexttile(3);
            % ax=subplot(2,3,3)%subplot(1,3,3)
                                    % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            imagesc(theta,xp_exp,R_exp-R_sim)
            % axis ('equal','xy', 'off');
            % colormap gray;
            xlabel('\theta (degrees)')
            ylabel('x''')
            % axis off
            title('Experiment - Simulation, radon Transform')
            colormap(ax3,cmap_div)
            colorbar('position',[0.94,0.63,0.01,0.25])        
            clim([-3, 3])

            ax4 = nexttile(4);
            % ax=subplot(2,3,3)%subplot(1,3,3)
                                    % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,Pat_exp_cor-Pat_sim)
            imagesc(theta,xp_exp,R_exp-R_sim)
            % axis ('equal','xy', 'off');
            colormap gray;
            xlabel('\theta (degrees)')
            ylabel('x''')
            % axis off
            title('Experiment - Simulation, radon Transform')
            colormap(ax4,cmap_div)
            colorbar('position',[0.34,0.16,0.01,0.25])        
            % clim([-3, 3])
            
            % set(gcf, 'Position' , [100 000 1500 900])
            % % ------------ Quotient Patterns                          --------------
            % ax6 = nexttile(6)%subplot(1,3,3)
            % imagesc(EBSD_geom_plot.x_screen,EBSD_geom_plot.y_screen,normalize(Pat_exp_cor)./normalize(Pat_sim))
            % axis equal;
            % axis xy;
            % colormap gray;
            % ylim([EBSD_geom_plot.y_gn_min,EBSD_geom_plot.y_gn_max])
            % clim([0.01,10])
            % caxis(log([0,10]));
            % xlim([EBSD_geom_plot.x_gn_min,EBSD_geom_plot.x_gn_max])
            % axis off
            % title('Division')
            % colormap(ax6,cmap_div)
            % 
            % ------------ refinement results (belonging to sim pattern --------------------

            % text_offsealeft_offset,text_start, sprintf('Eangs, out: phi1= %.2f, Phi= %.2f, phi2= %.2f',...
            %         eangsOut(1,patternNumber_loop)/degree,...
            %         eangsOut(2,patternNumber_loop)/degree,...
            %         eangsOut(3,patternNumber_loop)/degree), ...
            %         'units', 'normalized', 'fontsize',14)
            % text(text_left_offset,text_start-text_offset, sprintf('Eangs, in : phi1= %.2f, Phi= %.2f, phi2= %.2f', ....
            %         Eangs_in(1,patternNumber_loop)/degree,...
            %         Eangs_in(2,patternNumber_loop)/degree,...
            %         Eangs_in(3,patternNumber_loop)/degree),...
            %         'units','normalized','fontsize',14)

            %TODO make optional
            % text(0, -0.14, sprintf('LPR= %.4f , %.4f, %.4f', LPR(1),LPR(2),LPR(3)),'units','normalized','fontsize',14)
            
            % text(text_left_offset,text_start-text_offset*2, sprintf('PC= %.4f , %.4f, %.4f', ...
            %          PC(1,patternNumber_loop), ...
            %          PC(2,patternNumber_loop), ...
            %          PC(3,patternNumber_loop) ), ...
            %         'units', 'normalized', 'fontsize', 14);
            % text(text_left_offset,text_start-text_offset*3, sprintf('NCC= %.4f', ...
            %          PH),...
            %         'Units', 'normalized', 'fontsize',14);
            % if InputUser.refineLatticeParams ==1
            %     text(text_left_offset,text_start-text_offset*4, sprintf(' a/a0 = %.4f, b/b0= %.4f,c/c0= %.4f', ...
            %          LatticeRatio(1,patternNumber_loop),LatticeRatio(2,patternNumber_loop),LatticeRatio(3,patternNumber_loop) ),...
            %         'Units', 'normalized', 'fontsize',14);
            % end
            sgtitle(cleanNameLatex(InputUser.Phase_Input{phaseNumber}),'fontSize',18);
            exportgraphics(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults_Radon.png',patternNumber,InputUser.Phase_Input{phaseNumber})),'resolution','300');
            exportgraphics(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults_Radon.pdf',patternNumber,InputUser.Phase_Input{phaseNumber})),'resolution','300');
            saveas(f,fullfile(InputUser.imDir,...
                        sprintf('Pattern_%i_%s_RefinementResults_Radon.svg',patternNumber,InputUser.Phase_Input{phaseNumber})),"svg");
            
            
            
            % catch
               %     warning(sprintf('Could not plot Patterns for %i ', patternNumber)) 
               % end
        end
    end
end