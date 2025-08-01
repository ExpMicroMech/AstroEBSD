function plot_refinementHQ(PatternIn, Settings_Cor, PatternInfo, screen_int, pc_in, pc_best, eangs, eangs_best, Input_Data)
    c_num = 255;
    [cmap_div] = cbrewer('div', 'RdBu', 255, 'spline');
    center_color = ceil(c_num / 2);
    cmap_div(center_color, :) = [1, 1, 1];
    cmap_div(cmap_div < 0) = 0;
    
    gmatrix = conv_EA_to_G(eangs);
    gmatrix_ref = conv_EA_to_G(eangs_best);

    %% Generate images for plotting
    [EBSD_geometry] = EBSP_Gnom(PatternInfo, pc_in);
    [Pat_sim_eang] = EBSP_gen(EBSD_geometry, gmatrix, screen_int);
    [Pat_sim_eang] = EBSP_BGCor(Pat_sim_eang, Settings_Cor);

    [EBSD_geometry_ref] = EBSP_Gnom(PatternInfo, pc_best);
    [Pat_sim_eang_ref] = EBSP_gen(EBSD_geometry_ref, gmatrix_ref, screen_int);
    [Pat_sim_eang_ref] = EBSP_BGCor(Pat_sim_eang_ref, Settings_Cor);

    eangs_ref_i = conv_G_to_EA(gmatrix_ref) / degree;
    eangs_hough_i = conv_G_to_EA(gmatrix) / degree;

    [Pat_inBG] = EBSP_BGCor(PatternIn, Settings_Cor);
    %% Plot and export, so we can have a look at it later.
    f = figure('Position', [20, 20, 1200, 800]);

    subplot(2, 3, 1);
    pPattern(Pat_sim_eang, EBSD_geometry);
    axis off;
    axis equal;
    title('Sim_{initial}');
    text(0, -0.1, sprintf('Before Refinement: phi1= %.3f, Phi= %.3f, phi2= %.3f, PCx= %.1f, PCy= %.1f, DD=%.3f %s', ...
        eangs_hough_i(1), eangs_hough_i(2), eangs_hough_i(3), pc_in(1), pc_in(2), pc_in(3)), 'units', 'normalized');
    hold on;
    % scatter(0, 0, 100, [0 1 1], 'pentagram', 'filled');
    scatter(0, 0, 100, [0 1 1],'wo', 'LineWidth', 1, 'MarkerEdgeColor', 'k');  scatter(0, 0, 100, [0 1 1], 'kx');
    subplot(2, 3, 2);
    pPattern(Pat_inBG, EBSD_geometry);
    axis off;
    axis equal;
    hold on;
    % scatter(0, 0, 100, [0 1 1], 'pentagram', 'filled');
    scatter(0, 0, 100, [0 1 1],'wo', 'LineWidth', 1, 'MarkerEdgeColor', 'k');  scatter(0, 0, 100, [0 1 1], 'kx');
    title('Experiment');
    colormap(subplot(2,3,3), cmap_div);

    ax_par0 = subplot(2, 3, 3);
    diff_img = normalize(PatternIn(:)) - normalize(Pat_sim_eang(:));
    diff_img = reshape(diff_img, size(PatternIn));
    imagesc(EBSD_geometry.x_screen, EBSD_geometry.y_screen, normalize_radius(PatternIn, Settings_Cor) - normalize_radius(Pat_sim_eang, Settings_Cor));
    axis xy;
    axis off;
    axis equal;
    colormap(ax_par0, cmap_div);
    clim([-3, 3]);
    title('Sim_{initial} - Experiment');

    subplot(2, 3, 4);
    pPattern(Pat_sim_eang_ref, EBSD_geometry_ref);
    axis off;
    axis equal;
    title('Sim_{refined}');
    text(0, -0.1, sprintf('After Refinement: phi1= %.3f, Phi= %.3f, phi2= %.3f, PCx= %.1f, PCy= %.1f, DD=%.3f %s', ...
        eangs_ref_i(1), eangs_ref_i(2), eangs_ref_i(3), pc_best(1), pc_best(2), pc_best(3)), 'units', 'normalized');
    hold on;
    % scatter(0, 0, 100, [0 1 1], 'pentagram', 'filled');
    scatter(0, 0, 100, [0 1 1],'wo', 'LineWidth', 1, 'MarkerEdgeColor', 'k');  scatter(0, 0, 100, [0 1 1], 'kx');
    ax_par = subplot(2, 3, 5);

    diff_img = normalize(PatternIn(:)) - normalize(Pat_sim_eang_ref(:));
    diff_img = reshape(diff_img, size(PatternIn));
    imagesc(EBSD_geometry_ref.x_screen, EBSD_geometry_ref.y_screen, normalize_radius(PatternIn, Settings_Cor) - normalize_radius(Pat_sim_eang_ref, Settings_Cor));
    axis xy;
    axis off;
    axis equal;
    title('Sim_{refined} - Experiment, norm.');

    colormap(ax_par, cmap_div);
    clim([-3, 3]);
    ax_par2 = subplot(2, 3, 6);
    diff_img = normalize(Pat_sim_eang(:)) - normalize(Pat_sim_eang_ref(:));
    diff_img = reshape(diff_img, size(PatternIn));
    imagesc(EBSD_geometry_ref.x_screen, EBSD_geometry_ref.y_screen, normalize_radius(Pat_sim_eang, Settings_Cor) - normalize_radius(Pat_sim_eang_ref, Settings_Cor));
    axis xy;
    axis off;
    axis equal;
    title('Sim_{initial} - Sim_{refined}');

    colormap(ax_par2, cmap_div);
    clim([-3, 3]);
    colorbar('position', [0.92, 0.16, 0.02, 0.22]);
    
    output_folder = [Input_Data.image_folder, '\refinement_results'];

    % Create the output folder if it doesn't exist
    if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    end
    plot_name=fullfile(output_folder,Input_Data.image_name);
    
    exportgraphics(f, plot_name, 'Resolution', '600');
    fprintf('Plot saved to : %s', Input_Data.image_name );
    fprintf('\n');




    %% Save data to text file (working)
    % output_file = fullfile(Input_Data.image_folder, 'refinement_results.txt');
    % fid = fopen(output_file, 'a');
    % fprintf(fid, '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
    %     Input_Data.image_name, ...
    %     eangs_ref_i(1), eangs_ref_i(2), eangs_ref_i(3), ...
    %     pc_best(1), pc_best(2), pc_best(3));
    % fclose(fid);

output_file = fullfile(output_folder, 'refinement_results.txt'); 
fid = fopen(output_file, 'a');

% Remove file extension from the image name
[~, file_name_without_ext, ~] = fileparts(Input_Data.image_name);

fprintf(fid, '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
    file_name_without_ext, ... % Use file name without extension
    eangs_ref_i(1), eangs_ref_i(2), eangs_ref_i(3), ...
    pc_best(1), pc_best(2), pc_best(3));
fclose(fid);



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
    % EBSP2(r_grid<r_thresh)=  normalize(EBSP2(r_grid<r_thresh));
    EBSP2(r_grid<r_thresh)=  (EBSP2(r_grid<r_thresh)); %normaliztion turned off, TBB 2025-08-01


end


end
