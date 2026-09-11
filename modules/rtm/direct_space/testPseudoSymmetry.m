function [G_RefinedPseudo, regout_RPseudo] = testPseudoSymmetry(rtmConfigRefine,MultiPhaseStruct,Pat_Ref_r,...
                            PC, G_SO3, SettingsXCF)
    
    % numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
    %FOR DEBUGGING:
    % rotPseudo =cat(3,rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree)) 
    
    % rtmConfigRefine.options.rotPseudo = cat(3,rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(00*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree))%, ...
     % ..   rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(180*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(180*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(270*degree)*rtmConfigRefine.RTM.Rz(90*degree)) 
 
    rotPseudo = rtmConfigRefine.options.rotPseudo;
    

    numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
    G_RefinedPseudo = zeros(3,3,numPseudo);
    regout_RPseudo = zeros(4,numPseudo);
    
    % this one is fast, so lets do it again here
    [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                            PC);
%%
% rotPseudo = cat(3,rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(30*degree),...
                % rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree))
%
% close all
    for pseudoNum = 1:numPseudo
        % pseudoNum
        % G_in = G_SO3;
        % G_SO3_temp = rotPseudo(:,:,1)
        % G_in = rotPseudo(:,:,pseudoNum)*G_SO3_temp;
        G_in = rotPseudo(:,:,pseudoNum)*G_SO3; % From the front, as in EBSP_gen, we still transpose the rotation matrix (eq. inversion for orthonromal matrices) (rot*rot2)' = rot2'*rot'
        [G_RefinedP, regout_RP] = refine5(Pat_Ref_r, EBSD_geom,...
        PC, G_in,...
        SettingsXCF, MultiPhaseStruct.screen_int,...
        MultiPhaseStruct.screen_int.isHex,...
        rtmConfigRefine.RTM);
        G_RefinedPseudo(:,:,pseudoNum) = G_RefinedP;
        regout_RPseudo(:,pseudoNum) = regout_RP;
        % test wether correct
        % [ Pat_sim ] = EBSP_gen(EBSD_geom, G_RefinedP,MultiPhaseStruct.screen_int);
        [ Pat_sim ] = EBSP_gen(EBSD_geom, G_in,MultiPhaseStruct.screen_int);
        Pat_sim = EBSP_BGCor(Pat_sim,rtmConfigRefine.Settings_Cor);
        % figure
        % imagesc(Pat_sim)
        % % Pat_EBSP_BGcor = EBSP_BGCor(Pat_Ref_r.pat_in,rtmConfigRefine.Settings_Cor);
    % 
    %     figure
    %     imagesc(Pat_sim)
    %     axis equal
    %     axis xy
    %     colormap gray
    %     ncc= normxcorr2(Pat_Ref_r.pat_in,Pat_sim);
    %     max(ncc(:))
    %     title(sprintf('Pseudo matrix %i',pseudoNum))
    % %      figure
    % imagesc(Pat_Ref_r.pat_in-Pat_sim)
    % colorbar()
    % colormap (returnCmapDiv)
    % clim([-3 3])
    % %
        
    end
    % figure
    % imagesc(Pat_Ref_r.pat_in)
    % colormap gray
    % axis equal
    % axis xy
    


end