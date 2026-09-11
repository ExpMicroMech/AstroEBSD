function [G_RefinedPseudo, regout_RPseudo,pC_Out] = testPseudoSymmetryReal(rtmConfigRefine,MultiPhaseStruct,Pat_Ref_r,...
                            PC, G_SO3,optimisationOptions,varargin)
    
    % numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
    %FOR DEBUGGING:
    % rotPseudo =cat(3,rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree)) 
    
    % rtmConfigRefine.options.rotPseudo = cat(3,rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(00*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree))%, ...
     % ..   rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(180*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(180*degree)*rtmConfigRefine.RTM.Rz(90*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(270*degree)*rtmConfigRefine.RTM.Rz(90*degree)) 
 
    rotPseudo = rtmConfigRefine.options.rotPseudo;
    

    numPseudo = size(rtmConfigRefine.options.rotPseudo,3);
    G_RefinedPseudo = zeros(3,3,numPseudo);
    regout_RPseudo = zeros(numPseudo,1);
    pC_Out = zeros(numPseudo,3);
    
    % this one is fast, so lets do it again here
    [ EBSD_geom ] = EBSP_Gnom( rtmConfigRefine.Settings_Cor,...
                            PC);
    optimisationOptionsCopy = optimisationOptions;
    optimisationOptionsCopy.refinePC =0;
    optimisationOptionsCopy.maxIter =5; %just try to be quick for the symmetry test.
%%
% rotPseudo = cat(3,rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(30*degree),...
                % rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefisne.RTM.Rz(0*degree))
%
% close all
    for pseudoNum = 1:numPseudo
        % pseudoNum
        % G_in = G_SO3;
        % G_SO3_temp = rotPseudo(:,:,1)
        % G_in = rotPseudo(:,:,pseudoNum)*G_SO3_temp;
        G_in = rotPseudo(:,:,pseudoNum)*G_SO3; % From the front, as in EBSP_gen, we still transpose the rotation matrix (eq. inversion for orthonromal matrices) (rot*rot2)' = rot2'*rot'
        % [G_RefinedP, regout_RP] = refine5(Pat_Ref_r, EBSD_geom,...
        % % PC, G_in,...
        % SettingsXCF, MultiPhaseStruct.screen_int,...
        % MultiPhaseStruct.screen_int.isHex,...
        % rtmConfigRefine.RTM);
        [ Pat_sim ] = EBSP_gen(EBSD_geom, G_in,MultiPhaseStruct.screen_int);
        % % Pat_sim = EBSP_BGCor(Pat_sim,rtmConfigRefine.Settings_Cor);

        % test wether correct
        eangs = conv_G_to_EA(G_in);
          [pc_best, eangs_out, NCCBest,~,~,optiparams]=optimize_pc_ori_rtm(PC,...
                    eangs,Pat_Ref_r,rtmConfigRefine.Settings_Cor,...
                    MultiPhaseStruct.screen_int,optimisationOptions,varargin{:});
    
        G_RefinedPseudo(:,:,pseudoNum) = conv_EA_to_G(eangs_out);
        regout_RPseudo(pseudoNum) =NCCBest;
        % getMatchMetric(Pat_Ref_r,Pat_sim,'NCC') ;
        
        pC_Out(pseudoNum,:) = pc_best;
        [ Pat_sim ] = EBSP_gen(EBSD_geom, G_in,MultiPhaseStruct.screen_int);
        % 
        % figure
        
        % imagesc(Pat_sim)
        % colormap gray
        % title(sprintf('rot%i'),pseudoNum)
        % axis xy 
        % axis equal
        % axis off
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