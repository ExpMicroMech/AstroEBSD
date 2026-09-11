function P = fitPatternCentreProjective(PCX,PCY,PCZ,x,y)
%% Fitting a projective transformation between PCX,PCY,PCZ and xGrid, yGrid vectors
% P = fitPatternCentreProjective(PCX, PCY, PCZ, x, y) computes a projective transformation 
% matrix that relates the input projection cetnre coordinates (PCX, PCY, PCZ) to the specified 
% grid vectors (x, y). This function constructs a minimization matrix and solves for the 
% transformation parameters using linear algebra techniques.
% NOTE: The minimization matrix is only solvable, if at least 5
% combinations of PCX,PCZ,PCZ, x and y are known. For more the linear
% equation is over determined an the best fit is determined by the matlab
% command.
%
% Inputs:
%   - PCX: Pattern Centre X.
%   - PCY: Pattern Centre Y.
%   - PCZ: Pattern Centre Z.
%   - x: Beam Position on Sample x.
%   - y: Beam Position on Sample y.
%
% Outputs:
%   - P: The projective transformation matrix that aligns the projection
%   center with the sample coordinates.
%
% Usage:
%   	PCs = P(1:3,:)*[x;y;ones(size(x))]
% 
%     scale = P(10)*x + P(11)*y +1;
%     PCs = PCs./scale;
%     PCX = PCs(:,1:3);
%     PCY = PCs(:,1:3);   
%     PCZ = PCs(:,1:3);
% %% BUILD a Large LEQ from this.
% Also have a look at.
% Refined Calibration Model for Improving the
% Orientation Precision of Electron Backscatter
% Diffraction Maps Aimo Winkelmann et al. doi:10.3390/ma13122816
% + Estimating Projective Transformation Matrix
%%(Collineation, Homography)
% Zhengyou Zhang Microsoft Research Techical Report MSR-TR-2010-63
% Handbook of Image
% Processing and Computer
% Vision https://doi.org/10.1007/978-3-030-42374-2

    % lets start with only a small subset to make the matrix easier "graspable"
    % numPatsSel = 55;
    numPatsSel = numel(PCX);
    PCX_temp = PCX(1:numPatsSel);
    PCY_temp = PCY(1:numPatsSel);
    PCZ_temp = PCZ(1:numPatsSel);
    PC_temp = [PCX_temp;PCY_temp;PCZ_temp]; % [pc_x,y,z]x n_points (rowsxcols)
    PC_temp = PC_temp; %[pcx1 pcy1 pcz1 pcx2 pcy2 pcz2 ... pcxn pcyn pczn
    x_temp = x(1:numPatsSel);
    y_temp = y(1:numPatsSel);
    z_temp = ones(size(x_temp));
    gridVecs = [x_temp; y_temp; z_temp]; % 
    % gridVecs = gridVecs(:);     %[x1 y1 z1 x2 y2 z2 x3 y3 z3]....
    
    % matrix needs to look sth like this. How do we automate this ?
    % [x_temp(1) y_temp(1) z_temp(1) 0 0 0 0 0 0    -x_temp(1)*PCX(1) - y_temp(1)*PCX(1);...
     % 0 0 0 x_temp(1) y_temp(1) z_temp(1)  0 0 0  -x_temp(1)*PCY(1) - y_temp(1)*PCY(1);...
     % 0 0 0 0 0 0 x_temp(1) y_temp(1) z_temp(1)   -x_temp(1)*PCZ(1) - y_temp(1)*PCZ(1)]
    
     % 14    14     1     0     0     0     0     0     0
    % 0     0     0    14    14     1     0     0     0
    % 0     0     0     0     0     0    14    14     1
     %%
    minMatrixComplete = zeros(size(PCX_temp,2)*3,11);
    for n = 1: numPatsSel%: numPoints
        % n=1
        minimizationMatrix = blkdiag(gridVecs(:,n)',gridVecs(:,n)',gridVecs(:,n)');
        minimizationMatrixRightPart = -gridVecs(1:2,n)'.*PC_temp(:,n);
        % put them together like ( minimizationMatrix
        % minimizationMatrixRightPart); -> 3x11 matrix
        % minimizationMatrix = 
        % minimizationMatrix()
        minMatrixComplete(3*n-2:3*n,:) = [minimizationMatrix minimizationMatrixRightPart];
    end
    %%
    
    %
    pTL = lsqminnorm(minMatrixComplete,PC_temp(:));
    scale = pTL(10)*x + pTL(11)*y +1;
    %
    %
    P = [pTL(1:3)'; pTL(4:6)';pTL(7:9)'; pTL(10), pTL(11), 1];



end