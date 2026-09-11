function [pcGrid, rotGrid] =  setupSimulationGrid(gridVecs,gridVar1,gridVar2,pc,rot_orig)

    gridMatrix = zeros(size(gridVecs,1),size(gridVecs,1),6);
    % gridMatrix(:,:,1:3) = gridMatrix(:,:,1:3)+1
    [Grid1,Grid2] =  meshgrid(gridVecs(:,gridVar1), gridVecs(:,gridVar2));
    %NOTE: we go with generic variable names, as the code is supposed to be generic
    
    
    gridMatrix(:,:,gridVar1) = Grid1;
    gridMatrix(:,:,gridVar2) = Grid2;
    
    pcGrid = cat(3,pc(1)+gridMatrix(:,:,1),pc(2)+gridMatrix(:,:,2),pc(3)+gridMatrix(:,:,3));
    
    % rotNew = rotation.byAxisAngle(vector3d.Y,RotYGrid);
    % rotGrid =  cat(3,pc(1)+gridMatrix(:,:,1),pc(2)+gridMatrix(:,:,2),pc(3)+gridMatrix(:,:,3)) 
    % see how easy we can set this up
    % for rot grid, we need to combine rotations.
    rotx = rotation.byEuler(0,gridMatrix(:,:,4),0,'Bunge');
    roty = rotation.byAxisAngle(vector3d.Y,gridMatrix(:,:,5),'Bunge');
    rotz = rotation.byEuler(0,0,gridMatrix(:,:,6),'Bunge');
    rotNew = rotx.*roty.*rotz;
    rotGrid = rot_orig* rotNew;
    rotGrid = reshape(rotGrid, size(rotNew,1),size(rotNew,2));
end