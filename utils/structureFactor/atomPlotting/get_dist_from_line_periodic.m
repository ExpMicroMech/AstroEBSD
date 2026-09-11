function [distance,extendReal]=get_dist_from_line_periodic(normal,origin,points)
    normal=normal./norm(normal);
    normalRealSpace = vector3d(normal.x,normal.y,normal.z);
    corners = get_lattice_corners(normal,'n_cells',1)
    bounds=alphaShape(corners(1,:)'*1.05,corners(2,:)'*1.05,corners(3,:)'*1.05,Inf);
    distanceQuery = -1.8*max(corners(:)):0.005:1.8*max(corners(:));
    % distanceQuery = repmat(distanceQuery,[8,1]);
    originNew = vector3d(origin.x,origin.y,origin.z);
    zeroOrigin = Miller({0,0,0},{0,1,0},{1,0,0},{0,0,1},{1,1,0},{0,1,1},{1,0,1},{1,1,1},origin.CS,'uvw');
    additions =[[0 0 0];[0 0 1];[0 1 0];[1 0 0];[1 1 0];[0 1 1];[1 0 1];[1 1 1]];
    %%
    atomsReturn = LineEquation(zeroOrigin',normal,distanceQuery);
    %%
    inside = bounds.inShape(atomsReturn.x,atomsReturn.y,atomsReturn.z);
    %%
    % sum(inside)
    % VecLens = zeros(8,1);
    % for i = 1:8
    %     VecLen = (max(distanceQuery(i,inside(i,:))-min(distanceQuery(i,inside(i,:)))));
    %     try
    %         VecLens(i) = VecLen;
    %     catch
    %         VecLens(i) = NaN;
    %     end
    % end
    %%
    nNew = normal/max(abs(normal.uvw))
    % nNew =repmat(nNew,[8,1])
     uRepLow = lsqminnorm(nNew.u',origin(1).u+0)
     vRepLow = lsqminnorm(nNew.v',origin(1).v+0)
     wRepLow = lsqminnorm(nNew.w',origin(1).w+0)

     uRepHigh = lsqminnorm(nNew.u',origin(1).u+1)
     vRepHigh = lsqminnorm(nNew.v',origin(1).v+1)
     wRepHigh = lsqminnorm(nNew.w',origin(1).w+1)
     
     maxRepsLow = [uRepLow, vRepLow, wRepLow];
     maxRepsLow(maxRepsLow>20) =0;
     [distLow, closestPlaneLow] = min(abs(maxRepsLow))

      maxRepsHigh = [uRepHigh, vRepHigh, wRepHigh];
     maxRepsHigh(maxRepsHigh>20) =0;
     [distHigh, closestPlaneHigh] = min(abs(maxRepsHigh))
     % vecLengthorig    
     % maxReps = ceil(maxReps)
%%
    atomsReturn = LineEquationUVW(origin(1)',normal,distanceQuery);
    inside = bounds.inShape(atomsReturn.x,atomsReturn.y,atomsReturn.z);
    VecLen = (max(distanceQuery(inside)-min(distanceQuery(inside))));
    originZero = atomsReturn(inside);
    originZero = originZero(1);
    originNew = originZero + VecLen*normal;
    originNew = Miller(mod(originNew.u,1), mod(originNew.v,1), mod(originNew.w,1),originNew.CS,'uvw');
    %
     
     % uvRep = lcm(maxReps(,vRep);
     % lcm(uvRep,wRep);
    %

       figure
    scatter3(originZero.x,originZero.y,originZero.z);
    % while abs(originZero-originNew)>1e-2
    for i = 1:5%max(abs(maxReps));       %%
        atomsReturn = LineEquationUVW(originNew',normal,distanceQuery);
        inside = bounds.inShape(atomsReturn.x,atomsReturn.y,atomsReturn.z);
        VecLen = (max(distanceQuery(inside)-min(distanceQuery(inside))))
        %%
        originZeroNew = atomsReturn(inside)
        originZeroNew = originZeroNew(1)
        %%
        originNew = atomsReturn(inside)
        originNew = originNew(end)
        originNew = Miller(mod(originNew.u,1), mod(originNew.v,1), mod(originNew.w,1),originNew.CS,'uvw')
        %%
        hold on
        plot3([originNew.x;  originZeroNew.x], [originNew.y; originZeroNew.y],[originNew.z; originZeroNew.z ]);
        scatter3(originZeroNew.x, originZeroNew.y, originZeroNew.z);
        scatter3(originNew.x, originNew.y, originNew.z)
        %%
        abs(originZero-originNew)
        % plot(bounds)
       
    end
   %%
   plot(bounds)
    %%
    %%
    VecLen = max(VecLens)
    extendReal = findMaxExtend(VecLen*normal);
    extend = ceil(extendReal);
    OriginN = origin+(0:1:extend-1)*VecLen*normal;
    OriginN = Miller(mod(OriginN.u,1), mod(OriginN.v,1), mod(OriginN.w,1),OriginN.CS,'uvw');
    OriginN = reshape(OriginN,[],extend);
% nor    %%
    originNew = vector3d(OriginN.x,OriginN.y,OriginN.z);

    % originNew = vector3d(origin.x,origin.y,origin.z)
    pointsNew = vector3d(points.x, points.y, points.z);
    % i=4
    % returnVec = LineEquation(originNew,normal,distanceQuery);
%%
    distance = zeros([size(originNew,1),size(originNew,2),size(points,1)]);
    % atomPos = 1
    for atomPos = 1:size(origin,1)
        distance(atomPos,:,:) =  (norm((originNew(atomPos,:) - pointsNew) - (dot((originNew(atomPos,:)-pointsNew),normal) * normal)))' ;
    end
    distance = reshape(distance,size(origin,1),[]);
    %%
    % atomsReturn = LineEquation(OriginN',normal,distanceQuery(1,:));
        % figure
    % imagesc(distance)
    % colorbar
    %%
    %%https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane
    %%
    
end
