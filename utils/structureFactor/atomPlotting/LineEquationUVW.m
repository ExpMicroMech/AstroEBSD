function returnVec = LineEquationUVW(origin,direction,distance)
    
    direction= normalize(direction);
    % direction = vector3d(direction.x,direction.y,direction.z);
    % origin = vector3d(origin.x,origin.y,origin.z);
    returnVec = origin + distance*direction;    
end