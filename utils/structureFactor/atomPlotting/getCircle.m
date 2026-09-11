function [X,Y] = getCircle(numPoints,radius) 
    theta = linspace(0,2*pi,numPoints);
    [THETA,RADIUS] = ndgrid(theta,radius);
    % [THETAREal,RADIUSReal] = ndgrid(theta,rReal);
    [X,Y] = pol2cart(THETA,RADIUS);
end