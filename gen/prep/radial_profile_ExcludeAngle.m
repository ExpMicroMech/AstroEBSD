function [Tics, Average] = radial_profile_ExcludeAngle(data, radial_step, angularRangeToExclude)
    x = (1:size(data,2)) - size(data,2)/2 - 1;
    y = (1:size(data,1)) - size(data,1)/2 - 1;
    [X, Y] = meshgrid(x, y);
    
    angleMatrix = abs(atan(Y ./ X));
    
    Z_integer = round(abs(X + 1i*Y) / radial_step) + 1;

    radii = abs(X + 1i * Y);
    Tics    = accumarray(Z_integer(:), radii(:), [], @mean);
    % Average = accumarray(Z_integer(:), data(:), [], @nanmean);

    indicesToExclude1a = find(angleMatrix <= angularRangeToExclude);
    indicesToExclude2a = find((pi/2-angularRangeToExclude < angleMatrix) & (angleMatrix <= pi/2));

    indicesToExclude = [indicesToExclude1a; indicesToExclude2a];

    data2 = data; data2(indicesToExclude) = NaN;
    Average = accumarray(Z_integer(:), data2(:), [], @nanmean);
end