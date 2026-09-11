function [x_u,y_u,z_u]=restrict_to_unit_cell(x,y,z)
    x_u=x-floor(x)
    y_u=y-floor(y)
    z_u=z-floor(z)
end