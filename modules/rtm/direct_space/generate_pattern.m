function output_pattern=generate_pattern(screen_int,Settings_Cor,pcx,pcy,pcz,phi1,PHI,phi2)
    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz]
    % output_pattern=
    EBSP_Gnom( Settings_Cor,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen( EBSP_Gnom,gmatrix,screen_int);
end