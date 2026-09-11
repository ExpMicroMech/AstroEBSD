function rot_mat=get_h_rot(h)
% get_h_rot -> find a rotation that puts h to the center of a projection
% 
%   rot_mat = get_h_rot(h) returns rotationmatrix to get h in the centre
%   h is a MTEX Miller variable
%   rot_mat is a MTEX rotation
%
%
%   nach h.j. bunge 1988
    h = normalize(h);
    h_new = Miller(h(1).h,h(1).k,h(1).l,'hkl',h.CS);
    h_new = normalize(h_new);
    % Phi= acos(h_new(1).l/sqrt(h_new(1).h^2+h(1).k^2+h_new(1).l^2));
    % %%;
    % phi2=acos(h_new(1).k/sqrt(h_new(1).h^2+h_new(1).k^2));

    Phi= acos(h_new(1).z/sqrt(h_new(1).x^2+h(1).y^2+h_new(1).z^2));
    %%;
    phi2=acos(h_new(1).y/sqrt(h_new(1).x^2+h_new(1).y^2));
    %%
    %%
%     rot_x=rotation.byAxisAngle(vector3d.X,angle(vector3d.Z,h(1)));
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% %     rot_x=rotation.byAxisAngle(vector3d.X,55*degree);
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% %     rot_mat=rot_x*rot_z;
% disp(h)
% disp(phi2)

% disp(Phi)
    if isnan(phi2)

        rot_mat=rotation.byEuler(0,Phi,0);
    else
        rot_mat=rotation.byEuler(0,Phi,phi2,"Bunge");
    end
    
end