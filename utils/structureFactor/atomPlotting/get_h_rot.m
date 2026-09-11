function rot_mat=get_h_rot(h)
% nach h.j. bunge 1988
    m=h.xyz
    h=m(1)
    k=m(2)
    l=m(3)
    Phi= acos(l/sqrt(h^2+k^2+l^2));
    %%;
    phi2=acos(k/sqrt(h^2+k^2));
    %%
%     rot_x=rotation.byAxisAngle(vector3d.X,angle(vector3d.Z,h(1)));
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
% %     rot_x=rotation.byAxisAngle(vector3d.X,55*degree);
%     rot_z=rotation.byAxisAngle(vector3d.Z,angle(vector3d.X,h(1)));
%     rot_mat=rot_x*rot_z;
    if isnan(phi2)
        rot_mat=rotation.byEuler(0,Phi,0);
    else
        rot_mat=rotation.byEuler(0,Phi,phi2);
    end
    
 end