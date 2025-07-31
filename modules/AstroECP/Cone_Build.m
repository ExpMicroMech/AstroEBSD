function [bands,cone3d] = Cone_Build(n_hkl,Astar,lambda_p,rotmat,params)
%CONE_BUILD Calculate the Kossel cones 
% Used to plot on the surface of a sphere,
% and also on gnomonic projection plane
%
% TBB Dec 2023

% params=360; %how many points on the conic section to plot

%convert to cartesian, and measure the lattice spacing + bragg angle

n_cart=n_hkl*Astar*rotmat;
n_hat=n_cart/norm(n_cart);
d_hkl=1/norm(n_cart);

theta_hkl=asin(norm(n_cart)/2*lambda_p);

%calculate the cone opening ratio
%Note that I can't remember why this is 2/tan(theta) rather than 1/tan(theta)
%It is indeed 1/tan(theta) - theta formula on line 16 is fixed - Tianbi Feb 07 2024.
cone_a=1/tan(theta_hkl); 

%do the axis angle rotation


if abs(abs(n_hat(3))-1 ) > 1E-6
    v_ang=acos(n_hat(3));
    z_cone=[0 0 1];
    v_rot=cross(z_cone,n_hat);
    v_rot=v_rot/norm(v_rot);
    M_nhat = makehgtform('axisrotate',v_rot,v_ang);
else
    M_nhat=eye(4);
end


% the up and down cone transformations
M_up=makehgtform('translate',[0,0,0]);
M_down=makehgtform('translate',[0,0,0],'yrotate',pi);

%the combined transformations
T_up=M_nhat*M_up;
T_down=M_nhat*M_down;


%% do the central band

%the conic sections
param=linspace(0,2*pi,params);

%Conic ring of the plane
r_sphere=sqrt(1+cone_a^2);
RC_u=[cone_a*cos(param);cone_a*sin(param);(0*param+1);param*0+r_sphere]./r_sphere;

%the center line
RN=[cos(param);sin(param);0*param;0*param+1];

%apply the transform
bands.centre=T_up*RN;
bands.upper=T_up*RC_u;
bands.lower=T_down*RC_u;

%% Limit to +Z plotting
%limit to that stuff which falls +ve - good for gnomonic plotting
bands.centre_zp=bands.centre(:,bands.centre(3,:)>0);
bands.upper_zp=bands.upper(:,bands.upper(3,:)>0);
bands.lower_zp=bands.lower(:,bands.lower(3,:)>0);

bands.theta=theta_hkl;
bands.n_hat=n_hat;
bands.d_hkl=d_hkl;

%% Create the cones
if nargout == 2
    %create a cylinder / cone
    [X2,Y2,Z2]=cylinder([0 cone_a],50 );

    %normalize to the surface of the sphere
    X2=X2/r_sphere;
    Y2=Y2/r_sphere;
    Z2=Z2/r_sphere;

    %transform into the coordinate sspaces

    %use 4 index transformations
    R2_u=[X2(:)';Y2(:)';Z2(:)';X2(:)'*0+1];

    %transform
    R2_ut=T_up*R2_u;
    R2_dt=T_down*R2_u;

    %house keeping
    cone3d.X_up=reshape(R2_ut(1,:),size(X2));
    cone3d.Y_up=reshape(R2_ut(2,:),size(X2));
    cone3d.Z_up=reshape(R2_ut(3,:),size(X2));

    cone3d.X_down=reshape(R2_dt(1,:),size(X2));
    cone3d.Y_down=reshape(R2_dt(2,:),size(X2));
    cone3d.Z_down=reshape(R2_dt(3,:),size(X2));
end

end

