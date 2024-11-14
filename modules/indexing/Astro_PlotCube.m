home; clear; 
%close all; - uncomment for regualr use

t1=clock; %start a timer

% Plugin locations
InputUser.MTEX_loc='C:\Communal_MatlabPlugins\mtex-5.2.beta2';
InputUser.Astro_loc='C:\Users\tbritton\Documents\GitHub\AstroEBSD_v2';


 run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m'));
 run(fullfile(InputUser.MTEX_loc,'startup_mtex.m'));

%%
InputUser.Phase_Folder='E:\Ben\phases';
InputUser.Phases={'Si'};

%create the crystal
[ ~,~,~,~,~, RTI_info ] = Phase_Builder_RTM(InputUser.Phases(1),InputUser.Phase_Folder);
cs=loadCIF(RTI_info.cif_file);

%these are the plane centres to plot on the pattern
CVectors.HKL=[0 1 1;0 0 1];

% CVectors.HKL=[0 0 4;
%               0 6 2;
%               5 0 2;
%               3 -6 2;
%               5  -6 0;
%               5 -4 2];


CVectors.colors={'r','g','b','y','m','k'};

%geometry set up
% pc=[0.5 0.5 0.5];
pc=[0.5 0.5 2.2];
%this is for a EBSD set up
% MicroscopeData.TotalTilt=-25*pi/180;

%this is for an ECP set up
MicroscopeData.TotalTilt=0;

%% Load the interpolants and cs list
pTime('Building the Phase Interpolants',t1);
%create the proper one
[screen_int_full,RTM_info_full,cslist_full] = Build_Interps(InputUser);

%% Build the crystal orientation

%orientation matricies
RTM_setup.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM_setup.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM_setup.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%Euler angles in
%in degrees
phi1=45;
PHI=0;
phi2=0;

%NOTE need to validate that this is the Eulerangle convention - I guessed here
G_best=RTM_setup.Rz(phi1*pi/180)*RTM_setup.Rx(PHI*pi/180)*RTM_setup.Rz(phi2*pi/180);

%% Generate a pattern

Settings_Cor.size=256; %simulated pattern in in pixels
PatternInfo.size=[Settings_Cor.size,Settings_Cor.size];

%build the detector geometry
[ Mean_EBSD_geom ] = EBSP_Gnom( PatternInfo,pc);

%set up the detector orientation
U.S = RTM_setup.Rx(MicroscopeData.TotalTilt);

%take the crystal orientation from the reference structure to the detector
U.O=G_best*inv(U.S);

%from the interpolation - extract the diffraction pattern
[ Pat_Refined ] = EBSP_gen( Mean_EBSD_geom,U.O*U.S,screen_int_full,0 ); %generate the EBSP for this iteration


%% Plot the band overlays

%the unit cell parameters are determined from the CIF file

UCell.a=norm(cs.aAxis);
UCell.b=norm(cs.bAxis);
UCell.c=norm(cs.cAxis);

UCell.alpha=cs.alpha*180/pi;
UCell.beta=cs.beta*180/pi;
UCell.gamma=cs.gamma*180/pi;

%equation 1
UCell.f=sqrt(1.0-( cosd(UCell.alpha)*cosd(UCell.alpha)...
    +cosd(UCell.beta)*cosd(UCell.beta)...
    +cosd(UCell.gamma)*cosd(UCell.gamma))...
    +2.0*cosd(UCell.alpha)*cosd(UCell.beta)*cosd(UCell.gamma));

%equation 2
UCell.ax = UCell.a * UCell.f/sind(UCell.alpha);

UCell.ay = UCell.a * (cosd(UCell.gamma)-cosd(UCell.alpha)*cosd(UCell.beta))...
    /sind(UCell.alpha);

UCell.az = UCell.a * cosd(UCell.beta);

%equation 3
UCell.by = UCell.b * sind(UCell.alpha);
UCell.bz = UCell.b * cosd(UCell.alpha);

%equation 4
UCell.cz = UCell.c;

%equation 5
UCell.StructureMat=[UCell.ax , 0,  0;
    UCell.ay , UCell.by, 0;
    UCell.az , UCell.bz, UCell.cz];

U.Astar=inv(UCell.StructureMat);

U.Kstar     =   U.Astar*U.O*U.S;

%rotate the HKLs into the detector frame
sym_matrix=cs.matrix;
% sym_matrix=symmetry_c;
% c=CVectors.HKL*U.Astar*U.O*U.S;
HKL.D=CVectors.HKL*U.Astar*U.O*U.S;
% CVectors.HKL=CVectors.inp;

%reduce the colours
HKL.col=CVectors.colors(1:size(CVectors.HKL,1));

if size(sym_matrix,3)>1
    for n=2:size(sym_matrix,3)
        HKL.D1=CVectors.HKL*U.Astar*sym_matrix(:,:,n)*U.O*U.S;
        HKL.D=[HKL.D;HKL.D1];
%         HKL.inp=[CVectors.inp;CVectors.HKL];
        HKL.col=[HKL.col CVectors.colors(1:size(CVectors.HKL,1));];
    end
end


% now plot pattern
%HKL - [X,Y,Z]
HKL.X=transpose(HKL.D(:,1));
HKL.Y=transpose(HKL.D(:,2));
HKL.Z=transpose(HKL.D(:,3));

%HKL - needed for hessian calculations
HKL.r=sqrt(HKL.X.^2+HKL.Y.^2+HKL.Z.^2);
HKL.kai=atan2(HKL.Y,HKL.X);
HKL.theta=acos(HKL.Z./HKL.r);

%Hessian construction
Hess.R_Hesse=10; %radius of the Hessian
Hess.d_Hesse=tan(0.5*pi-HKL.theta);
Hess.alpha_Hesse=acos(Hess.d_Hesse./Hess.R_Hesse);

Hess.alpha1_hkl=HKL.kai-pi+Hess.alpha_Hesse;
Hess.alpha2_hkl=HKL.kai-pi-Hess.alpha_Hesse;

%[C1x,C1y] to [C2x,C2y] are the coords on the screen
Hess.C1x=Hess.R_Hesse.*cos(Hess.alpha1_hkl);
Hess.C1y=Hess.R_Hesse.*sin(Hess.alpha1_hkl);
Hess.C2x=Hess.R_Hesse.*cos(Hess.alpha2_hkl);
Hess.C2y=Hess.R_Hesse.*sin(Hess.alpha2_hkl);



figure;
%plot the simulated pattern
s2=subplot(1,3,2); pPattern(Pat_Refined,Mean_EBSD_geom);
xlim([Mean_EBSD_geom.x_gn_min,Mean_EBSD_geom.x_gn_max]);
ylim([Mean_EBSD_geom.y_gn_min,Mean_EBSD_geom.y_gn_max]);
axis xy;
colormap('gray')

%plot the bands
num_HKL=size(Hess.C1x,2);
colour_set=[0.45 0.23 0.8];
colors_set_label=[0 0 0];

for n=1:num_HKL
    if HKL.Z(n)>0 %if upper hemisphere
        hold on
        l1=plot([Hess.C1x(n) Hess.C2x(n)],[Hess.C1y(n) Hess.C2y(n)],'-','LineWidth',2);
        l1.Color=HKL.col{n};
        colour_set=[colour_set;l1.Color];
        
    end
end




%% Plot this in a 3D render
%{
%figure
% s1=subplot(1,3,1);  pPattern(PatternCor,Mean_EBSD_geom); title('Input') %VERY USEFUL FUNCTION!!!
% s2=subplot(1,3,2); pPattern(Pat_Refined,Mean_EBSD_geom);
% 
% xlim([Mean_EBSD_geom.x_gn_min,Mean_EBSD_geom.x_gn_max]);
% ylim([Mean_EBSD_geom.y_gn_min,Mean_EBSD_geom.y_gn_max]);
% axis xy;
% colormap('gray')
% 
% %plot the bands
% num_HKL=size(Hess.C1x,2);
% colour_set=[0.45 0.23 0.8];
% colors_set_label=[0 0 0];
% 
% for n=1:num_HKL
%     if HKL.Z(n)>0 %if upper hemisphere
%         hold on
%         l1=plot([Hess.C1x(n) Hess.C2x(n)],[Hess.C1y(n) Hess.C2y(n)],'-','LineWidth',2);
%         l1.Color=HKL.col{n};
%         colour_set=[colour_set;l1.Color];
%         
%     end
% end
s3=subplot(1,3,3);

%plot the source point position
scatter3(0,0,0,100,'kx');

hold on;

%plot the pattern
z_pts_screen=Mean_EBSD_geom.xpts_screen*0+1;
s1=surf(Mean_EBSD_geom.xpts_screen,Mean_EBSD_geom.ypts_screen,z_pts_screen,Pat_Refined,'EdgeColor','none','FaceAlpha',0.9);
colormap('gray');

%plot the planes
for n=1:num_HKL
    hold on
    
    %work out if they fall on the screen
    H1=[Hess.C1x(n),Hess.C1y(n),1];
    H2=[Hess.C2x(n),Hess.C2y(n),1];
    H1to2=H2-H1;
    
    Lx_min=(Mean_EBSD_geom.x_gn_min-H1(1))./H1to2(1);
    Lx_max=(Mean_EBSD_geom.x_gn_max-H1(1))./H1to2(1);
    
    Ly_min=(Mean_EBSD_geom.y_gn_min-H1(2))./H1to2(2);
    Ly_max=(Mean_EBSD_geom.y_gn_max-H1(2))./H1to2(2);
    
    P1x=Lx_min*H1to2+H1;
    P2x=Lx_max*H1to2+H1;
    P3y=Ly_min*H1to2+H1;
    P4y=Ly_max*H1to2+H1;
    P_all=round([P1x;P2x;P3y;P4y],4);
        x_min=round(Mean_EBSD_geom.x_gn_min,4);
        x_max=round(Mean_EBSD_geom.x_gn_max,4);
        y_min=round(Mean_EBSD_geom.y_gn_min,4);
        y_max=round(Mean_EBSD_geom.y_gn_max,4);
        
        
    P_ok=(P_all(:,1) >= x_min & P_all(:,1) <= x_max & ...
        P_all(:,2) >= y_min & P_all(:,2) <= y_max);

    P_final=P_all(P_ok,:);
    
    if size(P_final,1)>=2
        p1=patch([0 P_final(1,1) P_final(2,1) 0],[0 P_final(1,2) P_final(2,2) 0],[0 1 1 0],[1 1 1 1]);
        scatter3(P_final(1,1),P_final(1,2),P_final(1,3),10,'k','filled');
        scatter3(P_final(2,1),P_final(2,2),P_final(1,3),10,'k','filled');
        p1.FaceColor=HKL.col{n};
        p1.EdgeColor=0.7*p1.FaceColor;
        p1.FaceAlpha=0.7;        
    end
    

end

unit_cube_plot_rotmat(cs,U.O, U.S,0.02)

axis equal;


%plot the sample as a disc

r_sample=0.6;
theta=0:360;
sam_x=r_sample*sind(theta);
sam_y=r_sample*cosd(theta);
sam_z=sam_x*0+0;

sam_det=transpose([sam_x;sam_y;sam_z])*U.S;

p_sample=patch(sam_det(:,1),sam_det(:,2),sam_det(:,3),[1 1 1]);
p_sample.FaceAlpha=0.9;
scatter3(sam_det(1,1),sam_det(1,2),sam_det(1,3),'rx');
scatter3(sam_det(91,1),sam_det(91,2),sam_det(91,3),'gx');

xlim([-1.5 1.5])
ylim([-1.5 1.5]);
zlim([-1 2]);
axis on;
xlabel('X');
ylabel('Y');
zlabel('Z');


%% To plot a pole figure in the sample frame

CVectors.HKL_pf=[0 0 1];
HKL_pf.D=CVectors.HKL_pf*U.Astar*U.O;
if size(sym_matrix,3)>1
    for n=2:size(sym_matrix,3)
        HKL_pf.D1=CVectors.HKL_pf*U.Astar*sym_matrix(:,:,n)*U.O;
        HKL_pf.D=[HKL_pf.D;HKL_pf.D1];
    end
end

HKL_pf.D_sample=HKL_pf.D;


figure;
s1=subplot(1,2,1);

theta=0:360;
cir_x=sind(theta);
cir_y=cosd(theta);
plot(cir_x,cir_y);
hold on
plot([-1 1],[0 0]);
plot([0 0],[-1 1]);

% for n=1:num_HKL
for n=1:4
    Q=HKL_pf.D_sample(n,:)./norm(HKL_pf.D_sample(n,:)); %force to unit length - i.e. hits the sphere
    
    if Q(3) > 0 %find the pole on the other hemisphere
        S=[0 0 -1];
    else
        S=[0 0 1];
    end
    
    QS=S-Q; %find the vector QS
    mu=-S(3)/QS(3);
    M=S+mu*QS;
    
    scatter(M(1),M(2),'x');
    
end
axis equal
scatter(1,0,'r');
scatter(0,1,'g');
s1.XDir='reverse';
s1.YDir='reverse';

s2=subplot(1,2,2);
unit_cube_plot_rotmat(cs,U.O, eye(3),0.05)
axis on;
xlabel('X');
ylabel('Y');
zlabel('Z');
s2.XDir='reverse';
s2.YDir='reverse';

%% To plot a pole figure in the sample frame

CVectors.HKL_pf=[1 0 0];
HKL_pf.D=CVectors.HKL_pf*U.Astar*U.O;
if size(sym_matrix,3)>1
    for n=2:size(sym_matrix,3)
        HKL_pf.D1=CVectors.HKL_pf*U.Astar*sym_matrix(:,:,n)*U.O;
        HKL_pf.D=[HKL_pf.D;HKL_pf.D1];
    end
end

HKL_pf.D_sample=HKL_pf.D;


figure;
s1=subplot(1,2,1);

theta=0:360;
cir_x=sind(theta);
cir_y=cosd(theta);
plot(cir_x,cir_y);
hold on
plot([-1 1],[0 0]);
plot([0 0],[-1 1]);

% for n=1:num_HKL
for n=1:4
    Q=HKL_pf.D_sample(n,:)./norm(HKL_pf.D_sample(n,:)); %force to unit length - i.e. hits the sphere
    
    if Q(3) > 0 %find the pole on the other hemisphere
        S=[0 0 -1];
    else
        S=[0 0 1];
    end
    
    QS=S-Q; %find the vector QS
    mu=-S(3)/QS(3);
    M=S+mu*QS;
    
    scatter(M(1),M(2),'x');
    
end
axis equal
scatter(1,0,'r');
scatter(0,1,'g');
s1.XDir='reverse';
s1.YDir='reverse';

s2=subplot(1,2,2);
unit_cube_plot_rotmat(cs,U.O, eye(3),0.05)
axis on;
xlabel('X');
ylabel('Y');
zlabel('Z');
s2.XDir='reverse';
s2.YDir='reverse';

%}
