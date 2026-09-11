Settings_Cor.size=[256,256];
eangs=extractions.(tile_name).oris(p_cur,:);
eangs=conv_G_to_EA(gmatrix_ref);

pcx_in=double(PC_best(1));
pcy_in=double(PC_best(2));
pcz_in=double(PC_best(3));
var_pcx=0.1;
var_pcy=0.1;
var_pcz=0.1;

phi1_in=double(eangs(1));
PHI_in=double(eangs(2));
phi2_in=double(eangs(3));
phi1
var_phi1=2*degree;
var_PHI=2*degree;
var_phi2=2*degree;
%%
phi2+var_phi2



%%
prob = optimproblem('ObjectiveSense','maximize');
% 
% %%
% 
% 
% 
pcx = optimvar('pcx',1,'Type','continuous','LowerBound',pcx_in-var_pcx,'UpperBound',pcx_in+var_pcx);
pcy = optimvar('pcy',1,'Type','continuous','LowerBound',pcy_in-var_pcy,'UpperBound',pcy_in+var_pcy);
pcz = optimvar('pcz',1,'Type','continuous','LowerBound',pcz_in-var_pcz,'UpperBound',pcz_in+var_pcz);

phi1 = optimvar('phi1',1,'Type','continuous','LowerBound',phi1_in-var_phi1,'UpperBound',phi1_in+var_phi1);
PHI = optimvar('PHI',1,'Type','continuous','LowerBound',PHI_in-var_PHI,'UpperBound',PHI_in+var_PHI);
phi2 = optimvar('phi2',1,'Type','continuous','LowerBound',phi2_in-var_phi2,'UpperBound',phi2_in+var_phi2);

%%
x0=[pcx_in,pcy_in,pcz_in,phi1,PHI,phi2]
f=@(x)dotproduct_optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor)
prob.Objective(f)
%%

options = optimoptions('fminunc','Algorithm','quasi-newton');
options.Display = 'iter';
%%
f(x0)
[x] = fminunc(f,x0,options)
%%

figure
subplot(1,2,1)
imagesc(compare_pat)
axis xy; axis equal;
subplot(1,2,2)
imagesc(Pattern_in)
axis xy; axis equal;
% normdotprod(normalize(EBSP_single.ebsp_cor),normalize(Pat_sim_eang_ref))
%%
dotproduct_optimization(x0,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor)
%%

function ndp=dotproduct_optimization(x,screen_int,Settings_Cor,PatternInfo,Pattern_in)
    pcx=x(1);
    pcy=x(2);
    pcz=x(3);
    phi1=x(4);
    PHI=x(5);
    phi2=x(6);
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    ndp=normdotprod(normalize(Pattern_in),normalize(compare_pat));
    
end

function output_pattern=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2)
    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz];
    % output_pattern=
    geometry=EBSP_Gnom(PatternInfo,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen(geometry,gmatrix,screen_int);
end

function ntd=normdotprod(A,B)
    ntd=sum(dot(A,B),'all')/(sqrt(sum(dot(A,A),'all'))*sqrt(sum(dot(A,A),'all')));
end