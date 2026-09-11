function optimize_pc_ori(pc,eangs,EBSP_single,Settings_Cor)
Settings_Cor.size=256;
Settings_Cor.SquareCrop=0;
eangs=extractions.(tile_name).oris(p_cur,:);
eangs=conv_G_to_EA(gmatrix_pseudo);

pcx_in=double(pc_in(1));
pcy_in=double(pc_in(2));
pcz_in=double(pc_in(3));
var_pcx=max_var_pc;
var_pcy=max_var_pc;
var_pcz=max_var_pc;

phi1=double(eangs(1));
PHI=double(eangs(2));
phi2=double(eangs(3));
var_phi1=2*degree;
var_PHI=2*degree;
var_phi2=2*degree;

[EBSP_single.ebsp_cor,EBSP_single.PatternInfo]=EBSP_BGCor(PatternIn,Settings_Cor );

%%
% prob = optimproblem('ObjectiveSense','maximize');
% 
% %%
% 
% 
% 
% pcx = optimvar('pcx',1,'Type','double','LowerBound',pcx_in-var_pcx,'UpperBound',pcx_in+var_pcx);
% pcy = optimvar('pcy',1,'Type','double','LowerBound',pcy_in-var_pcy,'UpperBound',pcx_in+var_pcy);
% pcz = optimvar('pcz',1,'Type','double','LowerBound',pcz_in-var_pcz,'UpperBound',pcx_in+var_pcz);
% 
% phi1 = optimvar('phi1',1,'Type','double','LowerBound',phi1_in-var_phi1,'UpperBound',pcx_in+var_phi1);
% PHI = optimvar('PHI',1,'Type','double','LowerBound',PHI_in-var_PHI,'UpperBound',pcx_in+var_PHI);
% phi2 = optimvar('phi2',1,'Type','double','LowerBound',phi2_in-var_phi2,'UpperBound',pcx_in+var_phi2);
%%
x0=[pcx_in,pcy_in,pcz_in,phi1,PHI,phi2]
f=@(x)(dotproduct_optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor))
options = optimoptions('fmincon')%,'Algorithm','quasi-newton');
options.Display = 'iter';
%%
l=[pcx_in-var_pcx;pcy_in-var_pcy;pcz_in-var_pcz;phi1-var_phi1;PHI-var_PHI;phi2-var_phi2]
u=[pcx_in+var_pcx;pcy_in+var_pcy;pcz_in+var_pcz;phi1+var_phi1;PHI+var_PHI;phi2+var_phi2]
    % < pcx_in < pcx_in+var_pcx;
 % < x(2) < pcy_in+var_pcy;
% pcz_in-var_pcz < x(3) < pcz_in+var_pcz;
f(x0)
%%
% [x] = fmincon(f,x0,[],[],[],[],l,u,[],options)
[x] = fminsearch(f,x0)%,l,u)%,[],options)
%%
f(x0)

%%
extractions.(tile_name).rotation_data_refined(:,:,p_cur)=RTM.Rz(x(6))*RTM.Rx(x(5))*RTM.Rz(x(4));
extractions.(tile_name).refined_pc(:,p_cur)=[x(1),x(2),x(3)];
%%
% pPattern
%%
geometry=EBSP_Gnom(EBSP_single.PatternInfo,[x(1),x(2),x(3)]);
%
% compare_pat=generate_pattern(screen_int,Settings_Cor,EBSP_single.PatternInfo,x(1),x(2),x(3),x(4),x(5),x(6));
compare_pat=generate_pattern(screen_int,Settings_Cor,EBSP_single.PatternInfo,x(1),x(2),x(3),x(4),x(5),x(6));
compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
figure
subplot(1,3,1)
imagesc(geometry.x_screen,geometry.y_screen,compare_pat)
 axis image;axis xy;
colormap gray;
subplot(1,3,2)
imagesc(geometry.x_screen,geometry.y_screen,EBSP_single.ebsp_cor)
axis image;axis xy; 
colormap gray;
ax_par=subplot(1,3,3)
imagesc(geometry.x_screen,geometry.y_screen,normalize(EBSP_single.ebsp_cor)-normalize(compare_pat))
axis image;axis xy;
colormap gray;
  clim([-3,3])
        colormap(ax_par,cmap_div)
        % colorbar
% normdotprod(normalize(EBSP_single.ebsp_cor),normalize(Pat_sim_eang_ref))
%%
dotproduct_optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor)
%%
C=normxcorr2(compare_pat,EBSP_single.ebsp_cor);
%%
figure
imagesc(C)
[ypeak,xpeak] = find(C==max(C(:)));
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
    normcross=normxcorr2(normalize(Pattern_in),normalize(compare_pat));
    ndp=1-max(normcross(:));
    
end

function output_pattern=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2)
    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz];
    % output_pattern=
    %%
    geometry=EBSP_Gnom(PatternInfo,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen(geometry,gmatrix,screen_int);
end

function ntd=normdotprod(A,B)
    ntd=sum(dot(A,B),'all')/(sqrt(sum(dot(A,A),'all'))*sqrt(sum(dot(A,A),'all')));
end