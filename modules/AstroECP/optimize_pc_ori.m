function [pc_out,eangs,peakhight] =optimize_pc_ori(pc_in,eangs,PatternIn,Settings_Cor,screen_int,varargin)
% OPTIMIZE_PC_ORI
%
%Optimizes a orientation and patterncenter based on the normalised
%crosscorrelation normxcorr2
%[pc_out,eangs,peakhight]=
%OPTIMIZE_PC_ORI(pc_in,eangs,PatternIn,Settings_Cor,screen_int,varargin)
% - flag 'Mode', with either 'fmincon' or 'fminsearch' to select for
% optimisation method
%- Settings_Cor: define either max_var_pc or max_var_pc_x,max_var_pc_y and
%max_var_pc_z as field names to set the constraints for the refienement.
%- pc_in: 3x1 double, containing pcx, pcy, DD
%- eangs: 3x1 double, contianing Euleranglex in ZXZ convention in radiants
%- screen_int, contains interpolants of the masterpattern
%- Settings_Cor, struct, with cropping and correction infos, as well as
%refinement constraints.
%
%
%Outputs
%pc_out: refined pattern_center
%eangs: refined orientation
%peakhigh: normalised cross product between simulation and experiment after
%refinement
%
%

% Settings_Cor.size=256;
% Settings_Cor.SquareCrop=0;
% eangs=extractions.(tile_name).oris(p_cur,:);
% eangs=conv_G_to_EA(gmatrix_pseudo);
mode='fmincon'
for i=1:size(varargin)
    if(strcmp(varargin{i},'Mode'))
        mode=varargin{i+1}
    end

end
pcx_in=double(pc_in(1));
pcy_in=double(pc_in(2));
pcz_in=double(pc_in(3));
fields={'max_var_pc_x','max_var_pc_y','max_var_pc_y'}
if any(isfield(Settings_Cor,fields))
    try
    var_pcx=Settings_Cor.max_var_pc_x;
    var_pcy=Settings_Cor.max_var_pc_y;
    var_pcz=Settings_Cor.max_var_pc_z;
    catch
        error('you have spacefied a constrained for one pc variblae but not all of them')
    end
else
    var_pcx=Settings_Cor.max_var_pc;
    var_pcy=Settings_Cor.max_var_pc;
    var_pcz=Settings_Cor.max_var_pc;
end

phi1_in=double(eangs(1));
PHI_in=double(eangs(2));
phi2_in=double(eangs(3));
var_phi1=Settings_Cor.max_var_ori;
var_PHI=Settings_Cor.max_var_ori;
var_phi2=Settings_Cor.max_var_ori;

[EBSP_single.ebsp_cor,EBSP_single.PatternInfo]=EBSP_BGCor(PatternIn,Settings_Cor );
EBSP_single.PatternInfo.size=size(EBSP_single.ebsp_cor);
% disp(Settings_Cor)
%%
x0=[pcx_in,pcy_in,pcz_in,phi1_in,PHI_in,phi2_in];
f=@(x)(dotproduct_optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor));
options = optimoptions('fmincon');%,'Algorithm','quasi-newton');
options.Display = 'iter';
%%
l=[pcx_in-var_pcx;pcy_in-var_pcy;pcz_in-var_pcz;phi1_in-var_phi1;PHI_in-var_PHI;phi2_in-var_phi2];
u=[pcx_in+var_pcx;pcy_in+var_pcy;pcz_in+var_pcz;phi1_in+var_phi1;PHI_in+var_PHI;phi2_in+var_phi2];
    % < pcx_in < pcx_in+var_pcx;
 % < x(2) < pcy_in+var_pcy;
% pcz_in-var_pcz < x(3) < pcz_in+var_pcz;
% f(x0)
%%
if mode=='fmincon'
    [x,fval] = fmincon(f,x0,[],[],[],[],l,u,[],options);
elseif mode=='fminsearch'
    [x,fval] = fminsearch(f,x0)%,l,u)%,[],options)
else
    error('you have specified a none supported mode')
end
%%
pc_out=[x(1),x(2),x(3)];
eangs=[x(4),x(5),x(6)];
% f(x0)
peakhight=1-fval;

% %%
% extractions.(tile_name).rotation_data_refined(:,:,p_cur)=RTM.Rz(x(6))*RTM.Rx(x(5))*RTM.Rz(x(4));
% extractions.(tile_name).refined_pc(:,p_cur)=[x(1),x(2),x(3)];
% %%
% pPattern
% %%
% geometry=EBSP_Gnom(EBSP_single.PatternInfo,[x(1),x(2),x(3)]);
%%
[cmap_div]=cbrewer('div','RdBu',30);
cmap_div(cmap_div<0)=0;
% 
% geometry=EBSP_Gnom(EBSP_single.PatternInfo,[pcx_in,pcy_in,pcz_in]);
% % % %
% compare_pat=generate_pattern(screen_int,Settings_Cor,EBSP_single.PatternInfo,x(1),x(2),x(3),x(4),x(5),x(6));
% % compare_pat=generate_pattern(screen_int,Settings_Cor,EBSP_single.PatternInfo,pcx_in,pcy_in,pcz_in,phi1_in,PHI_in,phi2_in);
% % compare_pat=generate_pattern(screen_int,Settings_Cor,EBSP_single.PatternInfo,x(1),x(2),x(3),x(4),x(5),x(6));
% % % compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
% figure
% subplot(1,3,1)
% imagesc(geometry.x_screen,geometry.y_screen,compare_pat)
%  axis image;axis xy;
% colormap gray;
% subplot(1,3,2)
% imagesc(geometry.x_screen,geometry.y_screen,EBSP_single.ebsp_cor)
% axis image;axis xy; 
% colormap gray;
% ax_par=subplot(1,3,3)
% imagesc(geometry.x_screen,geometry.y_screen,normalize(EBSP_single.ebsp_cor)-normalize(compare_pat))
% axis image;axis xy;
% colormap gray;
%   clim([-3,3])
%         colormap(ax_par,cmap_div)
%         colorbar
% normdotprod(normalize(EBSP_single.ebsp_cor),normalize(Pat_sim_eang_ref))
%%
% dotproduct_optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor)
%%
% C=normxcorr2(compare_pat,EBSP_single.ebsp_cor);
%%
% figure
% imagesc(C)
% [ypeak,xpeak] = find(C==max(C(:)));
%%

function ndp=dotproduct_optimization(x,screen_int,Settings_Cor,PatternInfo,Pattern_in)
    % DOTPRODUCT_OPTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %
    pcx=x(1);
    pcy=x(2);
    pcz=x(3);
    phi1=x(4);
    PHI=x(5);
    phi2=x(6);
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    %%
    
    normcross=normxcorr2(Pattern_in,compare_pat);
    ndp=1-max(normcross(:));
    %%
    % imagesc(normcross)
    % colorbar()
    % %%
    % figure
    % imagesc(compare_pat)
    % figure
    % imagesc(Pattern_in)
    % ndp=1-normdotprod(normalize(Pattern_in),normalize(compare_pat));
end

function output_pattern=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2)
    % GENERATE_PATTERN
    % calculates the gnomonic projection based on dynamical simulation and
    % pattern geometry information
    %
    %

    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz];
    % output_pattern=
    %%
    geometry=EBSP_Gnom(PatternInfo,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen(geometry,gmatrix,screen_int);
%     output_pattern= EBSP_BGCor(output_pattern)
end

function ntd=normdotprod(A,B)
    ntd=sum(dot(A,B),'all')/(sqrt(sum(dot(A,A),'all'))*sqrt(sum(dot(A,A),'all')));
end
end