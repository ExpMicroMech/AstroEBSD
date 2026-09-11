function [m_unique,f_hkl_sorted] = calcStructureFactor(PhaseInput,lambda,InputUser)
% clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).
% run([InputUser.Astro_loc,'\start_AstroEBSD']);
% InputUser.Phase_Folder = Phase_Folder
% % Mg 	 
% TODO: at some point think about putting inputUser in a configuration file
%
%
cutoff = 0.2;
if isfield(InputUser, 'minFHKL')
    cutoff = InputUser.minFHKL;
end
alignEdax = 0;
if isfield(InputUser,'AlignEdax')
    alignEdax = InputUser.AlignEdax;
end
alignEMsoft = 0;
if isfield(InputUser,'AlignEMsoft')
    alignEMsoft = InputUser.AlignEMsoft;
end
dCutoff = 0.5;
if isfield(InputUser, 'minD')
    dCutoff = InputUser.minD;
end

maxInd = 10;
if isfield(InputUser, 'maxInd')
    maxInd = InputUser.maxInd;
end
%%
% InputUser.Phase_Input={'Nb6Ni7_em'};
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( PhaseInput,InputUser.Phase_Folder );
% [cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 
%%
[cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
if alignEdax == 1
    cs = crystalAlignEdax(cs);
end
if alignEMsoft == 1
    cs = crystalAlignEMSoft(cs);
end
% a1 = [5.4204, 2.1735, 1.2269, 2.3073];
% b1 = [2.8275, 79.2611, 0.3808, 7.1937];
% c = 0.8584;
%%
% )
%%
output=get_atom_positions(atoms,cs);
%%
% atoms_sym=unique(output.all_atoms);
% atoms_sym=atoms_sym(all(atoms_sym.hkl<0.998,2))
% % atoms_sym=atoms_sym(all(atoms_sym.hkl>-0.001,2))
% uvw=round(atoms_sym.uvw,8);
% uvw_unique=unique(uvw,'rows');
% uvw_unique=uvw_unique(all(uvw_unique>=0,2),:)
%%
atoms_out = output.all_atoms;
% [C,ia,ic] = unique(output.all_atoms.xyz,'rows')
[~,ia,~] = unique(round(atoms_out.hkl,4),'rows');
atoms_sym=atoms_out(ia);
insideUnitCell = ~any(atoms_sym.uvw>=0.998,2);
hkl_unique = atoms_sym(insideUnitCell);
% atoms_sym(ia)
% scatter3()
[a,b,c] = getAtomicFormFactorParameters   (output.elements(ia(insideUnitCell)));
%%
r_x=cat(1,hkl_unique.uvw(:,1));
r_y=cat(1,hkl_unique.uvw(:,2));
r_z=cat(1,hkl_unique.uvw(:,3));
%%
theta=[0:0.1:180]*degree;
f=zeros(size(r_x,1),1);
%%
% figure
% plot(theta/degree,f.^2)
% lambda=InputData.lambda_p;
% lambda = 1.5

% for i=1:size(r_x,1)
%     f(i, 1)= f_xray(a1, b1, c, theta(count),lambda);    
% end
% f_hkl=[f_hkl ,calc_f_hkl(f,h,k,l,r_x,r_y,r_z)];
% count=count+1;
%
d_track=[];
hkl=[];
m=[];
% maxInd = 10;
%% from Mtex
lattice_a = abs(cs.aAxis);
lattice_b = abs(cs.bAxis);
lattice_c = abs(cs.cAxis);
maxLattice = max([lattice_a,lattice_b,lattice_c]);
%%
hMax = floor(maxInd*lattice_a/maxLattice)
kMax = ceil(maxInd*lattice_b/maxLattice)
lMax = ceil(maxInd*lattice_c/maxLattice)
[h_grid,k_grid,l_grid]=meshgrid(-hMax:1:hMax,...
                        -kMax:1:kMax ,...
                        -lMax:1:lMax  );
filter = find([h_grid(:)+k_grid(:)+l_grid(:)]==0);
hkls = [h_grid(~filter), k_grid(~filter),l_grid(~filter)];
% hkls = hkls([h_grid(:)+k_grid(:)+l_grid(:)]~=0);
m = Miller(h_grid(:),k_grid(:),l_grid(:),cs);
m = m(2:end); %remove zero
%%
m_unique = unique(m);
m_unique = m_unique(m_unique.dspacing~=inf)
% m_unique.multiplicity
[~,dSpacingIndex] = sort(m_unique.dspacing,'descend');
m_unique = m_unique(dSpacingIndex);
%%
% dCutoff = 0.2;
m_unique  = m_unique(m_unique.dspacing>dCutoff); %cutoff small bands
m_unique = m_unique(lambda/(2*m_unique.dspacing)<1); %cutoff >180°
%%

hkl=[];
theta=[];
count = 1;
% for i=1:size(m_unique,1)
% d = m_unique(i).dspacing;
% disp(d);
    % if lambda/(2*d) <1;
        d = m_unique.dspacing;
        % n=check_order(m_unique(i));
        theta=asin(lambda./(2*d));
        % d_track=[d_track d];
        % for k=1:size(r_x,1)
        
        f = f_electron(a, b, c, theta,lambda);    
        % end
        h=m_unique.h;
        k=m_unique.k;
        l=m_unique.l;
        f_hkl=calc_f_hkl(f',h,k,l,r_x,r_y,r_z);
        % hkl=[sprintf("%i %i %i",int8(h),int8(k),int8(l))];
        if strcmp(m_unique.CS.lattice,'trigonal') || strcmp(m_unique.CS.lattice,'hexagonal')
            hkl = int8(m_unique.hkil);
        else
            hkl = int8(m_unique.hkl);
        end
        % count=count+1;
    % end
%     disp(hkl)
% end
%%
%%
% calc_f_hkl(f,0,0,12,r_x,r_y,r_z)

% n=check_order(m_unique(i));
%%
I=(abs(f_hkl).^2);
I_norm = I./max(I(:));
I_annotate = I_norm(I_norm>cutoff);

% figure
% cutoff = I./max(I(:))>0.
% scatter(2*theta/degree,I./max(I(:)))
% i_annotate =
hkl_annotate = hkl(I_norm>cutoff,:);
hkl_plot = [];
if size(hkl_annotate,2)==4
    for i = 1:size(hkl_annotate,1)
        hkl_plot = [hkl_plot, sprintf("(%i %i %i %i)",hkl_annotate(i,1),hkl_annotate(i,2),hkl_annotate(i,3),hkl_annotate(i,4))];
    end
else
    for i = 1:size(hkl_annotate,1)
        hkl_plot = [hkl_plot, sprintf("(%i %i %i)",hkl_annotate(i,1),hkl_annotate(i,2),hkl_annotate(i,3))];
    end
end
% text(2*theta(I_norm>cutoff)/degree,I_norm(I_norm>cutoff),hkl_plot);
[f_hkl_sorted,intensityIndex] = sort(abs(f_hkl),'descend');

m_unique = m_unique(intensityIndex);

%%
function s=s_crystal(R_r,R_k,f_j,f_k,lambda,theta)
    q = 4 * pi / lambda * sin( theta / 2);
end

%
function n=check_order(h)
    hkl=uint8(h.hkl);
    select=[1 2 3];
    %%
    [~,ind]=max(hkl);
    select(ind)=[];
    A=[hkl(ind),hkl(ind)];
    B=[hkl(select(2)),hkl(select(1))];
    n=min(gcd(A,B));

end

end

