%% In this script we want to analyse the interlayer spacings of various mu phases. 
% Therefore we load a cif file replicate all atoms acording to symmetry
% operations and calculate the minimum z spacing between different layers.
%%
clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|(ORII5_loc)|MTEX_loc).

phaseList = ["Nb47Ni53_joubert","Nb51Ni49_joubert","Nb53Ni47_joubert","Nb55Ni45_joubert","Nb60Ni40_joubert"]
phaseList = ["Nb7Co6_em_relaxed_ray", "Nb6Co7_em_relaxed_ray"]
phaseList = ["Nb462Ni534_matProj_connetable","Nb496Ni504_joubert","Nb500Ni500_Eremenko_ILSpacings","Nb518Ni482_joubert","Nb530Ni470_joubert","Nb7Co6_em_relaxed_ray", "Nb6Co7_em_relaxed_ray","Nb534Ni462_Hunt","Nb534Ni462_lomnytska","Nb534Ni462_matProj_kripyakevich","Nb553Ni447_joubert","Nb569Ni431_joubert"]
phaseList = ["NbNi_dt0.400_c26.7_a4.9","NbNi_dt0.250_c26.7_a4.9"]
phaseList = ["Nb6Ni7_relaxed_ray","Nb7Ni6_relaxed_ray"]
phaseList = ["Nb6Ni7_relaxed_nisa"]
% p
% subFolder = 'NbCoLowNb_LatticeMatched'
subFolder = "ILTestingsForPublication"
PhasesToRead = dir(fullfile(Phase_Folder,sprintf('phasefiles/%s/*.pha',subFolder)))'
InputUser.Phase_Input = cell(size(PhasesToRead)-1);
phaseList  =[]
for i=1:numel(PhasesToRead) % -1 for debugging
    InputUser.Phase_Input{i} = PhasesToRead(i).name(1:end-4);
    phaseList = [phaseList string(PhasesToRead(i).name(1:end-4))]
end
%%
InputUser.ExportFolder = "C:\Users\berners\sciebo\Paper6_OccupancyAndEffectsInNbNi\muphase_stuff - Kopie\Structures\ParameterFiles\ILMeasure"
a =[]
c = []
sites = struct()
phaseName = []
for phaseNum = 1: size(phaseList,2)
InputUser.Phase_Input={phaseList{phaseNum}};
InputUser.Phase_Folder = Phase_Folder
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,Phase_Folder );
[cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 
a = [a abs(cs.aAxis)]
c = [c abs(cs.cAxis)]
% sites(phaseNum).phase = phaseList{phaseNum}
% sites(phaseNum).
end
%%
outputTable = table(phaseList',a',c','VariableNames',["phase","a","c"])
writetable(outputTable, fullfile(InputUser.ExportFolder,sprintf("LatticeParams.csv")),'WriteRowNames',true)
writetable(outputTable, fullfile(InputUser.ExportFolder,sprintf("LatticeParams.xlsx")),'WriteRowNames',true)

%%
for phaseNum = 1: size(phaseList,2)
InputUser.Phase_Input={phaseList{phaseNum}};
InputUser.Phase_Folder = Phase_Folder
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,Phase_Folder );
[cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 

keySet = {['Nb'],['Ni']};
color =  { [1 1 0] , [0 1 1] };
%%
% InputUser.ExportFolder = "I:\001OR_Analysis_Paper\AdditionsToFigures"
% InputUser.ExportFolder = "C:\Users\berners\sciebo\muphase_stuff - Kopie\Structures"
InputUser.FigureFolder = fullfile(InputUser.ExportFolder,InputUser.Phase_Input{1})
%%

h=Miller({1,0,0},{1,1,0},{1,1,1},cs);
h = Miller({1,1,-2,4},{0,2,-2,2},{2,2,-4,0},{0,1,-1,6},{0,0,1},cs)
miller_number=1
% figure
rot_mat=get_h_rot(h(miller_number));
ori=orientation(rot_mat,h.CS)
% miller_number=3
% figure
% rot_mat=get_h_rot(h(miller_number));

% ori=orientation(rot_mat,h.CS)
% plotIPDF(ori,vector3d.Z)

%%
%%

outputs=get_atom_positions(atoms,cs)
atoms_out = outputs.all_atoms;
% [C,ia,ic] = unique(output.all_atoms.xyz,'rows')
[~,ia,~] = unique(round(atoms_out.uvw,4),'rows');
atoms_sym=atoms_out(ia);
labels_sym = outputs.labels(ia);
siteNum = outputs.siteNum(ia);
[elements,~,elementIndex] = unique(outputs.elements(ia))
colorsNew = outputs.colors_elements(ia,:);
insideUnitCell = ~any(atoms_sym.uvw>=0.998,2);
hkl_unique = atoms_sym(insideUnitCell);
label_unique = labels_sym(insideUnitCell);
colorsNew = colorsNew(insideUnitCell,:);
siteNum = siteNum(insideUnitCell);
elementIndex = elementIndex(insideUnitCell);
%%
[uniqueLabels,indFind,~] = unique(label_unique(:))
%%
InterlayerSpacings = zeros(size(indFind,1),size(indFind,1))
for i = 1:size(indFind,1)
    for j = 1:size(indFind,1)
        InterlayerSpacings(i,j) = min(abs(hkl_unique(label_unique == uniqueLabels(j)).z - hkl_unique(indFind(i)).z))
    end
end

% outputTable = table(uniqueLabels,closestAtom',Distance,siteNumAtom,elementIndexLoc','VariableNames',{'AtomPosition','closest Atom','Distance','SiteNum Self','Site Num closest'})

%%
InterlayerSpacingsTable = array2table(InterlayerSpacings,'VariableNames',uniqueLabels,'RowNames',uniqueLabels)
writetable(InterlayerSpacingsTable, fullfile(InputUser.ExportFolder,sprintf("%s_ILSpacings.csv",InputUser.Phase_Input{1} )),'WriteRowNames',true)
writetable(InterlayerSpacingsTable, fullfile(InputUser.ExportFolder,sprintf("%s_ILSpacings.xlsx",InputUser.Phase_Input{1} )),'WriteRowNames',true)
end