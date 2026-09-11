function outputs = getAtomicRadii(InputUser,pNum)
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(pNum),InputUser.Phase_Folder );
[cs,~,atoms]=loadCIF_atoms(RTM_info.cif_file); 

outputs=get_atom_positions(atoms,cs);
atoms_out = outputs.all_atoms;
% [C,ia,ic] = unique(output.all_atoms.xyz,'rows')
[~,ia,~] = unique(round(atoms_out.uvw,4),'rows');
atoms_sym=atoms_out(ia);
labels_sym = outputs.labels(ia);
elements_sym = outputs.elements(ia);
siteNum = outputs.siteNum(ia);
[elements,~,elementIndex] = unique(outputs.elements(ia))
colorsAll = outputs.colors_elements(ia,:);
insideUnitCell = ~any(atoms_sym.uvw>=0.998,2);
hkl_unique = atoms_sym(insideUnitCell);
label_unique = labels_sym(insideUnitCell);
element_unique = elements_sym(insideUnitCell);
colorsNew = colorsAll(insideUnitCell,:);
siteNum = siteNum(insideUnitCell);
elementIndex = elementIndex(insideUnitCell);
%%
xMin = -7;
xMax = 7;
yMin = -7;
yMax = 7;
zMin = -7;
zMax = 7;
[XP, YP, ZP] = meshgrid(xMin:1:xMax,yMin:1:yMax,zMin:1:zMax);
replicator =Miller(XP(:),YP(:),ZP(:),'UVW',cs);
%
hklNew = [];
for i = 1:size(replicator)
    hklNew  = [hklNew, hkl_unique+replicator(i)];
end
colorsRep = repmat(colorsNew,size(replicator));
labelsRep = repmat(label_unique,size(replicator));
elementRep = repmat(element_unique,size(replicator));
% radiiRep
siteNumRep = repmat(siteNum,size(replicator));
%
elementIndexRep = repmat(elementIndex, size(replicator));
% figure
% scatter3(hklNew(:).x,hklNew(:).y,hklNew(:).z);
% newVecs = hkl_unique+replicator(1);
[uniqueLabels,indFind,~] = unique(label_unique(:));
siteNumUnique = siteNum(indFind);
NeighbourDistance = [];
Distance = [];
closestAtom =[];
siteNumAtom =[];
elementIndexLoc=[];
NeighbourLabels = [];
NeighbourSites = [];
hklRep = [];
%%
for i = 1: size(indFind)
    neighbourDists = abs(hkl_unique(indFind(i))- hklNew(:));
    NeighbourDistance = [NeighbourDistance, neighbourDists];
    NeighbourLabels = [NeighbourLabels, labelsRep];
    NeighbourSites = [NeighbourSites, siteNumRep];
    hklRep = [hklRep, hklNew(:)];
    labelsDiff = labelsRep(neighbourDists~=0);
    siteNumDiff = siteNumRep(neighbourDists~=0);
    elementIndexDiff = elementIndexRep(neighbourDists~=0);
    [minVal,minInd] = min(neighbourDists(neighbourDists~=0));

    Distance = [Distance, minVal];
    closestAtom = [closestAtom, labelsDiff(minInd)];
    siteNumAtom = [siteNumAtom, elementIndexRep(minInd)];
    elementIndexLoc = [elementIndexLoc, elementIndex(indFind(i))];
    % Between = [label_unique] 
end
siteNumAtom = siteNumAtom';
Distance = Distance';
%%

    

%%
LGS = zeros(size(uniqueLabels,1),size(uniqueLabels,1));
for i =1:size(uniqueLabels,1)
    for j = 1:size(uniqueLabels,1)
    LGS(i,j)= sum([siteNumAtom(i)==j elementIndexLoc(i)==j]);
    end
end
LGS;
%%
% C = 0.2*ones(size(LGS)) + 0.8*diag(ones(size(LGS)));
% [x,stdx,mse] = lscov(LGS,Distance,C)
%%
dec = decomposition(LGS,'cod');
radiussolution = lsqminnorm(LGS,Distance);%,'warn','RegularizationFactor',1)
%%
pinv(LGS)*Distance

outputTable = table(uniqueLabels,closestAtom',Distance,siteNumAtom,elementIndexLoc','VariableNames',{'AtomPosition','closest Atom','Distance','SiteNum Self','Site Num closest'})
writetable(outputTable, fullfile(InputUser.FigureFolder , ...
                                sprintf('%s_SiteDistances.xls',InputUser.Phase_Input{pNum}) ) ); 
writetable(outputTable, fullfile(InputUser.FigureFolder , ...
                            sprintf('%s_AtomiRadii.csv',InputUser.Phase_Input{pNum}) ) ); 

% uniqueLabels
%%
radii = zeros(size(outputs.labels));
uniqueLabels_1 = unique(outputs.labels);
for i = 1:size(uniqueLabels_1,1)
    radii(strcmp(outputs.labels,uniqueLabels_1(i)))=radiussolution(i);
end

radiiNew = zeros(size(elementRep));
uniqueLabels_2 = unique(elementRep);
for i = 1:size(uniqueLabels_2,1)
    radiiNew(strcmp(elementRep, ...
        uniqueLabels_2(i)))=radiussolution(i);
end

radiiInside = zeros(size(element_unique));
uniqueLabels_3 = unique(element_unique);
for i = 1:size(uniqueLabels_3,1)
    radiiInside(strcmp(element_unique, ...
        uniqueLabels_3(i)))=radiussolution(i);
end

radiiUnique = zeros(size(uniqueLabels));
uniqueLabels_4 = unique(uniqueLabels);
for i = 1:size(uniqueLabels_4,1)
    radiiUnique(strcmp(uniqueLabels, ...
        uniqueLabels_4(i)))=radiussolution(elementIndexLoc(i));
end
%%
[uniqueElements,firstUniqueLabel,~] = unique(outputs.elements);
for i = 1:size(uniqueElements,1)
    radii(strcmp(outputs.elements,uniqueElements(i)))=radiussolution(i);
end
% radii()
outputs.hklNew = hklNew(:);
outputs.radii = radii;
outputs.radiiNew = radiiNew;
outputs.uniqueLabels = uniqueLabels;
outputs.uniqueElements = elementRep;
outputs.atomsUnique = hkl_unique(indFind);
outputs.NeighbourDistance = NeighbourDistance;
outputs.NeighbourLabels = NeighbourLabels;
outputs.NeighbourSites = NeighbourSites;
outputs.colorsNew = colorsRep;
outputs.outputTableClosestNeighbour = outputTable;
outputs.hklRep = hklRep;


outputs.labelsInsideUnitCell = label_unique;
outputs.elementsInsideUnitCell = element_unique;
outputs.atomsinsideUnitCell = hkl_unique;
outputs.radiiUnique = radiiUnique;

% find(unique(outputs.radii)
outputTable = table(outputs.elements(firstUniqueLabel), radii(firstUniqueLabel),...
            'VariableNames',{'Atom','Radius'});
writetable(outputTable, fullfile(InputUser.FigureFolder , ...
                                sprintf('%s_AtomiRadii.xls',InputUser.Phase_Input{pNum}) ) ); 
writetable(outputTable, fullfile(InputUser.FigureFolder , ...
                                sprintf('%s_AtomiRadii.csv',InputUser.Phase_Input{pNum}) ) );
outputs.outputTable = outputTable;
end