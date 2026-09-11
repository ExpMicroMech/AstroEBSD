function [densities, hkls_unique, atomPosition] = getClosePackedDirections(atoms,InputUser,varargin)
%%
cutoff = 0.85;
for i = 1:size(varargin,2)
    if strcmp(varargin{i}, 'cutoff')
        cutoff = varargin{i+1}
    end
end
% atoms = outputs{pNum}.atomicRadiiVars;
% hklNew = repmat(atoms.hklNew,5)
minRadius = min(atoms.outputTable.Radius);
maxRadius = max(atoms.outputTable.Radius);
[neighbour,origin] = find(atoms.NeighbourDistance<(2*maxRadius*1.4) & atoms.NeighbourDistance>0);
ind = sub2ind(size(atoms.NeighbourDistance),neighbour,origin);
neighbours = atoms.NeighbourLabels(ind);
dists = atoms.NeighbourDistance(ind);
hkls = atoms.hklRep(ind)
origins = atoms.NeighbourLabels(origin)
% origins = atoms.
% dSpacings = outputs{pNum}.dSpacing;
% dSpacings(dSpacings>0.7*minRadius);
hkls_unique = unique(atoms.atomsUnique(origin)-hkls,'stable')
% [~,ind_hkl_unique,~] = unique([hkls_unique.x,hkls_unique.y,hkls_unique.z],'rows');
% hkls_unique = hkls_unique(ind_hkl_unique);
% abs(hkls_unique/(2*maxRadius))
% hkls_unique = hkls_unique(abs(hkls_unique)<2*maxRadius)
densities =[]
for i = 1:size(hkls_unique,1)
    densities = [densities plotAtomsAlongDirection(atoms,hkls_unique(i),neighbour,origin,InputUser,varargin{:})'];
end
%%
% uniqueLabels = atoms.uniqueLabels;
[densRow,densCol] = find(densities>cutoff);
densities = densities(densities>cutoff);
atomPosition = atoms.uniqueLabels(densRow);
hkls_unique = hkls_unique(densCol);
end