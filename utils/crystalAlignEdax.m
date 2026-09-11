function csEdax = crystalAlignEdax(cs)
% this function converts the crystal symmetry to the edax convetion for hexagonal 
% systems, where a || x (compare a* || x for HKL)
    xyz = cs.axes.xyz;
    csEdax=crystalSymmetry(cs.pointGroup, ...
        [round(norm(xyz(1,:)),4) round(norm(xyz(2,:)),4) round(norm(xyz(3,:)),4)], ...
        'X||a','Z||c');
    csEdax.mineral = cs.mineral;
    %%

end
