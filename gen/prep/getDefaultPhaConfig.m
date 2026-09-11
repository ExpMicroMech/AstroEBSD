function [keys,vals]= getDefaultPhaConfig()
    keys = ["$Name", "$Atoms", "$Symmetry", "$Lattice", "$FFactor", "$MaxFam", "$MaxIndex", "$Efac","$sim_type"];
    vals = ["template","0,0,0","2","4,4,4,90,90,90","2,2,2,2","8","20","100","Emsoft"];
end