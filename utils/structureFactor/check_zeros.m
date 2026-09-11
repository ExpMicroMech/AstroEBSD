function atoms_checked=check_zeros(atom_cords)
%% check those atoms that are close to zero
    if ~any(atom_cords<=0.002)
        atoms_checked=atom_cords;
    elseif sum(atom_cords<=0.002)==3
        atoms_checked=[[0;0;0;],[0;0;1],[0;1;0],[1;0;0],[1;1;0],[1;0;1],[0;1;1],[1;1;1]];
    elseif sum(atom_cords<=0.002)==2
        atoms_checked = [];
        index=find(atom_cords<=0.002);
        for i=0:1
            for j=0:1
                atom_cords_in=atom_cords;
                atom_cords_in(index(1))=i;
                atom_cords_in(index(2))=j;
                atoms_checked=[[atoms_checked], [atom_cords_in]];
            end
        end
    else
        atoms_checked=[];
        index=find(atom_cords<=0.002);
        for i=0:1
                atom_cords_in=atom_cords;
                atom_cords_in(index(1))=i;
                atoms_checked=[[atoms_checked], [atom_cords_in]];
        end

%         atoms_checked=[[0;0;0;],[0;0;1],[0;1;0],[1;0;0],[1;1;0],[1;0;1],[0;1;1],[1;1;1]]

    end
    %%
%     disp(atom_cords);
end