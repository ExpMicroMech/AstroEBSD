function output=get_atom_positions(atoms,cs)
    


    elements=vertcat(atoms.atoms.atom_site_type_symbol);
    keySet_element=unique(string(elements));
    labels = []
    for i = 1:size(atoms.atoms,2)
        labels = [labels; convertCharsToStrings(atoms.atoms(i).atom_site_label)];
    end
    %%
    % labels=vertcat(atoms.atoms.atom_site_label);
    keySet_label=unique(string(labels));
    max_c_e=size(keySet_element,1);
    max_c_l=size(keySet_label,1);
    try
        M_label = containers.Map(keySet_label,get_colors_from_list(max_c_l));
    catch
        M_label = containers.Map(keySet_label,ones(max_c_l,1)*0.5); %just set it to gray
    end
    M_element=containers.Map(keySet_element,get_colors_from_list(max_c_e));
    % except
    % figure
    hold on
    all_atoms=[];
    atom_cord_symmetrised=[];
    labels=[];
    elements=[];
    colors_elements=[];
    colors_label=[];
    siteNum = [];
    %%
    for i=1:size(atoms.atoms,2)
        atom_cord=[atoms.atoms(i).atom_site_fract_x;atoms.atoms(i).atom_site_fract_y;atoms.atoms(i).atom_site_fract_z];
        for j=1:size(atoms.symops,2)
            atom_cord_sym=atoms.symops(j).matrix*atom_cord+atoms.symops(j).t_vec; %%https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html
%             if ~any(atom_cord_sym==1)
            atom_cord_sym=atom_cord_sym-floor(atom_cord_sym);%(everything gets to 0. sth)
            atoms_zero_checked=check_zeros(atom_cord_sym);
            for k=1:size(atoms_zero_checked,2)
                atom_vector=Miller(atoms_zero_checked(1,k),atoms_zero_checked(2,k),atoms_zero_checked(3,k),'UVW',cs);
                labels=[labels ;string(atoms.atoms(i).atom_site_label)];
                siteNum = [siteNum; i];
                elements=[elements ;string(atoms.atoms(i).atom_site_type_symbol)];

                mapped_color = M_element(atoms.atoms(i).atom_site_type_symbol);
                if iscell(mapped_color)
                    mapped_color=[mapped_color{:}];
                end
                colors_elements=[colors_elements; mapped_color];

                mapped_label =  M_label(atoms.atoms(i).atom_site_label);
                if iscell(mapped_label)
                    mapped_label=[mapped_label{:}];
                end
                colors_label=[colors_label; mapped_label];

                all_atoms=[all_atoms ;atom_vector];
            end
% %           
% 
% %             empties=double(atom_cord_sym==0)
%                 ctrl=size(find(atom_cord_sym==0))
%                 ind=find(atom_cord_sym==0)
%                 for i=1

%           

        end
         %%
    %     plot(atom_vector);
    %     scatter3(atom_cord(1),atom_cord(2),atom_cord(3),'black')
    %     hold on
    end
    output.labels=labels;
    output.elements=elements;
    output.colors_elements=colors_elements;
    output.colors_label=colors_label;
    output.all_atoms=all_atoms;
    output.siteNum = siteNum;
end

% 
% function [all_atoms,colors_label,colors_elements,elements,labels] =populate_edges(atom_cord_sym,atoms,i)
%      all_atoms=[];
% %     atom_cord_symmetrised=[];
%     labels=[];
%     elements=[];
%     colors_elements=[];
%     colors_label=[];
%     
% end