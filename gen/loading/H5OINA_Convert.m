function [ebsdtemp] = H5OINA_Convert(ebsd_data,header_data,file_dset,phase_names,phase_colors,ebsd_nonpm)
%H5OINA_CONVERT Convert the output from the HFOINA_Read function into a
%MTEX EBSD container
%
% Ben Britton - Feb 2024
%
% Adapted from load_h5oina


slice_name=['s' file_dset];

CS{1}='notIndexed';
num_phases=numel(header_data.(slice_name).phase);

for n=1:num_phases
    phaseN=1+n;
    Space_Group=double(header_data.(['s' file_dset]).phase{n}.Space_Group);
    Lattice_Dimensions=double(header_data.(['s' file_dset]).phase{n}.Lattice_Dimensions);
    Lattice_Angles=double(header_data.(['s' file_dset]).phase{n}.Lattice_Angles);

    if nargin > 5 && ~isempty(phase_names)
        Mineral=phase_names{n};
        CS{phaseN} = crystalSymmetry('SpaceId',Space_Group, ...
            Lattice_Dimensions,...
            Lattice_Angles,...
            'Mineral',Mineral,'Color',phase_colors{n});
    else
        Mineral=char(header_data.(['s' file_dset]).phase{1}.Phase_Name);

        CS{phaseN} = crystalSymmetry('SpaceId',Space_Group, ...
            Lattice_Dimensions,...
            Lattice_Angles,...
            'Mineral',Mineral);
    end
end

rc = rotation.byEuler(double(header_data.(slice_name).Specimen_Orientation_Euler(:)')*degree); % what definition? Bunge?

% set up EBSD data
rot = rc*rotation.byEuler(ebsd_data.(slice_name).Euler');
phase = ebsd_data.(slice_name).Phase;
opt=struct;

% read some fields
EBSD_fieldnames=fieldnames(ebsd_data.(slice_name));

num_fields=size(EBSD_fieldnames,1);
for n = 1: num_fields
    s=ebsd_data.(slice_name).(EBSD_fieldnames{n});
    if size(s,2) == 1 && size(s,1) == numel(phase)
        try
            opt.(EBSD_fieldnames{n})=double(s);
        catch
            try
                opt.(EBSD_fieldnames{n})=s;
            catch
            end
        end

    end
end

%not put the Euler angles in the options too
if isfield(ebsd_data.(slice_name),'Euler')
    opt.euler1=double(ebsd_data.(slice_name).Euler(1,:)');
    opt.euler2=double(ebsd_data.(slice_name).Euler(2,:)');
    opt.euler3=double(ebsd_data.(slice_name).Euler(3,:)');
end

try
    opt.x=opt.Beam_Position_X;
    opt.y=opt.Beam_Position_Y;
catch
    try
        opt.x=ebsd_nonpm.prop.x;
        opt.y=ebsd_nonpm.prop.y;
    catch
        error('Beam Positions Not Loaded - was this a pattern match file?')
    end
end

ebsdtemp = EBSD(rot,phase,CS,opt,'unitCell',calcUnitCell([opt.x,opt.y]));
ebsdtemp.opt.Header = header_data.(slice_name);

end

