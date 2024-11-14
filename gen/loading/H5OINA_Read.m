function [ebsd_data,header_data,h5_contents]=H5OINA_Read(file1_full)
%read the H5OINA file for EBSD data - will return header and EBSD_type data
%takes a filename as the input and returns the contents (or most of it) at the end
%due to memory requirements, will not read patterns to memory
%
% Ben Britton - Feb 2024
%
% Change Log
% March 2024 - updated to read pattern matched files



%read the contents of the file
[h5_contents] = get_h5_list(file1_full);

%go through each dataset for the structure and read the h5 data
num_datasets=numel(h5_contents.datasets);

ebsd_data=struct;
header_data=struct;

for n=1:num_datasets
    %see if it is 'EBSD/Data'
    data_str='/EBSD/Data/';
    data_srt_len=numel(data_str);
    data_setname=h5_contents.datasets{n};

    data_loc=strfind(data_setname,data_str);
    if ~isempty(data_loc)
        if isempty(strfind(data_setname,'Patterns')) %avoid loading processed or unprocessed patterns
            %read the pre data number
            data_str_pre=[];
            data_str_pre=data_setname(1:data_loc-1);
            data_str_pre=['s' data_str_pre]; %add a letter to start, because matlab does not like numerical variable names

            data_str_end=data_setname(data_loc+data_srt_len:end);
            %convert any spaces to underscore
            data_str_space=strfind(data_str_end,' ');
            data_str_end(data_str_space)='_';


            ebsd_data.(data_str_pre).(data_str_end)=h5read(file1_full,['/' data_setname]);
        end
    end


    %see if it is 'EBSD/Header'
    data_str='/EBSD/Header/';
    data_srt_len=numel(data_str);
    data_setname=h5_contents.datasets{n};

    data_loc=strfind(data_setname,data_str);
    if ~isempty(data_loc)

        %read the pre data number
        data_str_pre=[];
        data_str_pre=data_setname(1:data_loc-1);
        data_str_pre=['s' data_str_pre]; %add a letter to start, because matlab does not like numerical variable names

        data_str_end=data_setname(data_loc+data_srt_len:end);

        %check this does not have any '/' subsets
        endian=strfind(data_str_end,'/');
        if isempty(endian)
            %convert any spaces to underscore
            data_str_space=strfind(data_str_end,' ');
            data_str_end(data_str_space)='_';

            %write the header data
            header_data.(data_str_pre).(data_str_end)=h5read(file1_full,['/' data_setname]);
        end
    end

    %read phase data
    data_str='/Header/Phases/';
    data_srt_len=numel(data_str);
    data_setname=h5_contents.datasets{n};

    data_loc=strfind(data_setname,data_str);
    if ~isempty(data_loc)

        %read the pre data number
        data_str_pre=[];
        data_str_pre=data_setname(1:data_loc-1);
        data_str_pre=['s' data_str_pre]; %add a letter to start, because matlab does not like numerical variable names
        data_str_pre=data_str_pre(1:strfind(data_str_pre,'/')-1);

        data_str_end=data_setname(data_loc+data_srt_len:end);
        ploc=strfind(data_str_end,'/');

        phasenum=str2double(data_str_end(1:ploc-1));
        fieldname=data_str_end(ploc+1:end);
        %convert any spaces to underscore
        fieldname(strfind(fieldname,' '))='_';
        
        %write the phase data for this slice
        try
        header_data.(data_str_pre).phase{phasenum}.(fieldname)=h5read(file1_full,['/' data_setname]);
        catch
%             disp([data_str_pre ' phase ' int2str(phasenum) ' ' fieldname 'not loaded']);
        end

    end
end

end

function [sOut] = get_h5_list(file1_full)
%read the h5 structure and return all the group data
%based upon https://www.mathworks.com/matlabcentral/answers/1898315-list-all-datasets-and-groups-of-an-hdf5-file
%
% Ben Britton - Feb 2024

s = struct;
s.datasets = {};
s.groups = {};
fid = H5F.open(file1_full, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
[status sOut] = H5O.visit(fid, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', @getObjInfo, s);

    function [status, sOut] = getObjInfo(objId, name, sIn)

        objID = H5O.open(objId, name, 'H5P_DEFAULT');
        obj_info=H5O.get_info(objID);
        H5O.close(objID);

        switch(obj_info.type)
            case H5ML.get_constant_value('H5G_LINK')
                % fprintf('Object #%s is a link.\n', name);
            case H5ML.get_constant_value('H5G_GROUP')
                % fprintf('Object #%s is a group.\n', name);
                sIn.groups{end+1,1} = name;
            case H5ML.get_constant_value('H5G_DATASET')
                % fprintf('Object #%s is a dataset.\n', name);
                sIn.datasets{end+1,1} = name;
            case H5ML.get_constant_value('H5G_TYPE')
                % fprintf('Object #%s is a named datatype.\n', name);
        end
        status = 0; % keep iterating
        sOut = sIn;
    end

end