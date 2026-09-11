function scan = load_edax_binary(filename, nav_shape, lazy)
    % Read EBSD patterns from an EDAX binary UP1/2 file.
    %
    % Parameters
    % ----------
    % filename : str
    %     File path to UP1/2 file with "up1" or "up2" extension.
    % nav_shape : tuple
    %     Navigation shape, as (n map rows, n map columns), of the
    %     returned signal, matching the number of patterns in the file. If
    %     not given, this shape will be attempted to be determined from the
    %     file. If it could not be, the returned signal will have only one
    %     navigation dimension. If patterns were acquired in an hexagonal
    %     grid, the returned signal will have only one navigation dimension
    %     irrespective of this parameter's value.
    % lazy : bool
    %     Read the data lazily without actually reading the data from disk
    %     until required. Default is False.
    %
    % Returns
    % -------
    % scan : dict
    %     Data, axes, metadata and original metadata.
    
    % Open the file
    fid = fopen(filename, 'rb');
%     scan=[]
    [~, ~, ext] = fileparts(filename);
    ext = lower(ext(2:end));
    info.filename=filename
%     .file = file

%     raise ValueError("Only files with version 1 or >= 3, not 2, can be read")
       
    if strcmp(ext, 'up1')
        info.dtype = 'uint8';
    elseif strcmp(ext, 'up2')
        info.dtype = 'uint16';
    end
    fseek(fid, 0, 'bof');
    info.version = fread(fid, 1, 'uint32');
    if info.version == 2
        error('Only files with version 1 or >= 3, not 2, can be read');
    end





    scan=read_scan(info)

        
    function header = read_header(info)
        % Read and return header information.
        disp(':)')
        fid = fopen(info.filename, 'rb');
        %%
        fseek(fid, 4, 'bof');
%         [sx, sy, pattern_offset] = 
        ans=fread(fid, 3, 'uint32',0)
        sx=ans(1)
        sy=ans(2)
        pattern_offset=ans(3)
        %%
        file_info = dir(info.filename);
        file_size = file_info.bytes;
        if info.version == 1
            n_patterns = (file_size - pattern_offset) / (sx * sy * self.dtype(0).nbytes);
            n_patterns = int32(n_patterns);
            nx = n_patterns;
            ny = 1;
            dx = 1;
            dy = 1;
            is_hex = false;
        else
            %%
            pos=ftell(fid)

%%
            fseek(fid,pos +1,'bof');
%             ftell
            %%
%             junk=fread(fid,20)
            nx = fread(fid, 1, 'uint32', 0);
            ny = fread(fid, 1, 'uint32', 0);
%             ny2 = fread(fid, 1, 'uint32', 1);
            %%
            is_hex = fread(fid, 1, 'uint8', 1);
%%

            is_hex = logical(is_hex);
            if is_hex
                warning('Returned signal has one navigation dimension since an hexagonal grid is not supported');
                n_patterns = (file_size - pattern_offset) / (sx * sy * info.dtype(0).nbytes);
                n_patterns = int32(n_patterns);
                nx = n_patterns;
                ny = 1;
            else
                n_patterns = nx * ny;
            end
            %%
            fseek(fid,pos +10,'bof')
            dx = fread(fid, 1, 'float64', 0);
            dy = fread(fid, 1, 'float64', 0);
        end
        header = struct('sx', sx, 'sy', sy, 'pattern_offset', pattern_offset, ...
            'nx', nx, 'ny', ny, 'n_patterns', n_patterns, 'dx', dx, 'dy', dy, 'is_hex', is_hex);
    end
    
    function scan = read_scan(info)
        % Return a dictionary with scan information and patterns.
        header = read_header(info);
        fid=fopen(info.filename, 'rb')
        if ~exist('nav_shape', 'var') || isempty(nav_shape) || isempty(header)
            nav_shape = [];
        end
        if ~exist('lazy', 'var') || isempty(lazy)
            lazy = false;
        end
%         if ~header.is_hex
%             ny = header.ny;
%             nx = header.nx;
%             n_patterns = int32(ny * nx);
% %             if n_patterns ~= header.n_patterns
% %                 error(['Given `nav_shape` ', num2str(nav_shape), ' does not match the number of patterns in the file, ', num2str(header.n_patterns), '.']);
% %             end
%         else
            ny = header.ny;
            nx = header.nx;
%         end
        sy = header.sy;
        sx = header.sx;
        data_shape = [sx,sy,nx,ny];
        if ny == 1
            data_shape = data_shape(2:end);
        end
        ndim = length(data_shape);
        nav_dim = ndim - 2;
        fseek(fid, header.pattern_offset, 'bof');
        data = fread(fid, prod(data_shape), [info.dtype '=>' info.dtype]);
        data = reshape(data, data_shape);
        if lazy
            data = reshape(data, data_shape);
        end
        data=permute(data,[2 1 3 4]);
        units = repmat({'um'}, [1, ndim]);
        scales = [header.dy, header.dx, 1, 1];
        axes_names = {'y','x','dy','dx'};
        %%
        axes_names(end-nav_dim+1:end)
        %%
        zeros(1, ndim)
%%
        axes = struct('size', num2cell(data_shape), 'index_in_array', num2cell(1:ndim), ...
            'name', axes_names, 'scale', num2cell(scales), ...
            'offset', num2cell(zeros(1, ndim)), 'units', units);
%         [~, fname, ~] = fileparts(fid);
            [~,fname,~]=fileparts(filename)
        metadata = struct('General', struct('original_filename', fname, 'title', fname), ...
            'Signal', struct('signal_type', 'EBSD', 'record_by', 'image'));
        scan = struct('axes', axes, 'data', data, 'metadata', metadata);
    end
end
% end

