function [file1_array,set_proj] = read_sim(simulation_file,pattern_file)
%READ_SIM Summary of this function goes here
%   Detailed explanation goes here

%read the simulation info in
[sim1_data]=read_text(simulation_file);

%read some useful parameter
set_proj=find_val(sim1_data,'projection_mode');

if set_proj == 0 %stereogram
    % 0: double-hemisphere stereographic master projection, ncols is the number of pixels along each hemisphere dimension,
    % this will return 2*ncols x ncols array with the two hemispheres (upper and lower) arranged vertically
    % The x,y coordinates of both hemispheres agree, only z is +z (upper) and -z (lower),
    % i.e. this means the backside of the lower hemisphere is seen from the top
    set_yhig=find_val(sim1_data,'pattern_ncols');
    set_xwid=set_yhig*2;
    set_yhig=set_yhig;
else
    set_yhig=find_val(sim1_data,'pattern_ncols');
    set_xwid=find_val(sim1_data,'pattern_nrows');
end

%read the pattern data
[file1_array]=read_nums(pattern_file,set_xwid,set_yhig);
end

function val_out=find_val(string_cell,text_to_find)
    % try
    %find the element you want to read
    el_found1=arrayfun(@(x) strfind(x,text_to_find),string_cell);
    el_found2=logical(1-cellfun(@isempty,el_found1));
    ret_string=string_cell{el_found2};
    %extract value
    colon_pos=strfind(ret_string,':');
    comma_pos=strfind(ret_string,',');
    val_out=str2double(ret_string(colon_pos+1:comma_pos-1));
    % catch
        % val_out=NaN;
    % end

end

function [sim1_data]=read_text(simulation_file)
    %read the file into a cell
    sim1_fid=fopen(simulation_file);
    sim1_data=textscan(sim1_fid,'%s','Delimiter','\n');
    fclose(sim1_fid);
    sim1_data=sim1_data{1};
end

function [file1_array]=read_nums(pattern_file,xwid,yhig)
    %read the pattern from the file
    file1_fid=fopen(pattern_file);
    file1_data=textscan(file1_fid,'%f ',xwid*yhig);
    fclose(file1_fid);
    file1_array=rot90(reshape(file1_data{1},yhig,xwid));
end