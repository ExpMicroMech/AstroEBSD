function [data1,table_data]=header_read_tescan(ftif,variable)

try
if numel(ftif) == 1
    ftif=ftif{1};
end

dot_loc=strfind(ftif,'.');
ftif2=ftif;
ftif2(dot_loc(end))='-';
filename=[ftif2 '.hdr'];

opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "=";

% Specify column names and types
opts.VariableNames = ["MAIN", "VarName2"];
opts.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "MAIN", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "MAIN", "EmptyFieldRule", "auto");


% Import the data
table_data = readtable(filename, opts);

num_var=size(variable,2);
var_row=zeros(num_var,1);
data1=cell(num_var,1);
for n=1:num_var
    try
        var_row(n)=find(table_data.MAIN == variable(:,n));
        data1{n}=table_data.VarName2(var_row(n));
    catch
    end
end

%% Clear temporary variables
clear opts

catch
    data1 = 0;
    table_data = 0;
end

end