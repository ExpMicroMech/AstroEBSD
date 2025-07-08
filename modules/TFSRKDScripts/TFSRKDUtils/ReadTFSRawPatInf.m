function [rowIndex, columnIndex] = ReadTFSRawPatInf(filename, dataLines)


%% Input handling

% If dataLines is not specified, define defaults
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["rowIndex", "columnIndex", "Var3", "Var4"];
opts.SelectedVariableNames = ["rowIndex", "columnIndex"];
opts.VariableTypes = ["double", "double", "string", "string"];

% Specify file level properties
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
rowIndex = tbl.rowIndex;
columnIndex = tbl.columnIndex;
end