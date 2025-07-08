function dataType = ReadTFSRawPatHeader(filename, dataLines)

% this function returns the data type of raw TFS patterns.
% Which is useful for parsing the raw pattern files.

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, 3];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["CompressionMethod", "VarName2", "uint16", "VarName4"];
opts.VariableTypes = ["string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["CompressionMethod", "VarName2", "uint16", "VarName4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["CompressionMethod", "VarName2", "uint16", "VarName4"], "EmptyFieldRule", "auto");

% Import the data
info = readmatrix(filename, opts);
dataType = char(info(2,3))
end