function [ebsdtemp,ebsd_header] = bIDX_to_EBSD(TFS_DataLoc)
%% bIDX_to_EBSD 
% Converts a TFS IDX file to an MTEX EBSD data set
% Created by T. Ben Britton in November 2024
% Edited and annotated by Tianbi Zhang in November 2024

% This script takes a Kikuchi pattern, index it, reproject it onto the
% diffraction sphere and create a stereogram.

% This script accompanies the following article(s):
%
% "Comparison of Kikuchi Diffraction Geometries in 
% Scanning Electron Microscope"
% Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton 

% Requirements: 
% (1) MATLAB toolboxes: image processsing, statistics and machine learning,
% parallel computing
% (2) AstroEBSD package - this script is available as a part of the latest
% AstroEBSD distribution.
% (3) MTEX toolbox (https://mtex-toolbox.github.io/)

multiphase_full=fullfile(TFS_DataLoc,'\results','multiphases.idx');
sitedata_full=fullfile(TFS_DataLoc,'site.json');

% Load the site data

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = [" ", ","];

% Specify column names and types
opts.VariableNames = ["VarName1", "Info", "VarName3", "VarName4"];
opts.VariableTypes = ["string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "Info", "VarName3", "VarName4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Info", "VarName3", "VarName4"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(sitedata_full, opts);

% Convert to output type
Info = tbl.Info;
VarName3 = tbl.VarName3;

ebsd_header=struct;

for n=1:size(Info,1)
    try
        fieldname=Info{n};
        %remove ':'
        s=strfind(fieldname,':');
        if numel(s) >=1
            fieldname=fieldname(1:s(1)-1);
            field_num=str2num(VarName3(n));
            if ~isempty(field_num)
            ebsd_header.(fieldname)=field_num;
            end

        end
    catch
    end
end

try
    ebsd_header.StepSize=ebsd_header.MapStepSize*1E6; %in um
catch
    ebsd_header.StepSize=1;
end

%%
% Load the map data from the idx stream

% Open the file for reading in binary mode
% Open the file for reading
fileID = fopen(multiphase_full, 'r');

% Read all lines into a cell array
lines = {};
tline = fgetl(fileID);
while ischar(tline)
    lines{end+1} = tline; % Append each line to the cell array
    tline = fgetl(fileID);
end

% Close the file
fclose(fileID);

lines_data=lines(1,5:end);

s=str2num(lines_data{1});
s_full=zeros(size(lines_data,2),size(s,2));

for n=1:size(lines_data,2)
s_full(n,:)=str2num(lines_data{n});
end

pts_num=str2num(lines{4});

%read the field data and build the options array
opt=struct;
opt_names=lines{:,2};
space_loc=strfind(opt_names,' ');
space_loc=[0 space_loc size(opt_names,2)+1];
%read the number of rows in the loaded data
num_fields=size(s_full,2);
for n=1:num_fields
    opt_text=opt_names((space_loc(n)+1):(space_loc(n+1)-1));
    opt.(opt_text)=s_full(:,n);

    % opt.(opt_text)=reshape(opt.(opt_text),pts_num(2),pts_num(1));
end

%scale the steps
try
    opt.x_in=opt.x;
    opt.y_in=opt.y;
    opt.x=opt.x*ebsd_header.StepSize;
    opt.y=opt.y*ebsd_header.StepSize;
catch
end

opt.y=-opt.y; %flip the Y axis because TFS do an axis switch

rot=rotation.byEuler([opt.Phi1 opt.PHI opt.Phi2]*degree);
phase=opt.PhaseID;

%% Import the phase data
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["entries", "values",];
opts.VariableTypes = ["string", "string",];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["entries", "values"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["entries", "values"], "EmptyFieldRule", "auto");

phaselist_full=fullfile(TFS_DataLoc,'phases\active.json');

Phase_List = readtable(phaselist_full, opts);

Phases_ok=find(strcmp(Phase_List.entries,'File:'));

CS='notIndexed';
phase_alldata=struct;

for n=1:numel(Phases_ok)
    phase_filename=Phase_List(Phases_ok(n),2).values;
    phase_filename_full=fullfile(TFS_DataLoc,'phases',phase_filename);

    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 3);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = [" ", ",", ":"];

    % Specify column names and types
    opts.VariableNames = ["name", "value"];
    opts.VariableTypes = ["string", "string"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";

    % Specify variable properties
    opts = setvaropts(opts, ["name", "value"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["name", "value"], "EmptyFieldRule", "auto");

    % Import the data
    phase_data = readtable(phase_filename_full, opts);
    
    for s=1:size(phase_data,1)
        try
            phase_alldata(n).(phase_data(s,:).name)=phase_data(s,:).value;
        catch

        end

    end
    cs_spacegroup=phase_alldata(n).space_group_HM;
    cs_abc=[str2num(phase_alldata(n).a) str2num(phase_alldata(n).b) str2num(phase_alldata(n).c)];
    cs_angs=[str2num(phase_alldata(n).alpha) str2num(phase_alldata(n).beta) str2num(phase_alldata(n).gamma)];
    cs_mineral=phase_alldata(n).name;

    CS2=crystalSymmetry(cs_spacegroup,cs_abc,cs_angs*degree,'mineral',cs_mineral);
    CS={CS CS2};
    
end

%% Clear temporary variables
clear opts

%%

ebsdtemp = EBSD(rot,phase,CS,opt,'unitCell',calcUnitCell([opt.y,opt.x]));

end

