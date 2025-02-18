function [ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location,astro_location)
% LOAD_H5OINA_PM Loads h5oina data, and pattern matched data
% pattern matched data should be in the same folder and with the file name of the original
% but an additional '-PatternMatching' appended before the
%
% [ebsd_original,dataset_header,ebsd_patternmatched] = load_h5oina_pm(fname_slice_data,fname_folder,mtex_location)
%
% INPUTS
% fname_slice_data = h5oina name of 1 slice,
% excluding the folder and h5oina extension
% fname_folder = folder location of the h5oina file
% mtex_location = location of mtex 5.8.1
%
% OUTPUTS
% ebsd_original container of the original data
% ebsd_patternmatched = patternmatched container, if the pattern matching exists
% [if this does not exist, then it will just be a copy of the original]
% dataset_header = header data from the ebsd file
%
% optional
% h5_original = h5 original file path
% h5_patternmatch = h5_pattern match file path

try EBSD;
catch
    run(fullfile(mtex_location,"startup.m"));
end

%start astro if needed
try astro_loadcheck
catch
    run(fullfile(astro_location,"start_AstroEBSD.m"));
end

%do not change unless you really need to
fset_match_suffix='-PatternMatching';
fset_fileformat='.h5oina';

%h5 full files
h5_original=fullfile(fname_folder,[fname_slice_data fset_fileformat]);
h5_patternmatch=fullfile(fname_folder,[fname_slice_data fset_match_suffix fset_fileformat]);

%check file exists
if exist(h5_original,'file') ~= 2
    % error([h5_original ' can not be found'])
end

if exist(h5_patternmatch,'file') == 2
    pm_data=1;
else
    warning([h5_patternmatch ' can not be found - ignore if data not pattern matched'])
    pm_data=0;
end

%load the data
% Load data - pattern matched
% [ebsd_o,~,opt_o] = loadEBSD_h5oina_getopt(h5_original);
[ebsd_o,EBSDdata,opt_o,ebsdpm2,EBSDdata2]= loadEBSD_h5oina_getopt(h5_original);

if ~isempty(fieldnames(ebsdpm2))
    ebsd_pm=ebsdpm2;
elseif pm_data == 1
    [ebsd_pm,~,~,~] = loadEBSD_h5oina_PatMatch(h5_patternmatch,ebsd_o,opt_o) ;
end

% create header
dataset_header=ebsd_o.opt.Header;

ebsd_original=ebsd_o;

if pm_data == 1
    ebsd_patternmatched=ebsd_pm;
elseif ~isempty(fieldnames(ebsdpm2))
    ebsd_patternmatched=ebsdpm2;
else
    ebsd_patternmatched=ebsd_original;
end

end
