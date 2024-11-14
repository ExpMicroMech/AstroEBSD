function [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( Phases,Phase_Folder )
%PHASE_BUILD Build all the phases needed for the EBSD data
%Phases = a structure of phase file names
%Phase_Folder = folder where they exist

%% Versioning
%v1 - TBB 14/04/2017

% v2 - Alex Foden modified the code to work with RTM phase files
% Do not distribute. 09/04/2019

%v3 - TPM modified for integation into AstroEBSD for use with RTM and PCA
%analysis. Generalised and .pha no longer requires specified folder -
%defaults to Phase folder.

%v4 (tentative) - Tianbi Zhang modified the code to work with less
% symmetric structures by using directly the loadCIF function in MTEX
% 08/31/2023

%v5 - Lukas Berners modified the code, to load in the information about
%EMosft 'Masterpatterns' from the .pha files as well. If no information about
%the simulation type is found, it will be assumed, that a Bruker Dynamics
%.bin file should be loaded
% 2023_07_19 (JJJJ_MM_DD) 

%v6 - TBB modified the code to be able to load bwkd sterographic
%projections, as created from the bwkk python package / using the
%MapSweeper simulation engine

%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above

%% run code
num_Phases=size(Phases,2);

Phase_loc=fullfile(Phase_Folder,'phasefiles');
Bin_loc=fullfile(Phase_Folder,'dynamic_templates');
Cif_loc=fullfile(Phase_Folder,'cifs');

for num_P=1:num_Phases
    %read the phases
    Phase_File=fullfile(Phase_loc,[Phases{num_P} '.pha']);
    if ~exist(Phase_File) %#ok<EXIST>
        error(['This phase file does not exist: ' Phase_File])
    end

%% RTM extra information needed
delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(Phase_File,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
PhaseRaw = [dataArray{1:end-1}];

array_num=1:size(PhaseRaw,1);

%Extract the cif file
ix=array_num(contains(PhaseRaw,'$cif','IgnoreCase',true));
RTM_info.cif_file=fullfile(Cif_loc,PhaseRaw{ix+1});

%Extract the bin file name
 try
            %%
            ix=array_num(contains(PhaseRaw,'$dynamic','IgnoreCase',true));
            bin_file_name=PhaseRaw{ix+1};
            %%
            ix=array_num(contains(PhaseRaw,'$sim_type','IgnoreCase',true));
            if ~isempty(ix)
                simtype=PhaseRaw{ix+1};
            else
                simtype='default';
            end
            %%
            %%
            if strcmpi(simtype,'Emsoft')
                [s]=strfind(bin_file_name,'.h5');  
                if ~isempty(s)
                    bin_file_name=bin_file_name(1:s-1);
                end
                RTM_info.bin_file = fullfile(Bin_loc, [bin_file_name, '.h5']);
            elseif strcmpi(simtype,'bwkd')
                [s]=strfind(bin_file_name,'.bwkd');
                if ~isempty(s)
                    bin_file_name=bin_file_name(1:s-1);
                end
                RTM_info.bin_file = fullfile(Bin_loc, [bin_file_name, '.bwkd']);
            elseif strcmpi(simtype,'MapSweeper')
                [s]=strfind(bin_file_name,'.sdf5');  
                if ~isempty(s)
                    bin_file_name=bin_file_name(1:s-1);
                end
                RTM_info.bin_file = fullfile(Bin_loc, [bin_file_name, '.sdf5']);

            else
    
                [s]=strfind(bin_file_name,'.bin');
                if ~isempty(s)
                    bin_file_name=bin_file_name(1:s-1);
                end
                RTM_info.bin_file = fullfile(Bin_loc, [bin_file_name, '.bin']);
            end
            ix=array_num(contains(PhaseRaw,'$isHex','IgnoreCase',true));
 catch
     disp('We could not identify your Dynamic simulation :(')
end

%generate the unit cells
% This is a new version using loadCIF from MTEX toolbox.
[Crystal_UCell{num_P}]=Build_UCell(Phase_File, RTM_info.cif_file); %#ok<*AGROW>

% The older version of the function callout is the following:
% [Crystal_UCell{num_P}]=Build_UCell(Phase_File); %#ok<*AGROW>

%generate the plane families
[Crystal_Family{num_P} ] = Build_Reflec( Crystal_UCell{num_P} );
%generate the look up tables (LUTs)
[Crystal_LUT{num_P}] = Build_LUT(Crystal_Family{num_P});
%threshold for the LUT generation
Settings_LUT{num_P}.thresh_trig=1/Crystal_UCell{num_P}.efac;

%Extract the ishex information
ix=array_num(contains(PhaseRaw,'$isHex','IgnoreCase',true));
RTM_info.isHex=PhaseRaw{ix+1};
if ischar(RTM_info.isHex)
    RTM_info.isHex=str2double(RTM_info.isHex);
end

end
end

% function [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( Phases,Phase_Folder )
% %PHASE_BUILD Build all the phases needed for the EBSD data
% %Phases = a structure of phase file names
% %Phase_Folder = folder where they exist
% 
% %% Versioning
% %v1 - TBB 14/04/2017
% 
% % v2 - Alex Foden modified the code to work with RTM phase files
% % Do not distribute. 09/04/2019
% 
% %v3 - TPM modified for integation into AstroEBSD for use with RTM and PCA
% %analysis. Generalised and .pha no longer requires specified folder -
% %defaults to Phase folder.
% 
% %
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% % EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% % OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% % NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% % BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% % ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% % CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
% % Requirements:
% % MATLAB R2018a or above
% % MTEX version 5.2.beta2 or above
% 
% %% run code
% num_Phases=size(Phases,2);
% 
% Phase_loc=fullfile(Phase_Folder,'phasefiles');
% Bin_loc=fullfile(Phase_Folder,'masterpatterns');
% Cif_loc=fullfile(Phase_Folder,'cifs');
% 
% for num_P=1:num_Phases
%     %read the phases
%     Phase_File=fullfile(Phase_loc,[Phases{num_P} '.pha']);
%     if ~exist(Phase_File) %#ok<EXIST>
%         error(['This phase file does not exist: ' Phase_File])
%     end
%     
%     %generate the unit cells
%     [Crystal_UCell{num_P}]=Build_UCell(Phase_File); %#ok<*AGROW>
%     %generate the plane families
%     [Crystal_Family{num_P} ] = Build_Reflec( Crystal_UCell{num_P} );
%     %generate the look up tables (LUTs)
%     [Crystal_LUT{num_P}] = Build_LUT(Crystal_Family{num_P});
%     %threshold for the LUT generation
%     Settings_LUT{num_P}.thresh_trig=1/Crystal_UCell{num_P}.efac;
% 
% 
% %% RTM extra information needed
% delimiter = '';
% formatSpec = '%s%[^\n\r]';
% fileID = fopen(Phase_File,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% PhaseRaw = [dataArray{1:end-1}];
% 
% array_num=1:size(PhaseRaw,1);
% 
% %Extract the cif file
% ix=array_num(contains(PhaseRaw,'$cif','IgnoreCase',true));
% RTM_info.cif_file=fullfile(Cif_loc,PhaseRaw{ix+1});
% 
% %Extract the bin file name
% ix=array_num(contains(PhaseRaw,'$dynamic','IgnoreCase',true));
% bin_file_name=PhaseRaw{ix+1};
% 
% [s]=strfind(bin_file_name,'.bin');
% if ~isempty(s)
%     bin_file_name=bin_file_name(1:s-1);
% end
% RTM_info.bin_file = fullfile(Bin_loc, [bin_file_name, '.bin']);
% 
% 
% %Extract the ishex information
% ix=array_num(contains(PhaseRaw,'$isHex','IgnoreCase',true));
% RTM_info.isHex=PhaseRaw{ix+1};
% if ischar(RTM_info.isHex)
%     RTM_info.isHex=str2double(RTM_info.isHex);
% end
% 
% end
% end
% 