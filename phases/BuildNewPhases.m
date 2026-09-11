% L. Berners, 15.08. 2025
% %Deck for building new phase files from .h5 files from EMsoft simulations
%
% The function parses information contained in our main simulation folder
% and generates new .phase files
% could likely be ajdusted to other simulation softwares


% For cluster we still need to do the folders manually
% InputUser.MTEX_loc='/work/mz071159/mtex';
% InputUser.MTEX_loc='C:\Users\berners\MTEX_dev\mtex';
try EBSD;
catch
    run(fullfile(InputUser.MTEX_loc,"startup.m"));
end

InputUser.Phase_Folder = Phase_Folder;
% InputUser.Phase_Folder="/work/mz071159/AstroEBSD_v2/phases";
% InputUser.Astro_loc='/work/mz071159/AstroEBSD_v2';
% InputUser.Phase_Folder=Phase_Folder;
InputUser.Astro_loc=Astro_loc;

%%


try astro_loadcheck
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end
%%
% clear -regexp ^(?!.*(Phase_Folder)|(Astro_loc)|MTEX_loc).

% tic 
% rootDir = '/rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/Simulation_NbNiStackingFaults'
% <<<<<<< HEAD
% % rootDir = '/rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/Simulation_YNi_replotDinah'
% =======
% rootDir = '/rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/Simulation_YNi_dbwTesting'
% >>>>>>> 484f988a3ccf923b64d54f4ee0da980af23b21c9

% rootDir = '/hpcwork/mz071159/EMsoft_work/NbNi'
%;_LatticeMeasured'
% phaseFolder  = "/work/mz071159/AstroEBSD_v2/phases"
rootDir = 'I:\002Simulation_EBSD\EMSoftHPC\SimulationNb_Ni_StackingFaults'
rootDir = "I:\002Simulation_EBSD\EMSoftHPC\Simulation_YNi_replotDinah"
rootDir = "I:\002Simulation_EBSD\EMSoftHPC\C15MgAlCa"
% rootDir = 'I:\002Simulation_EBSD\SimulationNb_Oxides'
tic
buildNewPhaseFiles(rootDir,InputUser.Phase_Folder)
toc

% check if .h5 files are not copy again the second time

% tic
% buildNewPhaseFiles(rootDir,phaseFolder)
% to
%%
allPhases = buildPhaseLibrary(InputUser.Phase_Folder);
%%
load(fullfile(InputUser.Phase_Folder,"PhaseLibrary.mat"),"-mat")
figure
%%

% histogram(PhaseLibrary.Voltage) 