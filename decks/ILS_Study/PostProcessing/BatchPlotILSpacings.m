%% L. Berners 27.01. 2026.
% This script is supposed to run the PlotRefinementResutlts_staggered
% skripts more automatically for all analysis files
%

%%



%% Define Paths that the skript is supposed to run on.
rootpahtList = [%"NbNi_septemberbatch\51Nb-49Ni_v2\Area10", ...
                 "NbNi_septemberbatch\58Nb-42Ni_v2\Area7" ,...               
                    % "NbCoUBC\48Nb52Co_verticalCut\Area2",...
                 % "NbCoUBC\56Nb44Co_verticalCut\Area3",...
                % "NbCoUBC\49Nb51Co_WeiLuo_HT2\Area2" ...
                %];
]
rootrootpath = "J:\Users\Berners_Lukas\clara";


                % r"NbCoUBC\49Nb51Co_WeiLuo_HT2\Area2"]   a
mapFilterString = "Map*";               
% structureFilterString = "*HighBethe";
structureFilterString = "*3aILGrid*_01";% Change To :*3aILGrid* for IL happenings  % here you can specify how the "type" of structure is named.
additionalFilterString = "*EBSD+EDS Map*"

for rPNum =1 : numel(rootpahtList)
    mapDirs = dir(fullfile(rootrootpath,rootpahtList(rPNum), mapFilterString));
   
    for mapPNum = 1:numel(mapDirs)
        strucutre_dir = dir(fullfile(rootrootpath,rootpahtList(rPNum),mapDirs(mapPNum).name, structureFilterString));
        for structPNum = 1: numel(strucutre_dir)
            %%
            clean_dir = fullfile(strucutre_dir(structPNum).folder, strucutre_dir(structPNum).name, "cleaned");
            h5File = dir(fullfile(clean_dir,sprintf('%s.hdf5',additionalFilterString)));
            h5File = fullfile(clean_dir,h5File.name)
            disp("\n \n \n \n \n \n ")
            disp("##########################")
            disp("#########################")
            disp(" running next Files   ")
            %%
            try
                run PlotRefinementREsults_staggered_NbX_Template.m
               catch e
                warning("%s did not work",clean_dir)
                fileID = fopen('log_IL.txt','a');
                fprintf(fileID,'The identifier was:\n%s',e.identifier);
                fprintf(fileID,'There was an error! The message was:\n%s',e.message);
                fclose(fileID);
            end
            % error('end')
        end
    end
end