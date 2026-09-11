% Lukas Berners, 10.08.2025
% witht this deck we want to automatically build the documentation for the 
% AstroEBSD decks utilising the publish command from matlab
% the rest of the documentation will be build using the autobuild tools
% from the matlab sphinx domain
%

%% Get all files from the decks directory
deckScriptFolders= dir("../decks/*.*")
deckFolders = deckScriptFolders([deckScriptFolders.isdir])
generateSourceFolders(deckFolders,'decks')
%%
% deckScriptFolders= dir("../decks/**/*.*")
% deckFolders = deckScriptFolders([deckScriptFolders.isdir])
% generateSourceFolders(deckFolders,splitString)
%%
deckScriptFiles= dir("../decks/ILS_Study/**")
generateScriptDocs(deckScriptFiles(),'decks')
%%
%% Try the same for the Documentation 
deckScriptFolders= dir("source/UserGuide/**/*.*")
deckFolders = deckScriptFolders([deckScriptFolders.isdir])
generateSourceFolders(deckFolders,'UserGuide')
%%
% deckScriptFolders= dir("../decks/**/*.*")
% deckFolders = deckScriptFolders([deckScriptFolders.isdir])
% generateSourceFolders(deckFolders,splitString)
%%
deckScriptFiles= dir("source/UserGuide/**")
generateScriptDocs(deckScriptFiles(),'UserGuide')
%%
deckScriptFolders= dir("../modules/**/*.*")
deckFolders = deckScriptFolders([deckScriptFolders.isdir])
generateSourceFolders(deckFolders,'modules')
%%
deckScriptFiles= dir("../modules/**")
generateFunctionDocs (deckScriptFiles,'modules')
%% Function Definitions utils
deckScriptFolders= dir("../utils/**/*.*")
deckFolders = deckScriptFolders([deckScriptFolders.isdir])
generateSourceFolders(deckFolders,'utils')
%%
deckScriptFiles= dir("../utils/**")
generateFunctionDocs (deckScriptFiles,'utils')


% fid = fopen('source\index.rst','a')
% S = readlines('source\index.rst')
%%
function generateFunctionDocs(deckScriptFiles,splitString)
% try
    
for fileNum = 1:numel(deckScriptFiles)
   extension = strsplit(deckScriptFiles(fileNum).name,'.')   
   %%
   extension = extension{end}
   if strcmp(extension,'m')
        addRSTforModules(deckScriptFiles(fileNum),splitString)
   end

end

end

function addRSTforModules(folder,splitString)

    %%    
    [parentFolder, tocParentFolder] = getParentFolders(folder,splitString)
    lowestLevelFolder = strsplit(parentFolder,'/')
    lowestLevelFolder = lowestLevelFolder{end}
    try
        lowestLevelFolder = sprintf("%s_TOC__",lowestLevelFolder(1:end-2))
    end
    % strsplit(deckScriptFiles(fileNum).name,'.')   
    writeName = folder.name;
    htmlName = replace(writeName,'.m','.html')
    writeName = replace(writeName,'.m','')
    %heading
    writelines( folder.name,sprintf("%s/__%s__.rst",(parentFolder),writeName));
    writelines('================',sprintf("%s/__%s__.rst",(parentFolder),writeName),WriteMode='append') ;
    
    %include raw html files %NOTE: the doc says that this is potentially
    %insecure, but for now it is the easiest way for us.
    
    writelines('.. raw:: html'...
            ,sprintf("%s/__%s__.rst",parentFolder,writeName),WriteMode='append') ;

    writelines(sprintf('\t:file: %s',htmlName)...
            ,sprintf("%s/__%s__.rst",parentFolder,writeName), WriteMode='append') ;
    % exportFolder = sprintf()
    % tryToPublish(folder,parentFolder)
    % % folderName = sprintf("%s__%s__",tocParentFolder,writeName);
    % writelines(sprintf("\t%s", folderName),sprintf('%s/%s.rst',parentFolder,lowestLevelFolder),WriteMode='append');
    %%
end




%% handle .m files
function generateScriptDocs(deckScriptFiles,splitString)
% try
    
for fileNum = 1:numel(deckScriptFiles)
   extension = strsplit(deckScriptFiles(fileNum).name,'.')   
   extension = extension{end}
   if strcmp(extension,'m')
        addRSTforScripts(deckScriptFiles(fileNum),splitString)
   end

end

end

function addRSTforScripts(folder,splitString)

    %%    
    [parentFolder, tocParentFolder] = getParentFolders(folder,splitString)
    lowestLevelFolder = strsplit(parentFolder,'/')
    lowestLevelFolder = lowestLevelFolder{end}
    try
        lowestLevelFolder = sprintf("%s_TOC__",lowestLevelFolder(1:end-2))
    end
    % strsplit(deckScriptFiles(fileNum).name,'.')   
    writeName = folder.name;
    htmlName = replace(writeName,'.m','.html')
    writeName = replace(writeName,'.m','')
    %heading
    writelines( folder.name,sprintf("%s/__%s__.rst",(parentFolder),writeName));
    writelines('================',sprintf("%s/__%s__.rst",(parentFolder),writeName),WriteMode='append') ;
    
    %include raw html files %NOTE: the doc says that thisders= dir("../decks/**/*.*")
    % deckFolders = deckScriptFolders([deckScriptFolders.isdir])
    % generateSourceFolders(deckFolders,splitString) is potentially
    %insecure, but for now it is the easiest way for us.
    
    writelines('.. raw:: html'...
            ,sprintf("%s/__%s__.rst",parentFolder,writeName),WriteMode='append') ;

    writelines(sprintf('\t:file: %s',htmlName)...
            ,sprintf("%s/__%s__.rst",parentFolder,writeName), WriteMode='append') ;
    % exportFolder = sprintf()
    tryToPublish(folder,parentFolder)
    folderName = sprintf("%s__%s__",tocParentFolder,writeName);
    writelines(sprintf("\t%s", folderName),sprintf('%s/%s.rst',parentFolder,lowestLevelFolder),WriteMode='append');
    %%
end


function tryToPublish(folder,exportFolder)
    % exportFolder = sprintf("%s/%s.rst",exportFolder,writeName)
    try
        publish(fullfile(folder.folder,folder.name),'outputDir',exportFolder )
    catch
        warning("Publishing with code running was not possible, publish without running code")
        publish(fullfile(folder.folder,folder.name),'outputDir',exportFolder,'evalCode',false)
    end
end
%% try to copy images
function copyPublishedImages(exportFolder)
    %   only for UserGuid hard coded
    splittedFolder = strsplit(exportFolder,'source') % This will of course 
            % only work if source is not in any of the previous names...
    % try
        %% assume publish will always write to .png
       images = dir(fullfile(exportFolder,"html","*.png"))
       for fileNum = 1:numel(images)
            buildDir = fullfile(splittedFolder{1},'build','html',splittedFolder{2},'html')
            copyfile(fullfile(images(fileNum).folder,images(fileNum).name),fullfile(buildDir,fileNum.Name))
       end
    % end
end

%% Handle Tocs


function makeTocEntry(folder,tocParentFolder,splitString)
            folderName = sprintf("%s%s/%s",tocParentFolder,getExportFolderName(folder),getExportFolderName(folder));
            writelines(sprintf("\t%s", folderName), sprintf('source/%s/__%s_TOC__.rst',splitString,splitString),WriteMode='append');
end



function createTocRstFiles(folder,parentFolder)
    folderName = sprintf("%s/%s",parentFolder,getExportFolderName(folder));
    % signal that they are "private" and might not be read by .git
    mkdir(folderName);
    writeName = getExportFolderName(folder);
    writelines( folder.name,sprintf("%s/%s/%s.rst",parentFolder,writeName,writeName));
    writelines('================',sprintf("%s/%s/%s.rst",parentFolder,writeName,writeName),WriteMode='append') ;
    writelines(sprintf('.. include:: %s.rst',getTocFileFromFolderName(folder))...
            ,sprintf("%s/%s/%s.rst",parentFolder,writeName,writeName),WriteMode='append') ;
    writelines(".. toctree::", sprintf("%s/%s/%s.rst",parentFolder, writeName, getTocFileFromFolderName(folder)))
end




%% General Fucntions for folder hierarchy

function [parentFolder, tocParentFolder] = getParentFolders(folder,splitString)
    parentFolder = strsplit(folder.folder,splitString);
   %%
    % disp(folder.name)
    % parentFolder = ""
    % 
    % tocParentFolder = parentFolder; % just to suprres erors might lead to unexpected results
    if numel(parentFolder)>=2
        %%
     
        
        parentFolder =parentFolder{2}
        %
        
        try
            if strcmp(parentFolder(1),'/') || strcmp(parentFolder(1),'\')
                parentFolder = parentFolder(2:end)
                parentFolder2 = strsplit(parentFolder,'\')
                parentFolder = []
                for fileNum = 1:numel(parentFolder2)
                    inputParentFolder = parentFolder2{fileNum}
                    while strcmp(inputParentFolder(1:2),'__')
                        inputParentFolder = inputParentFolder(3:end)
                    end
                    while strcmp(inputParentFolder(1:2),'__')
                        inputParentFolder = inputParentFolder(1:end-3)
                    end
                    parentFolder=sprintf('%s/__%s__',parentFolder,inputParentFolder)
                end
                parentFolder = char(parentFolder)
                %
                parentFolder = parentFolder(2:end)
                %
                % parentFolder = sprintf('__%s__',parentFolder)
                
            end
        catch

        end
        tocParentFolder = sprintf("%s",parentFolder);
        tocParentFolder = tocParentFolder(2:end)
        parentFolder = sprintf('source/%s/%s',splitString, parentFolder) ;
        % parentFolder = parentFolder(1:end-1)
        % we only need the relative path behind decks as we 
        % want to make the directory in the /src
    else
        parentFolder = sprintf('source/%s',splitString)
        tocParentFolder = ""
    end
    %%
end



function generateSourceFolders(deckFolders,splitString)
    mkdir(sprintf("source/%s",splitString))
    writelines(".. toctree::",sprintf('source/%s/__%s_TOC__.rst',splitString,splitString))
    for folderNum = 1:numel(deckFolders)
        folder = deckFolders(folderNum);
        %%
        if ~ ( strcmp(folder.name,'.') | strcmp(folder.name, '..'))
            [parentFolder,tocParentFolder] = getParentFolders(folder,splitString)
            % Build the TOC for docs.rst
            % hint: the docs.rst includes a .rst file that only lists all
            % decks as docs
            % folderName = fullfile(parentFolder,folder.name);
            makeTocEntry(folder,tocParentFolder,splitString)
            createTocRstFiles(folder,parentFolder)
            
            %
        end
    
    end
end





function exportFolderName = getExportFolderName(folder)
    exportFolderName = sprintf("__%s__",folder.name);
end

function exportFolderName = getTocFileFromFolderName(folder)
    exportFolderName = sprintf("__%s_TOC__",folder.name);
end

