function buildNewPhaseFiles(rootDir,phaseFolder)
    %%
    nmlComparePat = ("ehistmin ="|"ekev ="|"ebinsize ="|"c1 ="|"c2 ="|"c3 ="|"dmin ="|"npx ="|"xtalname =");
    %% TODO: load compareKeys from a .config file

    files = dir(fullfile(rootDir,'*/**/*master*.out'));
    fNum = 1;
    %%
    % 
    % FinishedSimDirs = struct()
    % remove files, which did not run by comparint with the output form the
    % master pattern simulation
    fprintf("numel Files: %i", numel(files))
    while fNum <= numel(files)
        % fNum = 98;
        currentfile = files(fNum);
        outputs = readlines(fullfile(currentfile.folder,currentfile.name));        
        %%
        if any(contains(outputs,"Total execution time [s]"))
            files(fNum).folder = string(files(fNum).folder);
            fNum = fNum+1;

        else
           files(fNum) = [];
        end
    end
    %%
    
    fileNames = vertcat(files.folder);
    [fileNames,fileIndex] = sort(fileNames,'descend'); % sort files in descending older thus, 
    % newest elements are earlier, which are more likely to contain correct structures
    %%
    disp("--------Copying Files-------------")
    parsedFiles = [];
%%
    disp(numel(files))
    for fNum = 1:numel(files)
        %%
        if mod(fNum,5)==0
            disp(sprintf('running file %03i',fNum))
        end
        %%
        file = files(fileIndex(fNum));
        
        
        
        try
            phaseFile = dir(fullfile(file.folder,'*.pha'));
            [phaseKeys, phaseVals] = readAllPhaKeys(fullfile(phaseFile.folder,phaseFile.name));
        catch
            [phaseKeys, phaseVals] = getDefaultPhaConfig();
        end
        %
        %
        % find(strcmp(phaseKeys,"$Name"));
        % end
        % there should be only one, use try cach to be save
        try
            h5file = dir(fullfile(file.folder,'*.h5'));%.name
            h5fileName = fullfile(h5file.folder,h5file.name);
        catch
            warning(sprintf('no h5 file found in %s', file.folder))
            continue
        end
        
        %
        % nml keys
        [nmlKeys , nmlVals ] = parseNMLs(file,nmlComparePat);
        % allKeys = [phaseKeys,nmlKeys];

        xtalName = nmlVals(find(strcmp(nmlKeys,"$xtalname")));
        %%
        if numel(h5file)>=2

             h5Names = []
            for i = numel(h5file)
                h5Names = [h5Names,string(h5file(i).name)]
            end
            h5file = h5file(strfind(h5Names,xtalName));
            h5file = h5file(1); %if the xtal name is hte first one it should be fine as from sorting the correct order is correct
        end
%%

        [phaseKeys , phaseVals ] = checkMissingItems(phaseKeys,phaseVals,xtalName,h5file);
        %
        %

        for keyNum = 1:numel(nmlKeys)
            keyPosition =  strfind(phaseKeys,nmlKeys(keyNum));
            keyLines = find(~cellfun(@isempty,keyPosition));
            if ~isempty(keyLines)
                phaseKeys(keyLines) = [];
                phaseVals(keyLines) = [];
            end
        end
        %
        [nmlKeys, nmlVals] = checkMultiBin(nmlKeys,nmlVals);
%
        allKeys = [phaseKeys,nmlKeys];
        allVals = [phaseVals,nmlVals];
        % phaseKeys;
        % nmlKeys;
        % try
        Voltage = sprintf("_%.1fkV",str2num(nmlVals(find(strcmp(nmlKeys,'$Voltage')))));

        % catch
            % Voltage = sprintf("_%.1fkV",str2num(nmlVals(find(strcmp(nmlKeys,'$Voltage')))));
        % end
        dmin = sprintf("_dmin%.3f",str2num(allVals(find(strcmp(allKeys,'$dmin')))));
        c1 = str2num(allVals(find(strcmp(allKeys,'$c1'))));
        c2 = str2num(allVals(find(strcmp(allKeys,'$c2'))));
        c3 = str2num(allVals(find(strcmp(allKeys,'$c3'))));
        npx = str2num(allVals(find(strcmp(allKeys,'$npx'))));
        betheNames = "";
        if c1 ~=5
            betheNames = betheNames + sprintf("_%.1fc1",c1);
        end
        if c2 ~=10
            betheNames = betheNames + sprintf("_%.1fc2",c2);
        end
        if c3 ~=50
            betheNames = betheNames + sprintf("_%.1fc3",c3);
        end
        npxString = "";
        if npx ~=1000
            npxString = sprintf("_%inpx",npx) ;
        end
        %
        phaseFileName = sprintf("%s%s%s%s%s" ,xtalName, ...
                Voltage, dmin, npxString, betheNames);
        % phaseFileFullName = fullfile(phaseFolder,'phasefiles', sprintf('%s.pha',phaseFileName))
        
        if ~any(strcmp(parsedFiles,phaseFileName))
            fileID = fopen(fullfile(phaseFolder,'phasefiles', sprintf('%s.pha',phaseFileName)),'w');
            for keyNum = 1:numel(allKeys)
                fprintf(fileID,'%s\n',allKeys(keyNum));
                fprintf(fileID,'%s\n',allVals(keyNum));
            end
            fclose(fileID);
            %%
            try
                cifname = sprintf("%s.cif",xtalName);
                copyfile(fullfile(file.folder, cifname), fullfile(phaseFolder, "cifs", cifname))
            catch
                warning('could not copy cif file')
            end
                
            try
                if ~isfile(fullfile(phaseFolder,"dynamic_templates",h5file.name))
                    copyfile(h5fileName,fullfile(phaseFolder,"dynamic_templates",h5file.name));
                end
            catch
                warning(sprintf('Could not copy h5 file: %s',h5fileName))
            end
        end
        parsedFiles = [parsedFiles phaseFileName];
    end
end

function [keys, vals] = parseNMLs(folder,compareKeys)
%  
    keys = [];
    vals = [];
    nml_files = dir(fullfile(folder.folder,'*.nml'));
    allparNum = 1;
    for nmlNum = 1: numel(nml_files)
        nmlFile = nml_files(nmlNum);
        nmlContents = readlines(fullfile(nmlFile.folder, nmlFile.name));
    
        parPosition =  strfind(nmlContents,compareKeys);
        parLines = find(~cellfun(@isempty,parPosition));

        for parNum = 1: numel(parLines)
            nmlLine = nmlContents(parLines(parNum));
            matches = regexp(nmlContents(parLines(parNum)),["(\w+)\s*=\s*([\d.]+|[\S.]+)"],'tokens');
            if ~isempty(matches)
                key = string(matches{1}{1});  % Extracted word

                if strcmp(key,"xtalname")
                    value = string(matches{1}{2});
                    value = strrep(value,"'","");
                    value = strrep(value,".xtal","");%strip .xtal ending

                else
                    value = string(matches{1}{2});             % % Extracted number as a double
                end
                keys = [keys sprintf("$%s",key)]; % add it directly here, saves loop later
                vals = [vals value];
                allparNum = allparNum+1;
            end
        end
        %%
    end
    % 
end
% function checkKeyExists()
% end
function latticeParameters = getLatticeParametersEMsoftH5(h5file)
%%
    try
      latticeParameters=   h5read(h5file,'/CrystalData/LatticeParameters');
    catch
        latticeParameters = [0 0 0 0 0 0 ];
        warning(sprintf("could not read lattice parameters from %s \n returning Zeros", h5file))
    end
end
function isHex = checkHexH5(h5file)
    spaceGroup =  h5read(h5file,'/CrystalData/SpaceGroupNumber');
    if spaceGroup<=194 & spaceGroup>=143
        isHex = 1;
    else
        isHex = 0;
    end
end

function [phaseKeys, phaseVals] = checkMissingItems(phaseKeys,phaseVals,xtalName,h5file)

    h5fileName = fullfile(h5file.folder,h5file.name);
    latticeParameters = getLatticeParametersEMsoftH5(h5fileName);
    %%        
    NamePos = find(strcmp(phaseKeys,"$Name"));
    phaseVals(NamePos) = xtalName;
    LatticePos = find(strcmp(phaseKeys,"$Lattice"));
    phaseVals(LatticePos) = sprintf("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",...
        latticeParameters(1),latticeParameters(2),latticeParameters(3),...
        latticeParameters(4),latticeParameters(5),latticeParameters(6));


    %%
    cifPos = find(strcmp(phaseKeys,"$cif"));
    if isempty(cifPos)
        phaseVals(end+1) =sprintf("%s.cif",xtalName);
        phaseKeys(end+1) ="$cif";
        warning('guessing cif name from simulation information');
    end
    %%
    dynPos = find(strcmp(phaseKeys, "$dynamic"));
    if isempty(dynPos)
        phaseVals(end+1) = h5file.name;
        phaseKeys(end+1) = "$dynamic";
        warning('guessing h5 name from simulation information');
    end
    %%
    isHexPos = find(strcmp(phaseKeys,"$isHex"));
    if isempty(isHexPos)
        % ixHex = 
        phaseVals(end+1) = checkHexH5(h5fileName);
        phaseKeys(end+1) ="$isHex";
        warning('guessing isHex name from  h5 simulation information');
    end


end
% 
% function getPha() = 
% 
% end
function [nmlKeys, nmlValues] = checkMultiBin(nmlKeys,nmlValues)
%%    
    maxBin = str2num(nmlValues(find(strcmp(nmlKeys,'$ekev'))));
    minBin = str2num(nmlValues(find(strcmp(nmlKeys,'$ehistmin'))));
    binSize = str2num(nmlValues(find(strcmp(nmlKeys,'$ebinsize'))));
    nmlKeys(end+1) ="$Voltage";%nmlKeys
    if ((maxBin-minBin) >= binSize)
        nmlValues(end+1) =maxBin;
    else
        nmlValues(end+1) =minBin;
    end
    
    %%
end
