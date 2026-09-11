function makeNewH5(filePath,ASTROSKRIPTS_LOC)
    assumeParallelInstallation=0
    if nargin==1
        assumeParallelInstallation =1
    end
    try
         if ~isunix
             loadGitWindows
        end
    catch
        warning( "could not determine operating system")
    end
% ;
    try
        disp(ASTRO_LOC())
        [~, AstroGitHash] = system(sprintf('git -C "%s" rev-parse HEAD',ASTRO_LOC()));
        % [~, AstroIsClean] = system(sprintf('git -C "%s" rev-parse
        % HEAD',ASTRO_LOC())); % Add this later
        [~, AstroOrigin] = system(sprintf('git -C "%s" config --get remote.origin.url', ASTRO_LOC()));
        [~, AstroBranch] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', ASTRO_LOC()));
    catch
        warning("could not detect Astro EBSD git repository")
        AstroGitHash = ""
        AstroOrigin = ""
        AstroBranch = ""
        
    end
    try
        % disp(ASTROSKRIPTS_LOC())
        if assumeParallelInstallation==1
            [ASTROSKRIPTS_LOC,~,~] = fileparts(mfilename("fullpath"))
            ASTROSKRIPTS_LOC  = erase(ASTROSKRIPTS_LOC, "\AstroEBSD_v2\utils") %windows
            ASTROSKRIPTS_LOC  = erase(ASTROSKRIPTS_LOC, "/AstroEBSD_v2/utils") %linux
            ASTROSKRIPTS_LOC = fullfile(ASTROSKRIPTS_LOC,'astro_skripts')
            warning('You did not specify the astroskripts location, we assume, that you install astro_skripts in the same parent directory as AstroEBSD under the name astro_skripts')
        end
        % ASTROSKRIPTS_LOC = "../" % try to see whether relative path is working.
        [~,AstroSkriptsGitHash] = system(sprintf('git -C "%s" rev-parse HEAD',ASTROSKRIPTS_LOC()));
        [~, AstroSkriptsOrigin] = system(sprintf('git -C "%s" config --get remote.origin.url',ASTROSKRIPTS_LOC()));
        [~, AstroSkriptsBranch] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', ASTROSKRIPTS_LOC()));

    catch
        warning("could not detect Astro skripts repository")
        AstroSkriptsGitHash = ""
        AstroSkriptsOrigin = ""
        AstroSkriptsBranch =""
    end   
    % try % there should be no try / catch needed here as writeDataToH5 will do the error handling
        writeDataToH5(filePath,"/gitVariables/AstroSkriptsGitHash",...
                string(AstroSkriptsGitHash),'string'); %
        writeDataToH5(filePath,"/gitVariables/AstroSkriptsOrigin",...
                string(AstroSkriptsOrigin),'string'); %
        writeDataToH5(filePath,"/gitVariables/AstroSkriptsBranch",...
                string(AstroSkriptsBranch),'string'); %
        
                writeDataToH5(filePath,"/gitVariables/AstroBranch",...
                string(AstroBranch),'string'); %
        writeDataToH5(filePath,"/gitVariables/AstroGitHash",...
                string(AstroGitHash),'string'); %
        writeDataToH5(filePath,"/gitVariables/AstroOrigin",...
                string(AstroOrigin),'string'); %

    % catch
    %     warning('Datasets allready existing')
    % end
    
    % try
    % h5write(filePath,"/gitVariables/AstroSkriptsOrigin" , string(AstroSkriptsOrigin));
    % catch
    %     warning("could not create Astro Skripts Repo variables")
    % end    
   
    % try
    % h5write(filePath,"/gitVariables/AstroSkriptsBranch" , string(AstroSkriptsBranch));
    % catch
    %     warning("could not create Astro skripts Branch name variables")
    % end
    
    % try
    % h5write(filePath,"/gitVariables/AstroSkriptsGitHash" , string(AstroSkriptsGitHash));
    %   catch
    %     warning("could not create astro skripts hash variables")
    % end
    
    % try
    % h5write(filePath,"/gitVariables/AstroGitHash" , string(AstroGitHash));
    % catch
    %     warning("could not create astro  hash variables")
    % end

    % try
    % h5write(filePath,"/gitVariables/AstroOrigin" , string(AstroOrigin));
    % catch
    %     warning("could not create astro repo variables")
    % end   
    
    % try
    % h5write(filePath,"/gitVariables/AstroBranch" , string(AstroBranch));
    %   catch
    %     warning("could not create astro branch variables")
    % end 
end