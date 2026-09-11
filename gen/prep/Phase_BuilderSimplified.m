function Phases =  Phase_BuilderSimplified(phase_name,phase_folder)

        [phaseKeys, phaseVals] = readAllPhaKeys(fullfile(phase_folder,'phasefiles',sprintf("%s.pha",phase_name)));
        try
            for keyNum = 1: numel(phaseKeys)
                fieldName = strrep(phaseKeys(keyNum),'$',"");
                Phases.(fieldName) = phaseVals(keyNum);
            end
        catch
            warning('invalid field in Phase: %s,\n field: %s \n skiping phase',...
                phaseFiles(phaseNum).name,fieldName)
            
        end
end