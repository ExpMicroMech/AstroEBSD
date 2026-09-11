function phase_name = cleanNameLatex(phase_name)
    % startIndex = regexp(phase_name,'/');
    phase_name = erase(phase_name, '/');
    % for slashes=1:size(startIndex,2)
    %     phase_name=,startIndex(slashes)-2+slashes,'');
    % end 
    startIndex = regexp(phase_name,'_');
    for slashes=1:size(startIndex,2)
        phase_name=insertAfter(phase_name,startIndex(slashes)-2+slashes,'\');
    end 

end