function InputUser=setUpPhasesFromH5(InputUser,h5structure)
InputUser.Phase_Input = {};
countNum = 1
InputUser.Occupancy = [];
InputUser.dtSpacing = [];
pattern = 'dt([\d.]+).*?a\d+\.\d+(.*?)\d+\.\dkV';
for phaseNum = 1:numel(h5structure.Groups)
    phaseName = h5structure.Groups(phaseNum).Name;
    phaseName = phaseName(6:end); % drop the /RTM bit
    if strcmp(phaseName,'RTMPCA') || strcmp(phaseName,'tileInfo') || strcmp(phaseName,"RTMRefined")
        continue
    end
    InputUser.Phase_Input{countNum} = string(phaseName);

    
    matches = regexp(InputUser.Phase_Input{countNum} , pattern, 'tokens');

    if ~isempty(matches)
        InputUser.dtSpacing = [InputUser.dtSpacing str2num(matches{1}{1})]; % The number after dt
        occ = matches{1}{2};
        if strcmp(occ, '_')
            occ = 0.5;
        elseif strcmp(occ,'_3aNb_')
            occ = 1;
         elseif strcmp(occ,'_3aNi_')
            occ = 0;
        end
        InputUser.Occupancy = [InputUser.Occupancy occ];  % The text between a4.92 and 19.3kV
        
        % fprintf('Number after dt: %s\n', dt_number);
        % fprintf('Text between a4.92 and 19.3kV: %s\n', range_text);
    else
      disp('No matches found.');
    end
    
        countNum = countNum+1;


    % InputUser.h5PhaseInfo{countNum}.Name = InputUser.h5PhaseInfo{countNum}.cif(1:end-4)
end
%% setup phase Information
for i = 1:numel(InputUser.Phase_Input)
    InputUser.Phases = InputUser.Phase_Input;
    [~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input(i),InputUser.Phase_Folder );
     Phases = Phase_BuilderSimplified( InputUser.Phase_Input{i},InputUser.Phase_Folder );
    [cs,mineral,atoms]=loadCIF_atoms(RTM_info.cif_file);
    % cs.mineral = Phases.xtalname
    % load ebsd map
    rtmConfig = defaultRTMConfig();
    
    cs = crystalAlignEdax(cs);
    InputUser.CSList{i}=cs;
end
end