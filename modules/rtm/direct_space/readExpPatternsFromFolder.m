function Pat_exp = readExpPatternsFromFolder(InputUser,patternNumber)
    % for patternNumber = 1:numel(InputUser.PatternList)
            InputUser.PatternLoc = fullfile(InputUser.PatternFolder,InputUser.PatternList(patternNumber));
         
            % Pat_exp=double(flipud(rot90(imread(InputUser.PatternLoc))));
            Pat_exp =flipud(imread(InputUser.PatternLoc));
            if size(Pat_exp,3)>1
                Pat_exp = double(Pat_exp(:,:,1));
            end
            Pat_exp = double(Pat_exp);
    % end
end