function loadGitWindows()
try
    newPath = 'C:\Program Files\Git\bin'; % Adjust this path as necessary
catch
    warning('git was not loaded')
end
% Get the current PATH
currentPath = getenv('PATH');

% Add the new path to the existing PATH
updatedPath = [currentPath ';' newPath];

% Set the updated PATH for this session
setenv('PATH', updatedPath);

% Verify by checking if Git is accessible now
[status, gitVersion] = system('git --version');
disp(gitVersion);
% git version 2.37.3.windows.1
end