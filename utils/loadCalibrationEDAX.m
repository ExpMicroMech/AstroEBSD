function calibration = loadCalibrationEDAX(fname,varargin)
% from MTEX 6.0 (copied by Lukas Berners)
%
%
%
% read EDAX *.ang file
%
% Syntax
%
%   % change x and y values such that spatial and Euler reference frame 
%   % coincide, i.e., rotate them by 180 degree
%   ebsd = loadEBSD_ang(fname,'convertSpatial2EulerReferenceFrame')
%
%   % change the Euler angles such that spatial and Euler reference frame 
%   % coincide, i.e., rotate them by 180 degree
%   ebsd = loadEBSD_ang(fname,'convertEuler2SpatialReferenceFrame')
%
% Input
%  fname - file name
%
% Flags
%  convertSpatial2EulerReferenceFrame - 
%  convertEuler2SpatialReferenceFrame - 
%

% ebsd = EBSD;
calibration = struct();
%%
try
    %%
  % assertExtension(fname,'.ang');

  %maybe we need to introduce a nonIndexed phase
  % cs{1} = 'notIndexed';

  % read file header
  nh = 1000;
  while true
    hl = file2cell(fname,nh);

    % number of header lines
    nh = find(strncmp('#',hl,1),1,'last');
    
    if nh < length(hl) || nh>1e5, break; else, nh = nh + 1000; end
  end
  
  %phasePos = strmatch('# Phase ',hl);
  % % some ang files come with a line starting "# Phase index"
  % phasePos = find(~cellfun(@isempty, regexp(hl,['# Phase ' '\d'])));
  % if isempty(phasePos)
  %   phasePos = strmatch('# MaterialName ',hl)-1;
  % end
        
  % try
    % for i = 1:nh
        
      str = hl(1:nh);
            
      % load phase number
      calibration.steps = str2num(readByToken(str, '# TEM_PIXperUM'));
      calibration.pcx = readByToken(str,'# x-star',i);
      calibration.pcy = 1- readByToken(str,'# y-star',i);
      calibration.pcz = readByToken(str,'# z-star',i);
      calibration.WD  = str2num(readByToken(str,'# WorkingDistance'))       ;
      calibration.SampleTilt = str2num(readByToken(str,'# SampleTiltAngle'));
      calibration.CameraElevationAngle = str2num(readByToken(str, 'CameraElevationAngle'));
      calibration.CameraAzimuthAngle = str2num(readByToken(str, 'CameraAzimuthalAngle'));
      calibration.XSTEP = str2num(readByToken(str, '# XSTEP:'));
      calibration.YSTEP = str2num(readByToken(str, '# XSTEP:'));
    
% end
end
end
function value = readByToken(cellStr,token,default)
%%
  values = regexp(cellStr,[token '\s*(.*)'],'tokens');
  id = find(~cellfun(@isempty,values),1);
  if ~isempty(id)
    value = strtrim(char(values{id}{1}));
    
    if nargin > 2 && ~isempty(default) && isnumeric(default)
      value = str2num(value);
    end
    
  elseif nargin > 2
    value = default;
  else 
    value = [];
  end
  
end