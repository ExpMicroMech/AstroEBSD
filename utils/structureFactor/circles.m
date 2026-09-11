function c= circles(normals,theta,omega)
    for i = 1:length(normals)

      n = normals.subSet(i);

      % generate circles
      c = axis2quat(n,theta) * axis2quat(orth(n),omega)*n;

      % plot circles
%       h = [h,line(c,varargin{:},'hold','doNotDraw')]; %#ok<AGROW>

    end