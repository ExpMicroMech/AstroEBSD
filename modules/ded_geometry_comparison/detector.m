%% detector.m
% By Lukas Berners (2024)
classdef detector
  
  properties
    ncols  % width of the detector in pixel
    nrows  % height of the detector in pixel
    dist   % distance detector to scanning surface
    patterCenter  % (x,y) coordinates of the pattern center
    bounds % bounds of the detector in gnonomic projections [xmin, xmax, ymin, ymax]
    x      % coordinates of the pixels in the gnonomic projection
    y      %
    radius=0
  end
  
  properties (Dependent = true)    
    nodesS2 % the pixel positions as points on the sphere
    edges   %
    vertices
  end
  
  properties (Constant)
    proj = gnonomicProjection;
  end
  
  methods
    function det = detector(ncols, nrows, dist, pc)
    %%
      det.ncols = double(ncols);
      det.nrows = double(nrows);
      det.dist = double(dist);
      det.patterCenter = double(pc);
      
      det.bounds = double([([0,1] - pc(1))*det.ncols/det.nrows,...
        (pc(2) - [1,0])] ./ det.dist) ;
     
      % generate x and y values of the detector positions in the gnonomic projection
      [det.x,det.y] = ndgrid(linspace(det.bounds(1),det.bounds(2),det.ncols),...
        linspace(det.bounds(3),det.bounds(4),det.nrows));
      %%
    end
    
    function v = get.nodesS2(det)

      % project x,y values onto the sphere
      v = det.proj.iproject(det.x,det.y);
      
    end
        
    
    function e = get.edges(det)
      
      n = det.nodesS2;
      e = [n(1,1),n(1,end),n(end,end),n(end,1)];
            
    end
    
    function v = get.vertices(det)
      
      n = det.nodesS2;
      v = normalize([cross(n(1,1),n(1,end)),cross(n(1,end),n(end,end)),...
        cross(n(end,end),n(end,1)),cross(n(end,1),n(1,1))]);
            
    end
    
    function [x,y] = vec2xy(det,v)
      
      [x,y] = det.proj.project(v);
      
    end
    
    %%
    function mask = S2CutOffMask(det,delta)
      %%
      if nargin >= 1, delta = 0.1; end
      mask = S2FunHandle(@(v) evalMask(v));
      
      function value = evalMask(v)
        
        [x,y] = detector.proj.project(v);
        
%         if det.radius>=0
%         cs=[size(x)]
%         r_thresh=det.radius*cs(1)/2;
% 
%         [xgrid,ygrid]=meshgrid(1:cs(2),1:cs(1));
%         r_grid=sqrt((xgrid-cs(1)/2).^2+(ygrid-cs(2)/2).^2);        
%         temp=r_grid>=r_thresh
%         %%
%         x(r_grid>=r_thresh) = 0;
%         y(r_grid>=r_thresh) = 0;
%         end
        %cutoff = @(t) (1+erf(t ./ delta-2))./2;
        cutoff = @(t) max(0,min(1,t/delta));

        value = cutoff(x - det.bounds(1)) .* cutoff(y-det.bounds(3)) .*  ...
          cutoff(det.bounds(2)-x) .* cutoff(det.bounds(4)-y) .* (v.z > 0);
        %%
      end
%%
    end
    
    function pattern = simulatePattern(det,master,ori,flux,bg)
      % simulate a Kikushi pattern given a master pattern
      
      if nargin == 2, ori = rotation.id; end
      pattern = reshape(master.eval(ori \ det.nodesS2),det.nrows,det.ncols);
      
      % add some noise
      if nargin < 4, flux = 0; end
      if nargin < 5, bg = 0; end
      if nargin > 2 && flux>0, pattern = randp(flux*(pattern+bg))-flux*bg; end

      pattern = reshape(pattern,size(det.x)) ;
    end
    
    function [pHarm,plan] = pattern2Fun(det,pattern,varargin)

      plan = getClass(varargin,'struct',struct('mask',[],'S2G',[],'W',[]));
      if isempty(plan.mask)
        plan.mask = det.S2CutOffMask(get_option(varargin,'delta',0.3)); 
      end
      
      if check_option(varargin,'quadrature')
        if check_option(varargin,'ori')
            ori=get_option(varargin,'ori')
            %%
            if isempty(plan.S2G)
              [plan.S2G, plan.W] = quadratureS2Grid(2*get_option(varargin,'bandwidth',128));
              plan.maskDiscrete = plan.mask.eval(plan.S2G);
              plan.isInside = plan.maskDiscrete > 0.1;
              plan.maskDiscrete = plan.maskDiscrete(plan.isInside);
            end
%%
            plan.S2G_rot=ori*plan.S2G
            plan.S2G_rot=plan.S2G_rot'
            %%
        % detector positions of the quadrature grid
            [xQ,yQ] = det.vec2xy(plan.S2G(plan.isInside));
            v = zeros(size(plan.S2G));

            v(plan.isInside) = plan.maskDiscrete .* interp2(det.x.',det.y.',pattern.',xQ,yQ);
    %%
            pHarm = S2FunHarmonic.quadrature(plan.S2G_rot,v,'weights',plan.W,varargin{:});
%         end
        else
            if isempty(plan.S2G)
              [plan.S2G, plan.W] = quadratureS2Grid(2*get_option(varargin,'bandwidth',128));
              plan.maskDiscrete = plan.mask.eval(plan.S2G);
              plan.isInside = plan.maskDiscrete > 0.1;
              plan.maskDiscrete = plan.maskDiscrete(plan.isInside);
            end
            
        % detector positions of the quadrature grid
            [xQ,yQ] = det.vec2xy(plan.S2G(plan.isInside));
            v = zeros(size(plan.S2G));
    
            v(plan.isInside) = plan.maskDiscrete .* interp2(det.x.',det.y.',pattern.',xQ,yQ);
    
            pHarm = S2FunHarmonic.quadrature(plan.S2G,v,'weights',plan.W,varargin{:});
        end
      else

        maskDiscrete = plan.mask.eval(det.nodesS2);
        pHarm = S2FunHarmonic.approximation(det.nodesS2(:),maskDiscrete(:) .* pattern(:),'bandwidth',128,varargin{:});
       
      end
      %%
    end

    function [plan,v] = pattern2Fun_multi_pat(det,pattern,varargin)

      plan = getClass(varargin,'struct',struct('mask',[],'S2G',[],'W',[]));
      if isempty(plan.mask)
        plan.mask = det.S2CutOffMask(get_option(varargin,'delta',0.05)); 
      end
%         ori=orientation.by_euler
%       if check_option(varargin,'quadrature')
%         if check_option(varargin,'ori')
            ori=get_option(varargin,'ori');
            %%
            if isempty(plan.S2G)
%%
                [plan.S2G,plan.W] = quadratureS2Grid(2*get_option(varargin,'bandwidth',128));
%                 plan.S2G_rot=plan.S2G
                plan.S2G_rot=inv(ori)*plan.S2G;
              plan.S2G_rot=plan.S2G_rot';
%%            %
              plan.maskDiscrete = plan.mask.eval(plan.S2G_rot);
              plan.isInside = plan.maskDiscrete > 0.1;
              plan.maskDiscrete = plan.maskDiscrete(plan.isInside);
            end
%%
%             plan.S2G_rot=inv(ori)*plan.S2G
%             plan.S2G_rot=plan.S2G_rot'
            %%
            
        % detector positions of the quadrature grid
            [xQ,yQ] = det.vec2xy(plan.S2G_rot(plan.isInside));
            v = zeros(size(plan.S2G_rot));
            plan.v_dist=zeros(size(plan.S2G_rot));
            plan.v_dist_z=zeros(size(plan.S2G_rot));
%             plan.mult=zeros(size(plan.S2G_rot));
            % plan.mult=plan.mask.eval(plan.S2G_rot);
            %%
            plan.mult=plan.isInside*1; % makes it double
            plan.mult(plan.isInside)=plan.maskDiscrete.*plan.mult(plan.isInside);
            v(plan.isInside) =  plan.maskDiscrete.*interp2(det.x.',det.y.',pattern.',xQ,yQ);%
            
            plan.v_dist(plan.isInside)=sqrt((xQ).^2 + (yQ).^2);
            max_d=max(plan.v_dist(plan.isInside));
            plan.v_dist(plan.isInside)=plan.maskDiscrete.*(max_d-plan.v_dist(plan.isInside));
            plan.v_dist(~plan.isInside)=0;
%%
            plan.v_dist_z(plan.isInside)=sqrt((xQ).^2 + (yQ).^2+1);
            max_d_z=max(plan.v_dist_z(plan.isInside));
            plan.v_dist_z(plan.isInside)=plan.maskDiscrete.*(max_d_z-plan.v_dist_z(plan.isInside));
            plan.v_dist_z(~plan.isInside)=0;
            plan.flat=plan.isInside*1;
            plan.flat(plan.isInside)=plan.maskDiscrete.*plan.flat(plan.isInside);
%             plan.
            %%
            % plan.v_dist_z(plan.isInside)=max(plan.v_dist_z(plan.isInside))
            % plan.v_dist_z(~plan.isInside)=0;
    %%
%             pHarm = S2FunHarmonic.quadrature(plan.S2G,v,'weights',plan.W,varargin{:});
%         end
%         else
% %%
% figure
% plot(plan.S2G,v)
           
      %%
%     end
%       end
    end

  
  end
end