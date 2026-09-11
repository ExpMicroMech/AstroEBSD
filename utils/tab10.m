function [map,num] = tab10(N)
% Qualitative colormap from MatPlotLib, for plots using the line ColorOrder.
% In MatPlotLib 2 it is named VEGA10, for MatPlotLib 3 was renamed TAB10.
%
% Copyright (c) 2017-2024 Stephen Cobeldick
%
%%% Syntax:
%  map = tab10
%  map = tab10(N)
% [map,num] = tab10(..)
%
% For MatPlotLib 2.0 improved colormaps were created for plot lines of
% categorical data. The new colormaps are introduced here:
% <http://matplotlib.org/2.0.0rc2/users/dflt_style_changes.html>
% VEGA10/TAB10 is the default Line Color Order for MatPlotLib 2 and 3.
%
% MATLAB axes ColorOrder (note that this is NOT the axes COLORMAP):
% <https://www.mathworks.com/help/matlab/creating_plots/defining-the-color-of-lines-for-plotting.html>
%
%% Examples %%
%
%%% PLOT using matrices:
% N = 10;
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*sin(x+2*n*pi/N), X(:), 1:N);
% plot(X,Y, 'linewidth',4)
% colororder(tab10(N))
%
%%% PLOT in a loop:
% N = 10;
% colororder(tab10(N))
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*sin(x+2*n*pi/N), X(:), 1:N);
% for n = 1:N
%     plot(X(:),Y(:,n), 'linewidth',4);
%     hold on
% end
%
%%% LINE using matrices:
% N = 10;
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*cos(x+2*n*pi/N), X(:), 1:N);
% line(X(:),Y)
% colororder(tab10(N))
%
%% Input and Output Arguments %%
%
%%% Inputs (**=default):
% N = NumericScalar, N>=0, an integer to define the colormap length.
%   = []**, same length as MATLAB's inbuilt colormap functions.
%
%%% Outputs:
% map = NumericMatrix, size Nx3, a colormap of RGB values between 0 and 1.
% num = NumericScalar, the length of the defining colormap (number of nodes).
%
% See Also TAB20 TAB20B TAB20C
%          CIVIDIS INFERNO MAGMA PLASMA TWILIGHT TWILIGHT_SHIFTED VIRIDIS
%          MAXDISTCOLOR BREWERMAP LINES COLORORDER COLORMAP SET AXES
if nargin<1 || isnumeric(N)&&isequal(N,[])
	N = cmDefaultN();
else
	assert(isscalar(N)&&isfinite(N)&&isreal(N),...
		'SC:tab10:N:NotRealFiniteScalarNumeric',...
		'First input <N> must be a real finite numeric scalar.')
end
%
hex = '#1f77b4#ff7f0e#2ca02c#d62728#9467bd#8c564b#e377c2#7f7f7f#bcbd22#17becf';
raw = sscanf(hex,'#%2x%2x%2x',[3,Inf]).';
num = size(raw,1);
assert(numel(hex)==num*7,'SC:tab10:hex','HEX values incorrectly parsed.')
%
map = raw(1+mod(0:N-1,num),:) / 255;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tab10
function N = cmDefaultN()
% Get the colormap size from the current figure or default colormap.
try
	F = get(groot,'CurrentFigure');
catch %#ok<CTCH> pre HG2
	N = size(get(gcf,'colormap'),1);
	return
end
if isempty(F)
	N = size(get(groot,'DefaultFigureColormap'),1);
else
	N = size(F.Colormap,1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cmDefaultN
