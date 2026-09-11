function [ EBSP_out,r2] = EBSP_gen( EBSP_av,rotmat,screen_int,isHex )
%EBSP_GEN Generate a pattern for a specific orientation matrix

% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

%EBSP = ebsp coordinate system from the EBSD geometry setup
%rotmat = 3x3 rotation matrix
%screen_int = screens from a bin file reader
%ixHex = screen info from the bin file reader
%EBSP_out = pattern out


if nargin ~= 4 && isfield(screen_int,'isHex')
    isHex=screen_int.isHex;
end

lattice_scaling=[1 1 1];
for i=1:size(varargin,2)
    if strcmp(varargin{i},'lattice_scaling')
       lattice_scaling=varargin{i+1};
    end
end

if isfield(EBSP_av,'r')
    r=EBSP_av.r;
    r = [r(:,1)*lattice_scaling(1), r(:,2)*lattice_scaling(2), r(:,3)*lattice_scaling(3)];
    r= r./vecnorm(r,2,2);
    r2 =r*rotmat';
else
    r = [EBSP_av.xpts_screen(:), EBSP_av.ypts_screen(:), EBSP_av.ypts_screen(:)*0+1].*1./sqrt((EBSP_av.xpts_screen(:).^2+EBSP_av.ypts_screen(:).^2+1));
    r = [r(:,1)*lattice_scaling(1), r(:,2)*lattice_scaling(2), r(:,3)*lattice_scaling(3)];
    r= r./vecnorm(r,2,2);
    r2 = r*rotmat';
end

%sample the pattern from the interpolant
if isfield(screen_int,'pp') %use sterographic reprojection
    [i_data] = Stereo_Sample(r2(:,1),r2(:,2),r2(:,3),screen_int);
else %use cube reprojection
    [i_data] = Cube_Sample(r2(:,1),r2(:,2),r2(:,3),screen_int,isHex);
end

%reshape the output
EBSP_out=reshape(i_data,EBSP_av.size(1),EBSP_av.size(2));

end

