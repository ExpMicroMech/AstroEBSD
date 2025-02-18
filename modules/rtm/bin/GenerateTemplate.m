% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [GeneratedTemplatePat]=GenerateTemplate(phasenum,euler,PC_pat,screensize,MicroscopeTilt,InputUser,Phase_Folder)

%run('I:\TomMcA\GitHub\RTM_indexing\start_RTM')
%run('I:\TomMcA\GitHub\AstroEBSD\start_AstroEBSD')

%% Get phase number
%phasenum=1;
%component=1;

InputUser.Phase_Input=InputUser.Phases(phasenum);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,Phase_Folder);
[screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%generate the pattern for this phase
%xi=1;yi=1;
GMat_test=conv_EA_to_G([euler(1),euler(2),euler(3)]);

%Define all rotation matrices needed in the code
Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

Detector_tilt = Rx(MicroscopeTilt);
rottoplot=GMat_test*Detector_tilt;

PatternInfo.ScreenWidth=screensize;
PatternInfo.ScreenHeight=screensize;

[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);

[ GeneratedTemplatePat ] = EBSP_gen( EBSP_pat,rottoplot,screen_int,RTM_info.isHex ); %generate the EBSP for this iteration
end


