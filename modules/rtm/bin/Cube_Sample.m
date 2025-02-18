function [i_data,screen_number] = Cube_Sample(xjs,yjs,zjs,screen_int,isHex)
%Sample a diffraction cube for intepolation

% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHERe IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

% inputs -
% xjs, yjx, zjs = x y and z of sampling vectors
% screen_int = structure taken from MasterCube
% outputs -
% i_data = intensity
% screen_number = which screen was used from the cube
% if ~isfield(screen_int, "simSource")
%    screen_int.simSource='Bruker'

%v2 - Lukas Berners modified the code, to read in the reference patterns from
%the square Lambert projection as implemented by the Marc DeGraef group in
%EMsoft.  Callahan, Patrick G.; Graef, Marc de, Microscopy and
%microanalysis 19 (2013), 1255-1265 doi:10.1017/S1431927613001840
%
% 2023_07_19 (JJJJ_MM_DD)


if ~isfield(screen_int,'simSource')
    screen_int.simSource='Bruker'
end

if ~isfield(screen_int,'isHex')
    isHex=screen_int.isHex;
end
%%
%%
if strcmp(screen_int.simSource,'Bruker')
    magv = sqrt(xjs.^2 + yjs.^2 + zjs.^2);
    rjs=[xjs(:)./magv(:),yjs(:)./magv(:),zjs(:)./magv(:)];

    if isHex == 1
        preA=[0.5*sqrt(2)/sqrt(3),1/sqrt(2),1/sqrt(3);...
            -sqrt(2)/sqrt(3),0,1/sqrt(3);...
            0.5*sqrt(2)/sqrt(3) -1/sqrt(2) 1/sqrt(3)];
        %rotate [111] onto [001] and put Y//a2
        rjs=rjs*preA';
    end

    %find which screen they should sample
    [~,screen_number]=max(abs(rjs),[],2);
    absmax = rjs(sub2ind(size(rjs),1:size(rjs,1),screen_number'));
    signmax=sign(absmax);
    screen_number(signmax==-1)=screen_number(signmax==-1)+3;

    %sample the screens
    i_data=zeros(size(screen_number));

    %cycle through each screen
    for num_screen=1:6
        s_inds=find(screen_number==num_screen);
        rjs_screen=rjs(s_inds,:)./transpose(repmat(absmax(s_inds).*signmax(s_inds),[3,1]));

        switch num_screen %these cases have been built through trialing
            case 1 %x+
                i_data(s_inds)=screen_int.p1(rjs_screen(:,2),-rjs_screen(:,3));
            case 2 %y+
                i_data(s_inds)=screen_int.p2(-rjs_screen(:,3),rjs_screen(:,1));
            case 3 %z+
                i_data(s_inds)=screen_int.p3(rjs_screen(:,2),rjs_screen(:,1));
            case 4 %x-
                i_data(s_inds)=screen_int.p4(rjs_screen(:,2),rjs_screen(:,3));
            case 5 %y-
                i_data(s_inds)=screen_int.p5(rjs_screen(:,3),rjs_screen(:,1));
            case 6 %z-
                i_data(s_inds)=screen_int.p6(-rjs_screen(:,2),rjs_screen(:,1));
        end


    end

    if strcmp(screen_int.isHex,'1') || screen_int.isHex==1
        %
        %     disp('i')
        %     theta=30*2*pi/360
        %     rotz_30=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1] %%
        %     for now due to the symmetry elements, swapping x and y axes seems to
        %     be more straight forward in this example.
        theta=30*2*pi/360;
        rotz_30=[cos(theta) -sin(theta); sin(theta) cos(theta)] ;
        %
        XY=[xjs(:),yjs(:)];
        %  theta=30
        rotXY=XY*rotz_30;
        xjs=reshape(rotXY(:,1), size(xjs,1), []);
        yjs=reshape(rotXY(:,2), size(yjs,1), []);

        %     direction = [1 0 0];
        %     rotate(s,direction,25)
    end

elseif strcmp(screen_int.simSource,'EMsoft')

    pat_s=size(screen_int.p1.Values);
    %%
    ind_z=zjs>=0;
    [rjs_x_lamb,rjs_y_lamb]=project_sphere_to_lambert(xjs,yjs,zjs);
    %%
    rjs_x=(rjs_x_lamb/(sqrt(pi))+1)*(pat_s(1)-1)/2+1; %scaling
    rjs_y=(rjs_y_lamb/(sqrt(pi))+1)*(pat_s(2)-1)/2+1;
    % rjs_x=rescale(rjs_x,1,pat_s(1));
    % rjs_y=rescale(rjs_y,1,pat_s(2));
    % disp(max(rjs_
    sample_img=zeros(size(xjs));
    sample_img(ind_z)=screen_int.p1(rjs_x(ind_z),rjs_y(ind_z)); %interpolation
    sample_img(~ind_z)=screen_int.p2(rjs_x(~ind_z),rjs_y(~ind_z));
    %%
    % imshow(pat_s)
    %%
    i_data=sample_img;
    screen_number=2;
    % i_data=rescale(i_data,0,1);
    % ind_z=zjs>=0;
    % [rjs_x,rjs_y]=project_sphere_to_lambert(xjs(~ind_z),yjs(~ind_z),zjs(~ind_z));
    % rjs_x=rjs_x/(sqrt(pi))*(pat_s(1))/2+pat_s(1)/2;
    % rjs_y=rjs_y/(sqrt(pi))*(pat_s(2))/2+pat_s(2)/2;
    % sample_img(~ind_z)=screen_int.p2(rjs_x,rjs_y);
    %% for controlling the correct transformation, this may be activated
    %  figure
    % plot(rjs_x,rjs_y,'k.')
    % figure
    % h=pcolor(rjs(:,:,1),rjs(:,:,2),lamb_inv_n(:,:,4))
    % h.EdgeColor='none'
    % %%%%
elseif strcmp(screen_int.simSource,'bwkd')
    [i_data] = Stereo_Sample(xjs,yjs,zjs,screen_int);
    screen_number=3;
elseif strcmp(screen_int.simSource,'MapSweeper')
    [i_data] = Stereo_Sample(xjs,yjs,zjs,screen_int);
    screen_number=4;

end

end



function [ebsp_i] = Stereo_Sample(xjs,yjs,zjs,screen_int_s)
%STEREO_SAMPLE Sample from the stereographic projection

r_s=[xjs(:),yjs(:),zjs(:)];
r_s_norm=r_s./repmat(sqrt(dot(r_s,r_s,2)),1,3);

rn_x=r_s_norm(:,1);
rn_y=r_s_norm(:,2);
rn_z=r_s_norm(:,3);

s_inds_p=find(rn_z >= 0);
s_inds_m=find(rn_z  < 0);

pn_xp= rn_x(s_inds_p)./(1+rn_z(s_inds_p));
pn_yp= rn_y(s_inds_p)./(1+rn_z(s_inds_p));

pn_xm=-rn_x(s_inds_m)./(1-rn_z(s_inds_m));
pn_ym=-rn_y(s_inds_m)./(1-rn_z(s_inds_m));

ebsp_p=screen_int_s.pp(pn_xp,pn_yp);
ebsp_m=screen_int_s.pm(pn_xm,pn_ym);

ebsp_i=zeros(size(rn_x));
ebsp_i(s_inds_p)=ebsp_p;
ebsp_i(s_inds_m)=ebsp_m;


end

function [rjs_x,rjs_y]=project_sphere_to_lambert(xjs,yjs,zjs)
%% This function samples a vector from a sphere onto the lambert square projection, see publication from Calahan, deGraef 2013

% xjs=lamb_inv_n(:,:,1);
% yjs=lamb_inv_n(:,:,2);
% zjs=lamb_inv_n(:,:,3);
xypos=abs(xjs)<=abs(yjs);

%%
xc=xjs; %c meaning for computation
yc=yjs;
zc=zjs;
xc(xypos)=0; %this is necessary to not lose elements during computation
yc(xypos)=0;
zc(xypos)=0;
%% Coordinate transform
prefac=sign(xc).*sqrt(4*(1-abs(zc))); %modification to 4 instead of 2 to use the abs(zc) (not only northern hemisphere here)
rjs_x_plus=prefac*sqrt(pi)/2;
rjs_y_plus=prefac.*(2/sqrt(pi).*atan(yc./xc));
% rjs_y_plus(xypos)=0;
rjs_y_plus(isnan(rjs_y_plus))=0;
%%
% figure
% plot(rjs_x_plus(nan),rjs_y_plus(nan));
%%
% figure
% plot(xc,yc,'k.')
% hold on


%%
% xypos=abs(xjs)>abs(yjs);
xc=xjs; %c meaning for computation
yc=yjs;
zc=zjs;
xc(~xypos)=0; %this is necessary to not lose elements during computation
yc(~xypos)=0;
zc(~xypos)=0;
% plot(xc,yc,'g.')
%%
prefac=sign(yc).*sqrt(4*(1-abs(zc)));
rjs_x_minus=prefac.*(2/sqrt(pi).*atan(xc./yc));
rjs_y_minus=prefac*sqrt(pi)/2;
rjs_x_minus(isnan(rjs_x_minus))=0;
%%
rjs_y=rjs_y_plus + rjs_y_minus;
rjs_x=rjs_x_plus + rjs_x_minus;
end