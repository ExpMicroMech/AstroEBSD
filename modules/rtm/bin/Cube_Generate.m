function [screen_int,facedata] = Cube_Generate(BinFile,isHex,varargin)
%MasterCube - Read a bin file from Dynamics & Build Interpolants

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

%v2 - Lukas Berners modified the code, to read in the Master patterns from
%the square Lambert projection as implemented by the Marc DeGraef group in
%EMsoft. I also just decided, to give the
% 2023_07_19 (JJJJ_MM_DD)
%
%Callahan, Patrick G.; Graef, Marc de, Microscopy and
%microanalysis 19 (2013), 1255-1265 doi:10.1017/S1431927613001840
%
normalisation='mean_std';
for i=1:size(varargin,2)
    if strcmp(varargin{i},'mean_std')
        normalisation='mean_std';
    elseif strcmp(varargin{i},'min_max')
        normalisation='min_max';
    end

end
[~,~,ext]=fileparts(BinFile);
%%
if strcmp(ext,'.bin')||strcmp(ext,'')
    fileID = fopen(BinFile, 'r', 'ieee-le');
    % if fileID == -1, error('Cannot open file: %s', filename); end
    format = 'uint';
    Data = fread(fileID, Inf, format);
    %%
    fclose(fileID);
    %find out the simulation resolution
    cube_res=sqrt((size(Data,1)-12)/6)-1;

    %build the face data for the cube
    fd=zeros(cube_res+1,cube_res+1,6);
    for n=1:6
        dstart=(n-1)*(cube_res+1)^2+1+2*n-1;
        dend=n*(cube_res+1)^2+2*n-1;
        fd(:,:,n)=flipud(rot90(reshape(Data(dstart:dend),cube_res+1,cube_res+1)));
        fd(:,:,n)=rescale(fd(:,:,n));
    end

    %normalise

    %%
    if strcmp(normalisation,'min_max')
        fd=fd-min(fd(:));
        fd=fd/max(fd(:));
    else
        fd=fd-mean(fd(:));
        fd=fd./std(fd(:));
    end
    facedata=fd;
    %sort this data
    facedata(:,:,1)=fd(:,:,3); %x+
    facedata(:,:,2)=fd(:,:,5); %y+
    facedata(:,:,3)=fd(:,:,1); %z+



    if isHex == 1 %add in a hexagonal fix
        facedata(:,:,4)=rot90(fd(:,:,1),2);
        %     facedata(:,:,5)=rot90(fd(:,:,2),2);
        facedata(:,:,5)=rot90(fd(:,:,2),3); %corrected for Shuheng/Mg
        facedata(:,:,6)=(fd(:,:,3));
    else
        facedata(:,:,4)=fd(:,:,4); %x-
        facedata(:,:,5)=fd(:,:,6); %y-
        facedata(:,:,6)=fd(:,:,2); %z-
    end

    %build the interpolants
    [gx,gy]=ndgrid(linspace(-1,1,cube_res+1),linspace(-1,1,cube_res+1));
    screen_int.p1=griddedInterpolant(gx,gy,facedata(:,:,1),'cubic'); %x+
    screen_int.p2=griddedInterpolant(gx,gy,facedata(:,:,2),'cubic'); %y+
    screen_int.p3=griddedInterpolant(gx,gy,facedata(:,:,3),'cubic'); %z+
    screen_int.p4=griddedInterpolant(gx,gy,facedata(:,:,4),'cubic'); %x-
    screen_int.p5=griddedInterpolant(gx,gy,facedata(:,:,5),'cubic'); %y-
    screen_int.p6=griddedInterpolant(gx,gy,facedata(:,:,6),'cubic'); %z-

    screen_int.isHex=isHex;
    screen_int.simSource='Bruker';
    % screen_int2.p1=griddedInterpolant(gx,gy,1+0*facedata(:,:,1),'cubic'); %x+
    % screen_int2.p2=griddedInterpolant(gx,gy,2+0*facedata(:,:,2),'cubic'); %y+
    % screen_int2.p3=griddedInterpolant(gx,gy,3+0*facedata(:,:,3),'cubic'); %z+
    % screen_int2.p4=griddedInterpolant(gx,gy,4+0*facedata(:,:,4),'cubic'); %x-
    % screen_int2.p5=griddedInterpolant(gx,gy,5+0*facedata(:,:,5),'cubic'); %y-
    % screen_int2.p6=griddedInterpolant(gx,gy,6+0*facedata(:,:,6),'cubic'); %z-
elseif strcmpi(ext,'.h5')
    filepath=BinFile;
    em_h5_n=h5read(filepath,'/EMData/EBSDmaster/mLPNH'); %% we read both southern an northern hemisphere, to be in the future able to account for monoclinic / triclinic symmetries
    em_h5_s=h5read(filepath,'/EMData/EBSDmaster/mLPSH');
    %%
    if ndims(em_h5_n)==4 %% assume, this is the same for both
        em_h5_n=sum(em_h5_n,[3,4]);
        em_h5_s=sum(em_h5_s,[3,4]);%normalisation %sometimes the lambert seems to have 4 dims, if there is more than one atom type in the unit cell
    else
        em_h5_n=sum(em_h5_n,[3]);

        em_h5_s=sum(em_h5_s,[3]);
    end

    %% normalisation
    if strcmp(normalisation,'min_max')
        em_h5_n=em_h5_n-min(em_h5_n(:));
        em_h5_s=em_h5_s-min(em_h5_s(:));
        em_h5_n=em_h5_n./max(em_h5_n(:));
        em_h5_s=em_h5_s./max(em_h5_s(:));
    else
        em_h5_n=em_h5_n-mean(em_h5_n(:));
        em_h5_s=em_h5_s-mean(em_h5_s(:));
        em_h5_n=em_h5_n./std(em_h5_n(:));
        em_h5_s=em_h5_s./std(em_h5_s(:));
    end
    %% Build the interpolant
    screen_int.p1=griddedInterpolant(em_h5_n,'cubic');
    screen_int.p2=griddedInterpolant(em_h5_s,'cubic');
    screen_int.simSource="EMsoft";
    screen_int.isHex=isHex;
    facedata=0;

elseif strcmpi(ext,'.bwkd')

    simulation_file=BinFile;
    pattern_file=[BinFile(1:end-5) '.txt'];

    [file1_array,set_proj] = read_sim(simulation_file,pattern_file);
    xwid=size(file1_array,2);
    yhig=size(file1_array,1);
    if set_proj == 0
        yhig=yhig/2;
    end

    [gy,gx]=ndgrid(linspace(-1,1,yhig),linspace(-1,1,xwid));
    screen_int=struct;
    screen_int.pp=griddedInterpolant(gy,gx,flipud(rot90(file1_array(1:yhig,:),3)),'cubic'); %rz+
    screen_int.pm=griddedInterpolant(gy,gx,flipud(rot90(file1_array(yhig+1:end,:))),'cubic'); %rz-

    screen_int.isHex=0;
    screen_int.simSource='bwkd';
    screen_int.simulation_file=simulation_file;
    screen_int.pattern_file=pattern_file;
    
    %read the simulation info in
    [sim1_data]=read_text(simulation_file);
    screen_int.kV=find_val(sim1_data,'ed_kV');
    screen_int.sim1_data=sim1_data;

elseif strcmpi(ext,'.sdf5')
    screen_int=struct;
    upper=double(h5read(BinFile,'/Data/Master/Dynamical/Upper'));
    lower=double(h5read(BinFile,'/Data/Master/Dynamical/Lower'));
    xwid=size(upper,2);
    yhig=size(upper,1);
    [gy,gx]=ndgrid(linspace(-1,1,yhig),linspace(-1,1,xwid));
    screen_int.pp=griddedInterpolant(gy,gx,upper,'cubic'); %rz+
    screen_int.pm=griddedInterpolant(gy,gx,lower,'cubic'); %rz-

    screen_int.isHex=0;
    screen_int.simSource='MapSweeper';
    screen_int.simulation_file=BinFile;

else
    disp(['The bin file ' BinFile ' does not exist']);

    cube_res=256;
    facedata=zeros(cube_res+1,cube_res+1,6);
    [gx,gy]=ndgrid(linspace(-1,1,cube_res+1),linspace(-1,1,cube_res+1));
    screen_int.p1=griddedInterpolant(gx,gy,facedata(:,:,1),'cubic'); %x+
    screen_int.p2=griddedInterpolant(gx,gy,facedata(:,:,2),'cubic'); %y+
    screen_int.p3=griddedInterpolant(gx,gy,facedata(:,:,3),'cubic'); %z+
    screen_int.p4=griddedInterpolant(gx,gy,facedata(:,:,4),'cubic'); %x-
    screen_int.p5=griddedInterpolant(gx,gy,facedata(:,:,5),'cubic'); %y-
    screen_int.p6=griddedInterpolant(gx,gy,facedata(:,:,6),'cubic'); %z-

    screen_int.isHex=0;
end




end

function val_out=find_val(string_cell,text_to_find)
    % try
    %find the element you want to read
    el_found1=arrayfun(@(x) strfind(x,text_to_find),string_cell);
    el_found2=logical(1-cellfun(@isempty,el_found1));
    ret_string=string_cell{el_found2};
    %extract value
    colon_pos=strfind(ret_string,':');
    comma_pos=strfind(ret_string,',');
    val_out=str2double(ret_string(colon_pos+1:comma_pos-1));
    % catch
        % val_out=NaN;
    % end

end

function [sim1_data]=read_text(simulation_file)
    %read the file into a cell
    sim1_fid=fopen(simulation_file);
    sim1_data=textscan(sim1_fid,'%s','Delimiter','\n');
    fclose(sim1_fid);
    sim1_data=sim1_data{1};
end