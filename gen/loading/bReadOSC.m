function [ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadOSC(InputUser)
%BREADHDF5 Builds EBSD data from a HDF5 file
%%%
%% Put in a crop parameter for cropping the EBSD map, as percentage of the height
%%%
if isfield(InputUser,'breadhdf5_warnings_off')
    if InputUser.breadhdf5_warnings_off==1
        warn=0;
    end
else
    warn=1;
end

%create structures
MapData=struct;
MicroscopeData=struct;
PhaseData=struct;
EBSPData=struct;
InputUser.OSCFullFile=fullfile(InputUser.OSC_folder,InputUser.OSC_file);
disp(InputUser.OSCFullFile)
EBSPData.OSC_loc=InputUser.OSCFullFile;
if isfield(InputUser,'BinaryFile')
    EBSPData.edax_binary_loc=fullfile(InputUser.OSC_folder,InputUser.BinaryFile);
    
else
    [~,fname_osc,~]=fileparts(InputUser.OSC_file);
    try
        InputUser.BinaryFile = sprintf( '%s.up2', fname_osc );
        tempdir = fullfile(InputUser.OSC_folder, InputUser.BinaryFile);
        if isfile( tempdir );
            EBSPData.edax_binary_loc = tempdir;
        else
            error( 'Binary file does not exist' );
        end
    catch
        try
            %%
            InputUser.BinaryFile = sprintf('%s.up1',fname_osc);
            tempdir=fullfile(InputUser.OSC_folder,InputUser.BinaryFile);
            if isfile( tempdir );                
                EBSPData.edax_binary_loc=fullfile(InputUser.OSC_folder,InputUser.BinaryFile);
            else
                error( 'Binary file does not exist' );
            end
        catch
            warning('Binary file could not be found with same name as .ang / .osc file')
        end
    end
end
%%
[~,~,extension]=fileparts(InputUser.OSC_file);
if strcmp(extension,'.ang') & ~(isfield(InputUser,'CSList'));
    ebsd_map=loadEBSD_angAstro(InputUser.OSCFullFile,'convertEuler2SpatialReferenceFrame','setting 2');
else
    if ~isfield(InputUser, 'CSList')
        InputUser.CSList = InputUser.CS
    end
    ebsd_map=loadEBSD_angAstro(InputUser.OSCFullFile,InputUser.CSList,'convertEuler2SpatialReferenceFrame','setting 2');
end
%%

%%
[MapData,MicroscopeData,diff,check_alignment,ebsd_map]=setUpBeam(MicroscopeData,MapData,ebsd_map)
%%
% <<<<<<< HEAD
% try
% diff=(ebsd_map.prop.x(2)-ebsd_map.prop.x(1));
% if diff==0
%     diff=(ebsd_map.prop.y(2)-ebsd_map.prop.y(1));
% end
% catch
%         diff=(ebsd_map.x(2)-ebsd_map.x(1));
%     if diff==0
%         diff=(ebsd_map.y(2)-ebsd_map.y(1));
%     end
%     ebsd_map.prop.x=ebsd_map.x
%     ebsd_map.prop.y = ebsd_map.y
% end
% =======

pixels_shape=max(MapData.XBeam(:))*max(MapData.YBeam(:));
MapData.deleted='none';

disp('--- Before Cleaning ---');
disp(['Pixels Shape (expected total pixels): ', num2str(pixels_shape)]);
disp(['Size of ebsd_map (actual data points): ', num2str(size(ebsd_map, 1))]);

%%
% check MTEX Version, and restructure ebsd map depending on this.
V=getMtexVersion();
VNumber = regexp(V,'\-|\s|\.','split');
if str2num(VNumber{2})>=6 
ebsd_map.prop.x = ebsd_map.x;
ebsd_map.prop.y = ebsd_map.y;
end
    %%

if pixels_shape>size(ebsd_map,1)
    if check_alignment==0
        ebsd_map(ebsd_map.prop.x==max(ebsd_map.prop.x))=[];
        [MapData,MicroscopeData,diff,check_alignment,ebsd_map]=setUpBeam(MicroscopeData,MapData,ebsd_map);
        MapData.deleted='Column';
    else
        ebsd_map(ebsd_map.prop.y==max(ebsd_map.prop.y))=[];
        [MapData,MicroscopeData,diff,check_alignment,ebsd_map]=setUpBeam(MicroscopeData,MapData,ebsd_map);
        MapData.deleted='Row';
    end
end
differenceArray = ebsd_map(1:end-1).id-ebsd_map(2:end).id;
wrongIds = find(differenceArray<-1);
for ebsdPos=wrongIds
    ebsd_map(wrongIds:end).id=ebsd_map(wrongIds:end).id-1;
end
% ebsd_map = ebsd_map.gridify(:);
% ebsd_map = ebsd_map(:);
disp('--- After Cleaning ---');
disp(['Updated size of ebsd_map after cleaning: ', num2str(size(ebsd_map, 1))]);
disp(['Deleted section: ', MapData.deleted]);

% diff=(ebsd_map.prop.x(2)-ebsd_map.prop.x(1));
% if diff==0
%     diff=(ebsd_map.prop.y(2)-ebsd_map.prop.y(1));
% end
% >>>>>>> 0c153dd81ef229d7cb0e4524cbef9369053c5599
%%
EBSPData.ebsd_map=ebsd_map;
% disp(InputUser.CS)
if isfield(InputUser,'crop')
    diff=(ebsd_map.x(2)-ebsd_map.x(1));

    x_max=max(ebsd_map.x);
    y_max=max(ebsd_map.y);

   % xmin ymin dx dy
    fractions=InputUser.crop.fractions;
    steps=round(y_max/fractions/diff-diff)*diff;
    i=InputUser.crop.step_number;

    region=[0, steps*i, x_max,steps];     
    condition=inpolygon(ebsd_map,region);
    ebsd_map=ebsd_map(condition) ;
end
    
% %% this is only for testing right now
% condition=inpolygon(ebsd_map,region)
% ebsd_map=ebsd_map(condition)
%%


try
    MapData.PHI=ebsd_map.orientations.Phi;
    MapData.phi2=ebsd_map.orientations.phi2;
    MapData.phi1=ebsd_map.orientations.phi1;
catch
    if warn==1
        warning('Orientation data loading');
    end
end

try
   
    MapData.Phase=ebsd_map.phase;

catch
      warning('Phase data not loaded')
end
try
    MapData.MAD=ebsd_map.prop.fit;
   
    %     Microscope
% %     MapData.MADPhase=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name   '/EBSD/Data/MADPhase']);
%     MapData.NIndexedBands=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/NIndexedBands']);
%     MapData.RadonQuality=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name '/EBSD/Data/RadonQuality']);
%     MicroscopeData.MaxRadonBandCount=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MaxRadonBandCount']);
%     MicroscopeData.MADMax=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MADMax']);
catch
     if warn==1
    warning('Indexing of patterns not loaded')
    end
end

%%
    if isfield(ebsd_map.prop,'iq')
        MapData.RadonQuality=ebsd_map.prop.iq;
        MicroscopeData.RadonQuality=ebsd_map.prop.iq;
    elseif isfield(ebsd_map.prop,'imagequality')
        MapData.RadonQuality=ebsd_map.prop.imagequality;
        MicroscopeData.RadonQuality=ebsd_map.prop.imagequality;
    end
    %%


    if isfield(ebsd_map.prop,'ci')
        MapData.confidenceindex=ebsd_map.prop.ci;
        MicroscopeData.confidenceindex=ebsd_map.prop.ci;
    elseif isfield(ebsd_map.prop,'confidenceindex')
        MapData.confidenceindex=ebsd_map.prop.confidenceindex;
        MicroscopeData.confidenceindex=ebsd_map.prop.confidenceindex;
    end


    if isfield(ebsd_map.prop,'semsignal')
        MapData.semsignal=ebsd_map.prop.semsignal;
        MicroscopeData.SEMImage=ebsd_map.prop.semsignal;
    elseif isfield(ebsd_map.prop,'sem_signal')
        MapData.semsignal=ebsd_map.prop.sem_signal;
        MicroscopeData.SEMImage=ebsd_map.prop.sem_signal;
    end

       %% just load data from ebsd_map properties.
       ebsd_fieldnames = fieldnames(ebsd_map.prop)
       for ebsd_fNum = 1: numel(ebsd_fieldnames)
           MapData.props.(ebsd_fieldnames{ebsd_fNum}) = ebsd_map.prop.(ebsd_fieldnames{ebsd_fNum})
       end
% catch
%      if warn==1
%     warning('Indexing of patterns not loaded')
%     end
% end
%%
try
    diff=(ebsd_map.prop.x(2)-ebsd_map.prop.x(1));
    check_alignment=0;
    if diff==0
        diff=(ebsd_map.prop.y(2)-ebsd_map.prop.y(1));
        check_alignment=1;
    end
    %somehow for the .ang file this was differeent. 
    
    % else
    MapData.XBeam=uint16(ebsd_map.prop.x/diff+1);
    MapData.XSample=ebsd_map.prop.x+1; %note that 1 is added
    MapData.YBeam= uint16(ebsd_map.prop.y/diff+1);
    MapData.YSample=ebsd_map.prop.y+1; %%assume square steps
    MapData.PMap=transpose(1:numel(MapData.YSample));
    MicroscopeData.XSTEP=diff;
    MicroscopeData.YSTEP=diff;
    % end
catch
    if warn==1
    warning('Beam and sample data not loaded')
    end
end
%%
try
    %%
    % MicroscopeData.NCOLS=uint32(max(MapData.XBeam)-min(MapData.XBeam)+1);
    MicroscopeData.NCOLS=uint32(max(MapData.XBeam));%-min(MapData.XBeam)+1);
    MicroscopeData.NROWS=uint32(max(MapData.YBeam)-min(MapData.YBeam)+1);
    MicroscopeData.NPoints=size(ebsd_map.prop.x);
    if ndims(MicroscopeData.NPoints)>1
        MicroscopeData.NPoints=MicroscopeData.NPoints(1);
    end
catch
    if warn==1
    warning('Map scanning data not known')
    end
end

%%
MicroscopeData.CoordSystems='TRZM';
%%
try
    phase_info=ebsd_map.CSList;
    
    num_phases=size(phase_info{2:-1}); %first one will be unindexed
    for n=2:num_phases+1
        PhaseData(n).Name=phase_info{n}.mineral;
        PhaseData(n).LatticeConstants=[phase_info{n}.aAxis phase_info{n}.bAxis phase_info{n}.cAxis phase_info{n}.alpha phase_info{n}.beta phase_info{n}.gamma];
    end
    PhaseData(1).num_phases=num_phases;
catch
    if warn==1
    warning('Phase data not loaded')
    end
end

fns=fieldnames(MapData);
%%
if check_alignment==1;
        for i=1:size(fns,1)
            try
                temp=MapData.(fns{i});
                temp=reshape(temp,MicroscopeData.NROWS,MicroscopeData.NCOLS);
                temp=temp.';
                MapData.(fns{i})=temp(:);
            % warn if some field name cannot be rotated
            catch
                warning(sprintf("transformation for %s not possible!", fns{i} ));
            end
        end
end
% fieldnames=f
% for
% end

end

function [MapData,MicroscopeData,diff,check_alignment,ebsd_map]=setUpBeam(MicroscopeData,MapData,ebsd_map)
%         ebsd_map = ebsd_map.gridify
   %%
        try
            %%
            V=getMtexVersion();
            
            VNumber = regexp(V,'\-|\s|\.','split');
        
            if str2num(VNumber{2})<6   
                diff=(ebsd_map.prop.x(2)-ebsd_map.prop.x(1));
                check_alignment=0;
                if diff==0
                    diff=(ebsd_map.prop.y(2)-ebsd_map.prop.y(1));
                    check_alignment=1; %somehow for the .ang file this was differeent. 
                end
                MapData.XBeam=uint16(ebsd_map.prop.x/diff+1);
                MapData.XSample=ebsd_map.prop.x+1; %note that 1 is added
                MapData.YBeam= uint16(ebsd_map.prop.y/diff+1);
                MapData.YSample=ebsd_map.prop.y+1; 
            else 
                diff=(ebsd_map.x(2)-ebsd_map.x(1));
                check_alignment=0;
                if diff==0
                    diff=(ebsd_map.y(2)-ebsd_map.y(1));
                    check_alignment=1; %somehow for the .ang file this was differeent. 
                end
                MapData.XBeam=uint16(ebsd_map.x/diff+1);
                MapData.XSample=ebsd_map.x+1; %note that 1 is added
                MapData.YBeam= uint16(ebsd_map.y/diff+1);
                MapData.YSample=ebsd_map.y+1; 
                try 
                    MapData.XSampleCor = MapData.XSample+ebsd_map.prop.xoffset-1;
                    MapData.YSampleCor = MapData.YSample+ebsd_map.prop.yoffset-1;
                catch
                    warning('offset not read correctly from ebsd map. Setting same values as MapData.XSample / YSample')
                    MapData.YSampleCor = MapData.YSample-1;
                    MapData.XSampleCor = MapData.XSample-1;
                end

            end

           %%assume square steps
            MapData.PMap=transpose(1:numel(MapData.YSample));
            MicroscopeData.XSTEP=diff;
            MicroscopeData.YSTEP=diff;
        catch
        %     if warn==1
            warning('Beam and sample data not loaded')
            % end
        end
   end
%    function patternCentresGeometric()
% end
