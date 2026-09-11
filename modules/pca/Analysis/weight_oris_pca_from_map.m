function extractions=weight_oris_pca(tile,MapData,extractions)
% tile should contain the information about the x and y coordinates of the
% respective patterns. 
%MapData is supposed to contain the orientation data and the respective
%phases.
%
%In the future we would need to account for multi-phase alloys and then
%check for each component which phase is prevalent. In principle this would
%also need checking if this is correct and transformation between e.g.
%cubic and hexoganl phases. which should on the other hand allready nicely
%be covered by the hough. For other cases (eg. the orthorombic we could
%integrate that bit in the psudosymmetry algorithm
%
%


% MapData=MapInfo.MapData

names = fieldnames(tile);
match = {};
for i = 1:length(names)
  if ~isempty(regexp(names{i},'Tile\d*'))
      match{end+1} = names{i};
  end
end
t_size=size(match,2);
%% need for later but we do not need to loop over it
edges=[0:0.2:2*pi];
property=["phi1","PHI","phi2"];
if isfield(MapData,'g_num')
    
    tile_data=get_data(MapData,property);
    eulers=[];
    eulers(:,1)=tile_data.phi1;
    eulers(:,2)=tile_data.PHI;
    eulers(:,3)=tile_data.phi2;
    reshaping=0;
    rp=1;
    while reshaping==0
        try
        eulers=reshape(eulers,tile.dims(1,rp),tile.dims(2,rp),3);
        eulers=pagetranspose(eulers);
    %%
    
        eulers=reshape(eulers,double(tile.dims(1,rp))*double(tile.dims(2,rp)),3);
        reshaping=1;
        catch
            rp=rp+1;
        end
    end
    % eulers=reshape(eulers,tile.dims(1,1),tile.dims(2,1),3);
end
%% check which way round we have the datasdfasdf
if MapData.XBeam(1)-MapData.XBeam(2)==0
    stacking=1;
else
    stacking=2;
end

disp('start')
%%
% tile_name='Tile1';
% p_num=1;

 names = fieldnames(extractions);
    match = {};
    for i = 1:length(names)
        if ~isempty(regexp(names{i},'Tile'))
            match{end+1} = names{i};
        end
    end
%% Loop over all input patterns and extract the grid positions and intensities
% figure
for tnum=1:size(match,2)
% for tnum=1:t_size
    %%
    % tnum=2
    % tnum=2
    % t_name=sprintf('Tile%i',tnum);
    t_name=match{tnum};
%%
    if ~isfield(MapData,'g_num')
        tile_data=get_data_tiled(tile.colstart(tnum),tile.colfin(tnum),tile.rowstart(tnum),tile.rowfin(tnum),MapData,property);
        comps=reshape(tile.(t_name).filtered_map,tile.dims(1,tnum)*tile.dims(2,tnum),1);
         eulers=[];
        eulers(:,1)=tile_data.phi1;
        eulers(:,2)=tile_data.PHI;
        eulers(:,3)=tile_data.phi2;
    else
        comps=tile.(t_name).PCA_VM_num(:);
        eulers_2=eulers(:,:,1);

    end
    %
   % comps=comps(:)
    %
    % figure
    % eulers_2(comps~=1)=0
    % eulers_2(comps==1)=1
    % imagesc(eulers_2)
    %
    % figure
    % imagesc(eulers)
    %
     % figure
    % imagesc(comps)

    if stacking==1 &&  ~isfield(MapData,'g_num') 
    % % %     eulers=reshape(eulers)
        eulers=reshape(eulers,tile.dims(1,tnum),tile.dims(2,tnum),3);
        %
        eulers=pagetranspose(eulers);
        %
        eulers=reshape(eulers,double(tile.dims(1,tnum))*double(tile.dims(2,tnum)),3);
    end
    for c_num=1:tile.(t_name).NumPCA
        %
        % c_num=3
        % eulers=eulers(:);
        % clear eulers_sel;
        eulers_sel=[];
        eulers_sel=eulers(comps==c_num,:);
        % temp(:,2)=eulers(comps==c_num);
        % temp(:,3)=eulers(comps==c_num);
        % =temp;
        %%
        %
        % for now just check, wether the iqr is extremely high and hen raise a warning
        if max(iqr(eulers_sel)>0.05)
            warning('the scatter of your orientation data is large. It might be wrong')
        end
        for e_num=1:3
            N=histcounts(eulers_sel(:,e_num),edges);
            % disp([tnum,c_num,e_num])
            %
            [minind,maxind]=return_maxbins(N);
            %%
            % mean(eulers_sel(eulers_sel(:,e_num)>edges(minind) & eulers_sel(:,e_num)<edges(maxind))) %select only the eulers that are in the bin with the max count
            t.(t_name).oris(c_num,e_num)=mean(eulers_sel(eulers_sel(:,e_num)>edges(minind) & eulers_sel(:,e_num)<edges(maxind),e_num));
            
        end
    end
end
%%
% figure
% plot(eulers_sel)
% figure
% imagesc(reshape(eulers(:,1),37,37))
% % % axis xy;
% axis equal
% figure
% % imagesc(reshape(comps,37,37))
% imagesc(tile.(t_name).filtered_map)
% % axis equal;
%%

function MapData_sel=get_data_tiled(xmin,xmax,ymin,ymax,MapData,property)
    MapData_sel=struct();
    ind=MapData.XBeam>=xmin & MapData.XBeam<=xmax & MapData.YBeam>=ymin & MapData.YBeam<=ymax;
    if isStringScalar(property)
        try
            MapData_sel.(property)=MapData.(property)(ind);
        catch
            warning(sprintf('Property: %s not in input data',property))
        end
    else
        for i=1:numel(property)
            try
                MapData_sel.(property(i))=MapData.(property(i))(ind);
            catch
                warning(sprintf('Property: %s not in input data\n',property(i)))
            end
        end
    end
end
    fifunction MapData_sel=get_data(MapData,property)
    MapData_sel=struct();
    % ind=MapData.XBeam>=xmin & MapData.XBeam<=xmax & MapData.YBeam>=ymin & MapData.YBeam<=ymax;
    if isStringScalar(property)
        try
            MapData_sel.(property)=MapData.(property)(ind);
        catch
            warning(sprintf('Property: %s not in input data',property))
        end
    else
        for i=1:numel(property)
            try
                MapData_sel.(property(i))=MapData.(property(i));
            catch
                warning(sprintf('Property: %s not in input data\n',property(i)))
            end
        end
    end
end
    


end