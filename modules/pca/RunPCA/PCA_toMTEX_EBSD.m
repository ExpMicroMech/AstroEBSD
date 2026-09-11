function ebsd_map = PCA_toMTEX_EBSD(tile,extractions,CS_List,MicroscopeData,varargin)

names = fieldnames(tile);
match = {};
for i = 1:length(names)
  if ~isempty(regexp(names{i},'Tile\d*'))
      match{end+1} = names{i};
  end
end
t_size=size(match,2);

%%
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
if isfield(MicroscopeData,'SampleTilt')
    Detector_tilt = RTM.Rx(MicroscopeData.SampleTilt);
else
    Detector_tilt = RTM.RX(70*degree);
    warning(' We set your detector tilt to 70 °, as this is a default setting for EBSD analysis');

end
if isfield(MicroscopeData,'CameraEuler');
    detector_Euler=MicroscopeData.CameraEuler;
else
    detector_Euler=[ 0 0 0];
end
g_camera=RTM.Rz(detector_Euler(3))*RTM.Rx(detector_Euler(2))*RTM.Rz(detector_Euler(1));
pseudo_rot = RTM.Rx(0);
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'pseudo_rot')
        pseudo_rot=varargin{i+1};
end







%%
eulers = zeros(size(tile.map_reshaped, 1), size(tile.map_reshaped, 2), 3);
phs = zeros( size(tile.map_reshaped) );
lpr_x = zeros( size(tile.map_reshaped) );
lpr_y = zeros( size(tile.map_reshaped) );
lpr_z = zeros( size(tile.map_reshaped) );
phases = ones( size( (tile.map_reshaped) ) );
%
for tnum = 1:size(match,2)
    %
    tile_name = match{tnum}
    num_comp = size( extractions.(tile_name).patterns, 3);
    tile_cur = tile.(tile_name);
    zeros(tile_cur.height, tile_cur.width);

    mini_map_euler1 = zeros(tile_cur.height, tile_cur.width);
    mini_map_euler2 = zeros(tile_cur.height, tile_cur.width);
    mini_map_euler3 = zeros(tile_cur.height, tile_cur.width);

    mini_map_ph = zeros(tile_cur.height, tile_cur.width);
    
    mini_map_lpr_x = zeros(tile_cur.height, tile_cur.width);
    mini_map_lpr_y = zeros(tile_cur.height, tile_cur.width);
    mini_map_lpr_z = zeros(tile_cur.height, tile_cur.width);
    
    mini_maps = [];
    
    for p_cur = 1:num_comp
        input_ori = extractions.(tile_name).rotation_data_refined(:,:,p_cur);
        input_ori = pseudo_rot * input_ori * inv(Detector_tilt) * inv(g_camera');
        if isfield( extractions.(tile_name), 'LPR')
            LPR = extractions.(tile_name).LPR(:,p_cur);
        else
            LPR = [1 1 1];
        end

        index_rep = repmat( tile_cur.PCA_VM_num == p_cur, 1, 1 );
        ph = extractions.(tile_name).peak_hight(p_cur);
        input_ori=conv_G_to_EA(input_ori);

        mini_map_euler1(index_rep) = input_ori(1);
        mini_map_euler2(index_rep) = input_ori(2);
        mini_map_euler3(index_rep) = input_ori(3);

        mini_map_ph(index_rep) = ph;
        mini_map_lpr_x(index_rep) = LPR(1);
        mini_map_lpr_y(index_rep) = LPR(2);
        mini_map_lpr_z(index_rep) = LPR(3);
    end
    eulers(tile.rowstart(tnum) : tile.rowfin(tnum), tile.colstart(tnum) : tile.colfin(tnum), 1) = mini_map_euler1;
    eulers(tile.rowstart(tnum) : tile.rowfin(tnum), tile.colstart(tnum) : tile.colfin(tnum), 2) = mini_map_euler2;
    eulers(tile.rowstart(tnum) : tile.rowfin(tnum), tile.colstart(tnum) : tile.colfin(tnum), 3) = mini_map_euler3;

    phs(tile.rowstart(tnum) : tile.rowfin(tnum), tile.colstart(tnum) : tile.colfin(tnum) ) = mini_map_ph;

    lpr_x(tile.rowstart(tnum) : tile.rowfin(tnum) , tile.colstart(tnum) : tile.colfin(tnum) ) = mini_map_lpr_x;
    lpr_y(tile.rowstart(tnum) : tile.rowfin(tnum) , tile.colstart(tnum) : tile.colfin(tnum) ) = mini_map_lpr_y;
    lpr_z(tile.rowstart(tnum) : tile.rowfin(tnum) , tile.colstart(tnum) : tile.colfin(tnum) ) = mini_map_lpr_z;
%     eulers(tile.colstart(tnum):tile.colfin(tnum),tile.rowstart(tnum):tile.rowfin(tnum),2)=mini_map_euler2;
%     eulers(tile.colstart(tnum):tile.colfin(tnum),tile.rowstart(tnum):tile.rowfin(tnum),3)=mini_map_euler3;
end
%%
% CS_List = {'notIndexed',cs};
eulers = reshape(eulers(:),[],3);
% eulers=zeros(size(tile.map_reshaped,1),size(tile.map_reshaped,2),3);
% phs=zeros(size(tile.map_reshaped))
% figure
% scatter(tile.xloaded,tile.yloaded)
% [X,Y]=ndgrid(0:0.5/2-0.5,0:0.5:489/2-0.5);
[X,Y] = meshgrid(0:1:size( tile.map_reshaped, 1)-1, 0:1:size( tile.map_reshaped, 2)-1 );

prop.x = double(X(:)*MicroscopeData.XSTEP);
prop.y = double(Y(:)*MicroscopeData.YSTEP);%%
prop.ph = phs(:);

prop.component=tile.map_reshaped_original(:)

prop.lpr_x = lpr_x(:);
prop.lpr_y = lpr_y(:);
prop.lpr_z = lpr_z(:);
% Detector_tilt = RTM.Rx(MicroscopeData.SampleTilt);



rot = rotation.byEuler(eulers);
% rot=inv(g_camera')*inv(Detector_tilt)*inv(rot)
ebsd_map = EBSD(rot, phases(:), CS_List,prop);
%%
% ebsd_original = ebsd_map;
end