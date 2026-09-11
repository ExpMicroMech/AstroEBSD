function densities = plotAtomsAlongDirection(atoms,millerPlot,neighbour,origin,InputUser,varargin)
   % col=color{i};
   unitCellExtent = 7;
   % figure,ax = plt.subplots(2,2)
   % OFFSET in multiples of the lattice spacing
   offset = 0;
   plot = 1;
     for i=1:size(varargin(:))
        if strcmp(varargin{i},'offset')
            offset=varargin{i+1};
        end      
        if strcmp(varargin{i},'DoNotPlot')
            plot = 0;
    end
    [~,uniqueIndices,~] =  unique([atoms.hklNew.x,atoms.hklNew.y,atoms.hklNew.z],'rows');
    % [~,uniqueIndices,~] = unique([atoms.all_atoms.x,atoms.all_atoms.y,atoms.all_atoms.z],'rows')
    uniqueAtoms  = atoms.hklNew(uniqueIndices);

    
%
       
    n_vec = millerPlot;
    extend = findMaxExtend(n_vec);
    [z_diff]=get_dist_from_line(n_vec,atoms.atomsUnique, uniqueAtoms);
    %%
    % temp = n_vec/max(n_vec.uvw)
    % 1/temp.w
    %%
   
    %%
    radii = (atoms.radiiNew(uniqueIndices));%,[1,ceil(extend)]);
    % newColors = repmat(atoms.colorsNew(uniqueIndices,:),[1,1,ceil(extend)]);
    % extends = repmat(extend,size(atoms.radiiNew));
    % extends = extends .* (1:ceil(extend));
    % extends = extends(:);
    % uniqueAtoms = repmat(uniqueAtoms,[1,ceil(extend)]);
    % uniqueAtoms = uniqueAtoms(:);
    % newColors = shiftdim(newColors,2);
    % newColors = reshape(newColors,[],4);
    radii = reshape(radii,1,[]);
    %
     rInPlane = sqrt(((radii).^2)-(z_diff).^2);
     %
     rInPlane = real(rInPlane);
     %%
     densities= [];

     for latticeSite = 1:size(atoms.atomsUnique,1)
         try
         [rInPlaneRow, rInPlaneCol] = find(rInPlane(latticeSite,:)~=0);     
         % imagesc(real(rInPlane));
         % rInPlane(latticeSite,rInPlaneCol);
         % colorbar
         %%
         %

         % rInPlane = rInPlane(uniqueIndices);
        %%
        atomsVec = vector3d(atoms.atomsUnique(latticeSite).x,atoms.atomsUnique(latticeSite).y,atoms.atomsUnique(latticeSite).z);
        atomSelection =  uniqueAtoms(rInPlaneCol);

        atomSelectionVec = vector3d(atomSelection.x,atomSelection.y,atomSelection.z);
        distance =  getCoordOnLine(atomsVec,n_vec,atomSelectionVec);
        % distanc<e = distance + extends(rInPlaneCol)
        %%
        corners = get_lattice_corners(n_vec,'n_cells',unitCellExtent,'extend_neg'); %for now magic number :(
       

    %%
        minBound = 0.01;
        maxBound = 0.99;
        % if any(abs(n_vec.uvw<0.05))
            distanceQuery = -1.8*max(corners(:)):0.02:1.8*max(corners(:));
            atomsReturn = LineEquation(atoms.atomsUnique(latticeSite),normalize(n_vec),distanceQuery);
           % distanceQuery = -1/(1-minBound)*max(corners(:)):0.01:(10/maxBound)*max(corners(:));
           %  atomsReturn = LineEquation(atoms.atomsUnique(latticeSite), n_vec,distanceQuery)

        % else
        %     while minBound >=0 || maxBound <=1
        % 
                % distanceQuery = -1/(1-minBound)*max(corners(:)):0.01:(1/maxBound)*max(corners(:)); % maximum would be for cubic cells sqrt(3) -> ~1.73 so 1.8 should be big enough
        %         % bounds.inShape(distanceQuery)
        % 
        %         % atomsReturn = LineEquation(atomsVec,n_vec,distanceQuery);
        %         atomsReturn = LineEquation(atoms.atomsUnique(latticeSite),n_vec,distanceQuery);
        %         minBound = max(min(atomsReturn.uvw))
        %         maxBound = min(max(atomsReturn.uvw))
        % 
        %     end
        % end
        % uvwNew = mod(atomsReturn.uvw,1);
        %%
         bounds=alphaShape(corners(1,:)',corners(2,:)',corners(3,:)',Inf);
        % atomsReturn = Miller(uvwNew(:,1),uvwNew(:,2),uvwNew(:,3),atomsReturn.CS);
        range = cat(3,(distance - rInPlane(latticeSite,rInPlaneCol)),(distance + rInPlane(latticeSite,rInPlaneCol)));
        upper =(distanceQuery>squeeze(range(1,:,2))');
        lower =(distanceQuery<squeeze(range(1,:,1))');
        merged = (upper+lower)==0;
        merged = sum(merged);
        merged = merged >0;
        inside = bounds.inShape(atomsReturn.x,atomsReturn.y,atomsReturn.z);
        densities = [densities, sum(merged(inside)) /sum(inside)];

    %%
    if plot == 1
        figure
        [X,Y] = getCircle(30,rInPlane(latticeSite,rInPlaneCol));
        for cNum = 1:size(X,2)

            patch(X(:,cNum)+distance(cNum),Y(:,cNum),atoms.colorsNew(uniqueIndices(rInPlaneCol(cNum)),1:3));
        end

        axis equal
        hold on
        exportgraphics(gcf,fullfile(InputUser.FigureFolder,...
                sprintf("AtomsAlongDirection_%s_%s_%s.png", n_vec.CS.mineral,char(n_vec), atoms.uniqueLabels(latticeSite) ) ) );
        plot(distanceQuery(inside),inside(inside),'k','linewidth',2)
        axis xy;
        % lims = distanceQuery(inside);
        % shift = abs(atoms.atomsUnique(latticeSite));
        % limDist = max(lims)-min(lims);
        % xlim([(min(lims)-limDist) , (max(lims)+limDist)]' );
        exportgraphics(gcf,fullfile(InputUser.FigureFolder,...
                sprintf("AtomsAlongDirection_%s_%s_%s_wLims.png", n_vec.CS.mineral,char(n_vec), atoms.uniqueLabels(latticeSite) ) ) );
         end
     catch
        warning('Direction probably outside unit Cell')
        densities = [NaN];
     end
     end
     close all
end