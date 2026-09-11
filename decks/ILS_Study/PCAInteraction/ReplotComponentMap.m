%% L. Benrers 26.08.2025
% replotting the component map to make it more easily visible which pattern
% belongs to where.
%
%
% PCAFolder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\51Nb-49Ni_v2\EBSD\Area10\NbNi_Nb-49Ni_V2_september_Area 10_EBSD+EDS Map 1\_EBSD\vt0,05"
PCAFolder = "J:\Users\Berners_Lukas\clara\YNi2MaDinah\YNi2_S_D_repolish\Area3\_EBSD\vt0,1"
PCAFolder = "J:\Users\Berners_Lukas\clara\C15_2_comparisonExperiment\Area1\_EBSD\vt0,1"
tiles = dir(fullfile(PCAFolder,"\Tile*"))
% tiles = tiles(tiles.isdir)
maxCnum = 10;
[cmap] = cbrewer('qual', 'Paired',maxCnum+2,'cubic')
cmap(12,:)=[ 0.5 0.5 0.5]
cmap(1,:)=[ 0.5 0.5 0.5]
cmap(cmap<0)=0
cmap(cmap>1)=1
%%
for tNum = 1:numel(tiles)
   file= dir(fullfile(PCAFolder, tiles(tNum).name, "**/_PCA_.mat"))
   load(fullfile(file(1).folder,file(1).name),'VMOutput')
   %%

   if  VMOutput.NumPCA>20
        nTiles = ceil(VMOutput.NumPCA/maxCnum)
        f= figure
        f.Position = [100 100 670*ceil(nTiles/2) 1200]
        h= tiledlayout(2,ceil(nTiles/2),"TileSpacing",'tight')
       
   else
       f = figure
       h = tiledlayout(2,1,"TileSpacing",'tight')
       nTiles = 2
       f.Position = [100 100 1400 600]
   end
   for layoutNum = 1:nTiles
       nexttile(layoutNum)
       imagesc(VMOutput.filtered_map)
       colorbar()
       colormap(cmap)
       clim([(layoutNum-1)*maxCnum,(layoutNum)*(maxCnum)+2])
       axis off
       axis equal
       % axis xy
       
    nexttile(layoutNum)
    
   end
   exportgraphics(f,fullfile(file(1).folder,sprintf("_rPCA_wColorbar.png")))
   end
%%

    % cd(InputUser.ResultsDir)
    cmap1=cbrewer('qual','Dark2',80);
    
    %deal with a norming issue
    cmap1=cmap1-min(cmap1(:));
    while max(cmap1(:) > 1) %norm the range
        cmap1=cmap1./max(cmap1(:));
    end
    
    cmap=repmat(cmap1,1,1);

%%

%%
file= dir(fullfile(PCAFolder, "Results.mat"))
load(fullfile(file(1).folder,file(1).name),'tile')

%%

cmap1=cbrewer('qual','Set1',80);

%deal with a norming issue
cmap1=cmap1-min(cmap1(:));
while max(cmap1(:) > 1) %norm the range
    cmap1=cmap1./max(cmap1(:));
end

cmap=repmat(cmap1,100,1);
f=figure
imagesc(medfilt2(tile.map_reshaped)); colormap(cmap); axis off;
axis equal;
% % colorbar()
% % clim([1,cmap1./max(cmap1(:))])
exportgraphics(f,fullfile(file(1).folder,sprintf("_rPCA_wColorbar.png")))