function cmap_div = returnCmapDiv
warning ('off','all');
 [cmap_div]=cbrewer('div','RdBu',128);
%  [cmap_div] = cbrewer('')
cmap_div(cmap_div<0)=0;
warning ('on','all');
end