function mp = MTEX_mb_fix(mb_length)
% a quick function to help MTEX plots and fix their lengths if you want

f=gcm;
mp=getappdata(f.currentAxes,'mapPlot');
mp.micronBar.length = mb_length; % change length - in um

end
