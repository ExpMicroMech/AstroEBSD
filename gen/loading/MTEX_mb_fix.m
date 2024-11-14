function mp = MTEX_mb_fix(mb_length)
% a quick function to help MTEX plots and fix their lengths if you want
% 
% Ben Britton - Feb 2024

f=gcm;
mp=getappdata(f.currentAxes,'mapPlot');
mp.micronBar.length = mb_length; % change length - in um

end
