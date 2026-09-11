function MTEXversion = getMtexVersion()
    try
      fid = fopen(fullfile(mtex_path ,'VERSION'),'r');
      MTEXversion = fgetl(fid);
      fclose(fid);
      fprintf(' %s  ', MTEXversion);
      % regexp(MTEXversion,expression,'split')
    catch
      MTEXversion = 'MTEX';
    end
end