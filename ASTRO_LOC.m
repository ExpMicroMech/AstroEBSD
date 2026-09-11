function Astro_loc=ASTRO_LOC()
   
   currentFileName = mfilename('fullpath');
   Astro_loc = fileparts(currentFileName);
   % Astro_loc = Astro_loc(1)
end