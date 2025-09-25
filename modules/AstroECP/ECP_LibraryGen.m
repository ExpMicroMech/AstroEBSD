function [library_G,template_library] = ECP_LibraryGen(screen_int,RTM_info,RTM,SettingsXCF,PC_in)
%ECP_LIBRARYGEN Generate the pattern matching library

[ text_out ] = pNameClock;
fprintf('%s: Pattern library generating.\n',text_out); 
cs_phase=loadCIF(RTM_info.cif_file);
[ library_G ] = SO3_rotmat_gen( cs_phase,RTM.Sampling_Freq);
% generate the library
[EBSD_geom ] = EBSP_Gnom( RTM, PC_in); %
[ template_library ] = Library_Gen(EBSD_geom,screen_int,RTM_info.isHex,library_G,2,SettingsXCF);

[ text_out ] = pNameClock;
fprintf('%s: Pattern library generated.\n',text_out); 


end