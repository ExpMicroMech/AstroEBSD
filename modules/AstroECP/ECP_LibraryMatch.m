function [G_Refined_SO3,PH_SO3,G_SO3,PH_Library]=ECP_LibraryMatch(ECP_Pat,template_library,library_G,screen_int,Settings_CorX,SettingsXCF,RTM,PC_in)
% Match against the library

[ text_out ] = pNameClock;
fprintf('%s: Pattern match started.\n',text_out); 

[ pattern1 ] = EBSP_BGCor( ECP_Pat.ECP_Pat_BG, Settings_CorX);

[EBSD_geom ] = EBSP_Gnom( RTM, PC_in); %

% [ pattern1 ] = EBSP_BGCor( thisProcdPattern, []);
[Pat_Ref_r,XCF_data_fill] = refine_prep(pattern1,SettingsXCF,RTM);

usePar=1; %use the parallel loop for the library testing
[G_SO3,PH_Library,~,~]=fLibraryTest(template_library,library_G,Pat_Ref_r.FFT,SettingsXCF,XCF_data_fill,usePar);
%refine based upon the SO3 found solution
[G_Refined_SO3,regout_RSO3] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,G_SO3,SettingsXCF,screen_int,screen_int.isHex,RTM);
PH_SO3=regout_RSO3(4);
[ text_out ] = pNameClock;
fprintf('%s: Pattern match completed.\n',text_out); 


end
