function [PC_pattern_Astro,PatternInfo] = PC_OI_to_Astro(PC_pattern_OI,EBSP_size)
%PC_OI_TO_ASTRO Converts OI pattern centre to AstroEBSD pattern centre

PatternInfo.ScreenWidth=EBSP_size(1);
PatternInfo.ScreenHeight=EBSP_size(2);

%as per https://github.com/oinanoanalysis/h5oina/blob/master/H5OINAFile.md#ebsd-data
% Pattern Center X		H5T_NATIVE_FLOAT	(size, 1)	Pattern center X position scaled to the width of the image. This means that an X value of 0.5 is in the middle on the horizontal axis of the image. The origin is in the bottom left corner.
% Pattern Center Y		H5T_NATIVE_FLOAT	(size, 1)	Pattern center Y position scaled to the width of the image. Note that for a non-square image a Y value of 0.5 is not in the center of the vertical axis of the image. The origin is in the bottom left corner.
% Detector Distance		H5T_NATIVE_FLOAT	(size, 1)	Detector distance scaled to the width of the image.

%for Which Way is Up?
% PCx is measured from the left border of the EBSP in units of the pattern width (parallel and in the same direction to Xd in Fig. 2A). This is parallel to our sample tilt axis.
% PCy is measured from the top border of the EBSP in units of the pattern height (parallel but in the opposite direction to Yd in Fig. 2A).
% DD as the detector distance L normalised with respect to the pattern height.

PC_pattern_pixels_OI=PC_pattern_OI*PatternInfo.ScreenWidth;
PC_pattern_Astro=[PC_pattern_pixels_OI(1)/PatternInfo.ScreenWidth,1-(PC_pattern_pixels_OI(2)/PatternInfo.ScreenHeight),PC_pattern_pixels_OI(3)/PatternInfo.ScreenHeight];

end

