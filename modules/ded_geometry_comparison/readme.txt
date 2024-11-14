Accompanying scripts for:

"Comparison of Kikuchi Diffraction Geometries in Scanning Electron Microscopy"
by Tianbi Zhang, Lukas Berners, Jakub Holzer, T. Ben Britton

------

Contents: 

(1) Reconstruction_from_pattern: takes a gnomonic pattern, reproject it on a stereogram (based on simulation) and inflate the pattern to obtain the full diffraction sphere. 
	(*) Patterns must be pre-indexed - alternatively, you may use decks/TKD_HDR_stereo_reproject
	(**) Stereogram is plotted on a 2001x2001 regular grid. This is chosen for higher quality of the reconstruction.

(2) Reconstruction_band_profiling: takes a gnomonic pattern, reproject it on a stereogram (based on simulation) and inflate the pattern to obtain the full diffraction sphere. A simulated stereogram is also created. Then the stereograms are approximated by spherical harmonic functions, and band profiles of {002}, {022} and {111} are extracted.
	(*) Patterns must be pre-indexed.
	(**) Stereogram is approximated and will appear more blurry than that reconstructed from (1). Thus it is only used for band profiling.
	(***) The class defined in detector.m is a dependency.

(3) Get_RKD_Pattern, v1 and v2
	(*) These scripts use a few custom-written dependencies, which are also included in this bundle.
	(**) v1 works better on smaller datasets (~400 patterns) and v2 on larger datasets. However, neither is optimized to cope with the maximum allowable array size of MATLAB.

(4) Reproject_to_pattern: generate a gnomonic Kikuchi pattern from a diffraction sphere.

(5) Plot_band_profiles: plot the band profiles obtained from (2) and generate a figure akin to Figure 5(a) on the paper.

------

Notes:

(1) To process the raw .h5 EBSD and TKD patterns, use scripts in /modules/ded_general/
	(*) To process the on-axis TKD patterns using the exposure fusion method, use decks/TKD_HDR_pattern_astro which incorporates functions in /modules/ded_tkd

(2) The following were not performed within a single script, but manually:
	(i) Averaging multiple reconstructed stereograms (both regular grid and spherical harmonics)
	(ii) Saving the band profiles for plotting.
	(iii)

We exported the reconstructed stereograms and averaged them manually.
