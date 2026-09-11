# CompareMisalignmentETC

This Folder contains multiple scripts to plot the Delta NCC plots and extract the data from refinments of multiple crystal structures.

While BatchPlotILSpacings is supposed to the differences of the Interlayer Spacing template refinements.
Here you can specify the folders in a similar manner, as in the rtm_deck_loop.sh skripts.

The BatchUpdateMaps.m is used to "clean" the scatter of refinmeent from the template matching and add some variables to the hdf5 file, that will be used in the python skripts.

## CompareMisalignment
- These two .m skripts are used to generate the Figure 7 in the manuscript.

## GnomonicDifferences.m 
- Can be used to plot the difference plot of two reprojected master patterns (also contains some additional workflows to compare the radon transform and other similarity metrics such as SSIM)

## PlotRefinementResultsFromH5FileNbNi.m
- Plot refinement results from the h5 files (e.g. interlayer-spacing NCC values)
- This is mostly replaced by the python skripts now.

## RTM_PCA_CommericalHough_map_OnlyFromH5.m
- Can be used, to plot the refinement results and the reprojections, that are used to extract the band profiles. (by default we turn these plottings of in the cluster based skripts. So if you want to have these results you need to run this deck.)



