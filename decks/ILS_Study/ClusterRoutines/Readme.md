# Routines to run matlab PCA extraction and pattern matching on a HPC cluster.

You have three different folders, for different puposes (the folder names should be self explanatory)

## PCPlaneFitting
The shell skript in PCPlaneFitting does template matching on a sub-grid of an ebsd-map, that later can be a plane fitted onto (affine projection)
NOTE, that here we have two matlab files that you need.
The most basic parameters like file names, set up of the grid etc. can be controlled from the .sh file. However, for more detailed parameters please also have a look at the .m files

## PCA
Again this is just a shell wrapper, to run the pca as a batch job. This essentially helps , if the normal pca skripts do not have enough memory on your local machine -> you can go to a hpc, where there are typically nodes with lots of RAM.
Here you can also specify basic parameters in the .sh file. For more details have a look at the .m files

## Template Matching
Here are 3 shell skripts that call a matlab template file. The main file is RTM_DeckLoop.sh. This can run over multiple EBSD maps. (that you need to specify in the shell skript.) The "empty" variables in RTM_DECK_individual_TEMPLATE.sh will get populated by this skript. Again. Basic parameters can be set in the main shell skript. For details go to the (..)_mapscripted.m file.

The RTM_DECK_individual.sh file is used to run on single files instead of multiple folders
