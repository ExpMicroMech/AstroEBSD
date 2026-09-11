# Analysis Workflow:

## Simulation procedure

- Install EMSoft 
  
  - [Home · EMsoft-org/EMsoft Wiki · GitHub](https://github.com/EMsoft-org/EMsoft/wiki)[Home · EMsoft-org/EMsoft Wiki · GitHub](https://github.com/EMsoft-org/EMsoft/wiki)
  
  - **[Anmelden · GitLab](https://git.rwth-aachen.de/imm_git/emsoftfork)** (some small bug fixes here. Uli also wrote a installation skript. (in IMMG_git/EBSDSim)

- Use python "wrapper" to run multiple simulations:
  
  - [Anmelden · GitLab](https://git.rwth-aachen.de/sfb1394/ebsdsim/-/tree/modify_python_automation?ref_type=heads)
    
    - There is a runMultipleScans.sh shell script that can loop over multiple structure files.
      
      
    - ### Build AstroEBSD Phase Library

- Move all Output structures to one folder with a unique name. (They way this is typicall set up on the hpc cluster the end up in some folder in hpcwork/<your_TIM>/<EMSoftWork>).

- Go to  a recent clone of AstroEBSD 
  
  - In BuildPhases.m change the rootDir to the folder you just chose to move all your simulations in.
    
    - This will automatically build .pha files  in phases/phasefiles and further copy the correct cif to phases/cifs. It will also move the most recent .h5 in a specific subfolder to dynamic_templates and reference the dynamic template and the cif in the respective .pha files.
    
    - You might need to change the Astro and Mtex locations (running on the hpc we have these specified in each skript individually, while on Windows we typicaly just use the startAstroAndMtex.m skript.)
  
  - In AstroEBSD/phases run buildLibrary.sh (just source buildLibrary.sh on Linux)
  
  - On Windows you can just run the BuildNewPhases.m

- COPY all phases in the phasefiles folder, that you want to use for pattern matching to one folder with a unique name. We have adjusted the pattern matching skripts in such a way, that they can just load all phases from this folder. It is important to __COPY__ these files only, as Astro will need to acces the phase still from the base phasefile folder. This is implemented this way, to ensure compatibility between versions.

## Pattern Matching Workflow

### PCA

- In principle the refinment can also run on complete maps. However, to reduce the amount of data analysed and increase the signal to noise ratio of EBSD patterns, we use the PCA approach as described in https://www.sciencedirect.com/science/article/pii/S0304399119302773?via%3Dihub 
- To run the PCA go to astro_skripts (separate git repository) (https://github.com/potsci/astro_skripts)
- astro_skripts/decks/PCA
- pca_deckMulti.sh
- Here you can specify Locations of the EBSD files you want to run the PCA on. (The .osc and .up2 or .up1 files have to have the same name and be in the same folder)
- You can also specify the most common paramaters of the PCA,
- such as the crop factor and the variance retained.
  - You have to test which crop-factor works best. (It is reasonably fast with about 2GB of Patterns in RAM).
  - variance typically works well with variance ="0.1" (might be adjusted to smaller values if there are lots of orientaiton gradients.)
- You can specify, whether a realBG is substracted and whether this is from the same map or not (if from same map set mapBG=1)
- Other parameters can be modified in PCA_ded_Mg_remote_cluster.m (which we use the shell script to make a copy of). These paramters generally do not need to be touched.
- Leave pieces with \__VariableName\_\_ untouched. They will be overwritten by the .sh skritp.

### PC refinement.

- Go to decks/PCA/minimization

- run FitPatternCenter.sh (sbatch FitPatternCenter.sh)

- Here you need to specify the phaseFolder (the one created under [Build AstroEBSD Phase Library](#Build AstroEBSD Phase Library)). 

- Further specify the EBSD map that you want to use for calibration 
  
  - This should be a fairly large map that contains multiple grains. 
  
  - For our study we went with roughly 1x1 mm.
  
  - Select a number of grid-poins to run the refinement on (more than 10 seems to work great). 
  
  - pSerachSize allows to find with in a certain circumference a region that has the highest IQ value pattern.

- Leave the rest of parameters untouched.
  
  - If fine tuning needs to be done, this can be by modifying PatternCalibrationDeck and IndexingConfigFile.m hoewever. Take care to not crash any parameters.

- This will run on all structures that are in the specified folder. Later you can choose the one that best minimizes all PC values to a plane.

### Actual refinment of map to investigate (combine pca and PC refinement.)

 As this is designed to be working on multiple  maps, we have also written a neat shell script that handles mutliple ebsd map and structure setups.

- go to decks/PCA

- The skript you want to run is RTM_DeckLoop.sh
  
  - This adjust parameters in the RTM_DECK_indivual_Template.sh skript (a copy will be generated for each map and Phase folder combination)
  
  - Specify the PhaseFolders and folders where you ran the PCA on (EBSDFolder, sublevel to SampleFolder)
  
  - Specifine the location of the File you used for PC refinmenet.
  
  - The example inputs should be mostly self explanatory.
  
  - Further Down you can specify, which variance retained your PCA analysis had. This is necessary, if you want to swithc the vt of the pCA
  
  - rotationMatrixPseudo can be used to distinguish between different pseudo symmetries. just have rtmConfigRefine.RTM.Rz(0*degree) there if you do not want to use this.
  
  - Don't touch the loop this will just loop over each ebsd map it finds and runt he pattern matching.

- you wil get a .h5 file as an output that contains details about the refinment and the results:

- ![](C:\Users\berners\AppData\Roaming\marktext\images\2026-01-23-15-30-20-image.png)

Each structure will be inside RTM and then you have mutliple folders for each tile (1-4)

- Depending on which kind of analysis you use the name of the folders in here may vary.

- This file can be used by the various analysis skripts that we have built:

- NOTE: There is a python skript to remove the params folders, as these slow down loading quite a bit. (right now it sits in lukas branch of the sfb1394 rdm repo. But we will move this one #TODO) splitAstroH5Files.ipynb

## ## Analysis Workflow
