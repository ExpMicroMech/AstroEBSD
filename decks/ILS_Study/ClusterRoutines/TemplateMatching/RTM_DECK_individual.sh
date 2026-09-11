#!/usr/bin/zsh

### Job name
#SBATCH --job-name=RTMOPtimization

### File / path which STDOUT will be written to, %J is the job ID
#SBATCH --output=multiRTMCopies/logs/MATLAB_SERIAL.%J
#SBATCH -e multiRTMCopies/logs/MATLAB_ERROR_%j.err

### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --nodes=1
#SBATCH --time=23:59:00
##SBATCH --time=00:59:00
#SBATCH --mem=128000M
##SBATCH --mem=20000M
##SBATCH --cpus-per-task=1
#SBATCH --account=rwth1308
##SBATCH --partition=devel
#SBATCH --ntasks=50
##SBATCH --ntasks=4
### Change to the work directory, if not in submit directory
#cd "/work/mz071159/astro_skripts/decks/PCA"

### Load the required modules
module load MATLAB/2023b



## WARNING: SOFAR WE ONLY SPECIFY THE ROTATION MATRICES DIRECTLY IN THE MATLAB SCRIPT.



##
# Define the original MATLAB script
original_script="RTM_PCA_CommercialHough_mapscripted.m"
scriptFolder="multiRTMCopies"


# Define replacement strings
## slashes need to be escaped !!
## Folder setups
# phaseSubFolder="NbCo_MatchingPhases_LPMeasured" #copy the pha files you want to use into a folder in your Astro EBSD phases/phasesfiles/(newFolder) location
phaseSubFolder="Ni"
## make sure to escape slashes use syntax: \/
# EBSDFolder="\/hpcwork\/mz071159\/EBSD\/NbCoUBC_48Nb52Co_repolishOPS_newCut_vert\/Area2\/Map1"
# EBSDFolder="\/hpcwork\/mz071159\/EBSD\/LaNi2_LaNi2 vibropol_Area 10_EBSD+EDS\/Map4" 
EBSDFolder="\/hpcwork\/mz071159\/EBSD\/NiCu\/Area6"
EBSDMap="NiCuA01_DC_1.2a_Rest_sphericalIndentation_Area 6_EBSD+EDS Map 5" #
# EBSDMap="NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area 2_EBSD+EDS Map 1"
# PatternFolder="\/hpcwork\/mz071159\/EBSD\/NbCoUBC_48Nb52Co_repolishOPS_newCut_vert\/Area2\/Map1\/_EBSD\/vt0,1"
PatternFolder="\/hpcwork\/mz071159\/EBSD\/NiCu/Area6"
pccalibfile="Ni"
## Mode setups , run full tile or not, use initial eulers or not, use fitted pc or not
PCA_Mode="Full_Map" # alternatively: Full_Map / VMoutput / Tiles , whether to load all tiles or only roun the first one. For now it is hard coded a bit.
inputEulerMode="commercial" # commercial or BLANK
fittedPC=1 # whether or not the PC fit was determined before "-1"-> geometric , "0" -> let skript decide, ">1" ->assign pattern based on files in h5 file.
selectPCAuto=0 #22 for 51Nb-49Ni_v2 Area 6 #If -1, set a static pc based on commercial data.  if 0, it will select the data based on best ncc and deviation. Otherwise, it will select the number from the h5 file. insereted
PCResolution=514 # select from which intital PC resolution to choose. Typical values: 514,256,128,54. Warning: This might not be in your dataset.
refineLPR=0 #if 0 do not refine LPR
refineFFT=0  
plotOnly=0 #only replot old results or run again ?
refineOnly=1 #only do refinement
gradients=0
maxiter=300
## pcREfinementFromotherFile
pcrefineshfit=1
pcrefineFolder="\/hpcwork\/mz071159\/EBSD\/NiCu\/Area1" 
pcRefineH5File="NiCuA01_DC_1.2a_Rest_sphericalIndentation_Area 1_EBSD+EDS Map 1_reindex_PCRefinement___JOBID__gradients"
pcRefineOSCFile="NiCuA01_DC_1.2a_Rest_sphericalIndentation_Area 1_EBSD+EDS Map 1_reindex"
# /rwthfs/rz/cluster/hpcwork/mz071159/EBSD///NbNi_Nb-49Ni_V2_september_Area 8_EBSD+EDS Map 1/PatternCentrePlane

# Define the new MATLAB script name
rotationMatrixPseudo="rtmConfigRefine.RTM.Rz(0*degree)" #rtmConfigRefine.RTM.Rz(60*degree)" # for trigonal pseudosymmetry
# rotationMatrixPseudo="rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree)" ##for cubic - tetragonal pseudosymmetry
# rotationMatrixPseudo=""
new_script="multiRTMCopies/RTM_PCA_CommercialHough_map${SLURM_JOB_ID}.m"
script_name="RTM_PCA_CommercialHough_map${SLURM_JOB_ID}"



# Copy the original script to a new one
cp "$original_script" "$new_script"
newFolder="multiRTMCopies"
# MultiPCABashCopies







# Use sed to replace placeholders with specified strings
# Folder setups
sed -i "s/__EBSDFOLDER__/$EBSDFolder/g" "$new_script"
sed -i "s/__EBSDMAP__/$EBSDMap/g" "$new_script"
sed -i "s/__PATTERNFOLDER__/$PatternFolder/g" "$new_script"
sed -i "s/__PHASESUBFOLDER__/$phaseSubFolder/g" "$new_script"

# What to run 
sed -i "s/__PLOTONLY__/$plotOnly/g" "$new_script"
sed -i "s/__LPRMODE__/$refineLPR/g" "$new_script"
# sed -i "s/__PLOTONLY__/$plotOnly/g" "$new_script"
sed -i "s/__PCFITTED__/$fittedPC/g" "$new_script"
sed -i "s/__SELPCAUTO__/$selectPCAuto/g" "$new_script"
sed -i "s/__PCAMODE__/$PCA_Mode/g" "$new_script"
sed -i "s/__EULERMODE__/$inputEulerMode/g" "$new_script"
sed -i "s/__ROTATIONMATRIXPSEUDO__/$rotationMatrixPseudo/g" "$new_script"
sed -i "s/__MAXITER__/$maxiter/g" "$new_script"
sed -i "s/__PCCALIBFILE__/$pccalibfile/g" "$new_script"
sed -i "s/__REFINEONLY__/$refineOnly/g" "$new_script"
sed -i "s/__REFINEFFT__/$refineFFT/g" "$new_script"
sed -i "s/__GRADIENTS__/$gradients/g" "$new_script"

# PC refinement from other file
sed -i "s/__PCREFINESHIFT__/$pcrefineshfit/g" "$new_script"
sed -i "s/__PCREFINEFOLDER__/$pcrefineFolder/g" "$new_script"
sed -i "s/__PCREFINEANALYSIS__/$pcRefineH5File/g" "$new_script"
sed -i "s/__PCREFINEOSC__/$pcRefineOSCFile/g" "$new_script"





# cd $newFolder

matlab -nodisplay -nodesktop -nosplash -logfile multiRTMCopies/logs/job.${SLURM_JOB_ID}.txt <<EOF
cl=parcluster('local')
% objArray=findJob(cl,'running')
%instead go to /home/.../.matlab/local_cluster_jobs and delete there
% delete(objArray)
% objArray = findJob(cl,'pending')
% delete(objArray)
tmpdirforpool=${SLURM_JOB_ID}
%## tmpdirforpool=tmpdir
mkdir(tmpdirforpool)
cl.JobStorageLocation=tmpdirforpool
cl.NumWorkers=str2num(getenv('SLURM_NTASKS'))
parpool(cl,str2num(getenv('SLURM_NTASKS')))
cl; poolobj=gcp;
disp(['%-- parpool working with (' num2str(poolobj.NumWorkers) ') workers --%']);
disp(['%-- parpool working here (' poolobj.Cluster.JobStorageLocation ')  --%']);
% now run the 'dot-m-file.m' script
disp("${new_script}")

run("${new_script}")
delete(gcp)
quit();
EOF