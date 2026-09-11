#!/usr/bin/zsh

### Job name
#SBATCH --job-name=PCA_Matlab

### File / path which STDOUT will be written to, %J is the job ID
#SBATCH --output=logs/MATLAB_SERIAL.%J
##SBATCH -e logs/MATLAB_ERROR_%j.err
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --nodes=1
#SBATCH --time=5:02:00
#SBATCH --mem=128000M
#SBATCH --cpus-per-task=12
#BATCH --account=rwth1308
###SBATCH --partition=devel
#SBATCH --ntasks=1


### Change to the work directory, if not in submit directory
#cd "/work/mz071159/astro_skripts/decks/PCA"

### Load the required modules
module load MATLAB/2023b


##
# Define the original MATLAB script
original_script="PCA_deck_Mg_remote_cluster.m"
scriptFolder="MultiPCABashCopies"


# Define replacement strings
phaseName="Nb518Ni482_joubert_19.3kV_dmin0.020"
## slashes need to be escaped !!
EBSDFolder="\/hpcwork\/mz071159\/EBSD\/MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3" 
# EBSDFolder="\/hpcwork\/mz071159\/EBSD\/58Nb-42Ni_v2_septemberbatch\/Area7\/Map1"
#/rwthfs/rz/cluster/hpcwork/mz071159/EBSD/MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3

# EBSDFolder="\/hpcwork\/mz071159\/EBSD\/HfNi\/"
EBSDMap="MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3"
# EBSDMap="Scan12"
cropfactor="2";
variance="0.1";
realBG="0" ## will be set to 1 internally, if mapBG = 1
mapBG="0"
bgFolder="\/hpcwork\/mz071159\/EBSD\/LaNi2_LaNi2 vibropol_Area 10_EBSD+EDS\/Map5"
bgMap="LaNi2_LaNi2 vibropol_Area 10_EBSD+EDS Map 5"

# Define the new MATLAB script name
new_script="MultiPCABashCopies/PCA_deck_Mg_copy${SLURM_JOB_ID}.m"
script_name="PCA_deck_Mg_copy${SLURM_JOB_ID}"

# Copy the original script to a new one
cp "$original_script" "$new_script"
newFolder="MultiPCABashCopies"
# MultiPCABashCopies

# Use sed to replace placeholders with specified strings
sed -i "s/__PHASENAME__/$phaseName/g" "$new_script"
sed -i "s/__EBSDMAP__/$EBSDMap/g" "$new_script"
sed -i "s/__EBSDFOLDER__/$EBSDFolder/g" "$new_script"
sed -i "s/__CROPFACTOR__/$cropfactor/g" "$new_script"
sed -i "s/__VARIANCE__/$variance/g" "$new_script"
sed -i "s/__USEREALBG__/$realBG/g" "$new_script"
sed -i "s/__USEMAPBG__/$mapBG/g" "$new_script"
sed -i "s/__BGFOLDER__/$bgFolder/g" "$new_script"
sed -i "s/__BGMAP__/$bgMap/g" "$new_script"
# cd $newFolder

matlab -nodisplay -nodesktop -nosplash -logfile MultiPCABashCopies/logs/job.${SLURM_JOB_ID}.txt <<EOF
cl = parcluster('local')
tmpdirforpool = ${SLURM_JOB_ID}
## tmpdirforpool = tmpdir
mkdir(tmpdirforpool)
cl.JobStorageLocation = tmpdirforpool
cl.NumWorkers = str2num(getenv('SLURM_NTASKS'))
parpool(cl,str2num(getenv('SLURM_NTASKS')))
cl; poolobj = gcp;
disp(['%-- parpool working with (' num2str(poolobj.numworkers) ') workers --%']);
disp(['%-- parpool working here (' poolobj.Cluster.JobStorageLocation ')  --%']);
% now run the 'dot-m-file.m' script
disp("${new_script}")

run("${new_script}")
delete(gcp)
quit();
EOF