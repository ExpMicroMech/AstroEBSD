#!/usr/bin/zsh

### Job name
#SBATCH --job-name=Fit_PatternCentre

### File / path which STDOUT will be written to, %J is the job ID
#SBATCH --output=MultiPCFITScripts/logs/MATLAB_SERIAL.%J
#SBATCH -e MultiPCFITScripts/logs/MATLAB_ERROR_%j.err
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --nodes=1
#SBATCH --time=15:00:00
####SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=5200M
#SBATCH --cpus-per-task=1
####SBATCH --ntasks=4
#SBATCH --ntasks=20


### Change to the work directory, if not in submit directory
#cd "/work/mz071159/astro_skripts/decks/PCA/minimization"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berners@imm.rwth-aachen.de


###SBATCH --partition=devel
#SBATCH --partition=c23mm
##SBATCH --account=rwth1308
###SBATCH --constraint Rocky9  

### Load the required modules
module load MATLAB/2023b


original_script="PatternCalibrationDeck.m"
scriptFolder="MultiPCFITScripts"

original_config="IndexingConfigFile.m"

# Define replacement strings

phaseName='""' # gets overwritten if phaseFolder is not empty.
# if you want to use multiple phases, place them in a folder 
#inside AstroEBSD/phases/phasefiles/subfolder and specify the subfolder here
# otherwise put in '""' -> empty string
# phaseFolder='"YNi_FitPhases_29.3kV"'
phaseFolder='"FeSi"' ## needs to have single and doulbe quotes.
# phaseName="'Nb518Ni482_joubert_19.3kV_dmin0.020'"
# phaseName="'Nb496Ni504_joubert_19.5kV_dmin0.020'"
## slashes need to be escaped !!
# EBSDFolder="51Nb-49Ni_septemberbatch\/NbNi_Nb-49Ni_V2_september_Area 8_EBSD+EDS_Map_1" 
# EBSDFolder="58Nb-42Ni_v2_septemberbatch\/Area6"
EBSDFolder="FeSi\/MP1_S1_def0.7\/Area1"
# EBSDFolder="NbCoUBC_48Nb52Co_repolishOPS_newCut_vert\/Area1"
# EBSDFolder="49Nb-51Co_WeiLuo_HT2/Area1"
# EBSDFolder="Ni_VirtualECCI\/Area1"
# EBSDFolder="saaa
EBSDMap="FeSi_Macropillar1_Slice1_def0.7_Area 1_OIM Map 1"
# EBSDMap="NbCoUBC_49Nb51Co_WeiLuo_HT2_Area 1_EBSD+EDS Map 2"
gridPoints=10; # Will result in 10 points -> 10x10 grid.
kernelSize=3;
pSerachSize=5; ### this shoudl basically lead to a fixed grid.
useGradients=1;



# Define the new MATLAB script name
new_config="MultiPCFITScripts/IndexingConfigFile${SLURM_JOB_ID}.m"
config_name="IndexingConfigFile${SLURM_JOB_ID}"

new_script="MultiPCFITScripts/PatternCalibrationDeck${SLURM_JOB_ID}.m"
script_name="PatternCalibrationDeck${SLURM_JOB_ID}"
jobId="${SLURM_JOB_ID}"
# Copy the original script to a new one
cp "$original_script" "$new_script"

cp "$original_config" "$new_config"
newFolder="MultiPCABashCopies"
# MultiPCABashCopies

# Use sed to replace placeholders with specified strings
 
sed -i "s/__EBSDFOLDER__/$EBSDFolder/g" "$new_config"
sed -i "s/__EBSDMAP__/$EBSDMap/g" "$new_config"
sed -i "s/__PHASENAME__/$phaseName/g" "$new_config"
sed -i "s/__PHASEFOLDER__/$phaseFolder/g" "$new_config"
sed -i "s/__PSEARCHSIZE__/$pSerachSize/g" "$new_config"
sed -i "s/__KERNELSIZE__/$kernelSize/g" "$new_config"
sed -i "s/__GRIDPOINTS__/$gridPoints/g" "$new_config"
sed -i "s/__USEGRADIENTS__/$useGradients/g" "$new_config"
## Do the actual script
sed -i "s/__EBSDMAP__/$EBSDMap/g" "$new_script"
sed -i "s/__INDEXCONFIGFILE__/$config_name/g" "$new_script"
sed -i "s/__JOBID__/$jobId/g" "new_script"
### start non-interactive batch job
matlab -nodisplay -nodesktop -nosplash -logfile logs/job.${SLURM_JOB_ID}.txt <<EOF
# disp(['%-- parpool working with  workers --%']);

cl = parcluster('local')
tmpdirforpool = ${SLURM_JOB_ID}
## tmpdirforpool = tmpdir
mkdir(tmpdirforpool)
cl.JobStorageLocation = tmpdirforpool
cl.NumWorkers = str2num(getenv('SLURM_NTASKS'))
parpool(cl,str2num(getenv('SLURM_NTASKS')))
cl; poolobj = gcp;
disp(['%-- parpool working with (' num2str(poolobj.NumWorkers) ') workers --%']);
disp(['%-- parpool working here (' poolobj.Cluster.JobStorageLocation ')  --%']);
% now run the 'dot-m-file.m' script
disp("${new_script}")

run("${new_script}")
delete(gcp)

quit();
EOF