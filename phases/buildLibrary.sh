#!/usr/bin/zsh

### Job name
#SBATCH --job-name=Fit_PatternCentre

### File / path which STDOUT will be written to, %J is the job ID
#SBATCH --output=logs/MATLAB_SERIAL.%J
#SBATCH -e logs/MATLAB_ERROR_%j.err
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --mem-per-cpu=3200M
###SBATCH --cpus-per-task=1
#SBATCH --ntasks=3


### Change to the work directory, if not in submit directory
#cd "/work/mz071159/astro_skripts/decks/PCA/minimization"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berners@imm.rwth-aachen.de

 
#SBATCH --partition=devel
####SBATCH --partition=c23mm
####SBATCH --account=rwth1308
###SBATCH --constraint Rocky9  

### Load the required modules
module load MATLAB/2023b

### start non-interactive batch job
matlab -nodisplay -nodesktop -nosplash -logfile logs/job.${SLURM_JOB_ID}.txt <<EOF
# disp(['%-- parpool working with  workers --%']);
# run BuildNewPhases.m

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
run BuildNewPhases.m
delete(gcp)

quit();
EOF