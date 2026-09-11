#!/usr/local_rwth/bin/zsh

####Job Name
#SBATCH --job-name=__jobname__
####Outputfile
#SBATCH -o __filepath__%j.out
#SBATCH -e __filepath__%j.err
#SBATCH --account=__acc__

###execution time
#SBATCH --time=__hours_master__

##Memory need / node
#SBATCH --mem=__memory_master__

#number of threads
#SBATCH --cpus-per-task=__num_cpu_master__

##notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=__email__
#SBATCH --constraint Rocky9  



## reg shell commands
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# module load intel-compilers/2024.2.0 
# module load intel/2024a
# module load imkl-FFTW/2024.2.0
# module load GCCcore/13.3.0
# module load LLVM/18.1.8
# module load Clang/18.1.8
# module load GCC/13.3.0
# module load CUDA/12.6.3
# module load OpenMPI/5.0.3
# module load CMake/3.31.7

# # # module load Python/3.9.6
# module load HDF5/1.14.5
# module load tbb/2021.13.0
export PATH=~/EMsoftBuild/Release/Bin:$PATH
# export PATH=~/EMsoftBuildRelease//Release/Bin:$PATH
export PATH=~/EMsoftSDK/Debug:$PATH
# export PATH=~/EMsoftSDK/Release:$PATH
export PATH=/home/mz071159/EMsoftSuperbuild/tbb-2020.1-lin/tbb:$PATH
# export PATH=/rwthfs/rz/SW/UTIL/clang/11.1.0/lib:$PATH
# export PATH=/rwthfs/rz/SW/UTIL/clang/11.1.0/lib64:$PATH
# export PATH=/usr/local_rwth/sw/gcc/11.2.0/bin/gfortran:$PATH
# export PATH=/usr/local_rwth/sw/gcc/11.2.0/lib64:$PATH
# export LD_LIBRARY_PATH=/usr/local_rwth/sw/gcc/11.2.0/bin/gfortran:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/rwthfs/rz/SW/UTIL/clang/11.1.0/lib64:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/usr/local_rwth/sw/gcc/11.2.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/etc/OpenCL/vendors:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/mz071159/EMsoftSDK/Debug/fftw-3.3.8/lib64:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH
srun   EMEBSDmaster __nml_file__

