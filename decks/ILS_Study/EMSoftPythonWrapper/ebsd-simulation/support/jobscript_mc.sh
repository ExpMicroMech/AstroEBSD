#!/usr/local_rwth/bin/zsh

####Job Name
#SBATCH --job-name=__jobname__
####Outputfile
#SBATCH -o __filepath__%j.out
#SBATCH -e __filepath__%j.err
#SBATCH --account=__acc__

###execution time
#SBATCH --time=__hours_mc__

##Memory need / node
#SBATCH --mem=__memory_mc__

#number of threads
#SBATCH --cpus-per-task=__num_cpu_mc__
#SBATCH --gres=gpu:1
#SBATCH --partition=c23g

##notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=__email__
##Linux system
##SBATCH --constraint Rocky9  


## reg shell commands
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
##### export CUDA_VISIBLE_DEVICES=<GUID>
# echo $CUDA_VISIBLE_DEVICES$ >> __filepath__%j.out
####export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## module load CUDA/11.8.0
## module load GCC/12.3.0
## #module load Clang/16.0.6
## module load LLVM/14.0.3
# module load intel-compilers/2024.2.0 
# module load intel/2024a
# module load impi/2021.13.0
# module load imkl-FFTW/2024.2.0
# module load GCCcore/13.3.0
# module load LLVM/18.1.8
# module load Clang/18.1.8
# module load GCC/13.3.0

# # module load OpenMPI/5.0.3
# module load CMake/3.31.7

# # module load Python/3.9.6
# module load HDF5/1.14.5
# module load tbb/2021.13.0
# module load intel-compilers/2024.2.0 
# # module load intel/2024a
# module load impi/2021.13.0
# module load imkl-FFTW/2024.2.0
# module load CUDA/12.6.3
export PATH=~/EMsoftBuild/Release/Bin:$PATH
# export PATH=~/EMsoftBuildRelease//Release/Bin:$PATH
export PATH=~/EMsoftSDK/Debug:$PATH
# export PATH=~/EMsoftSDK/Release:$PATH
export PATH=/home/mz071159/EMsoftSuperbuild/tbb-2020.1-lin/tbb:$PATH
# export OPENBISPAT="$pat-berners-231205084303489x2EE3A5E38D34E2E58B41BB31C33379D3"
# export PATH=~/EMsoftSDK:$PATH
# export PATH=/rwthfs/rz/SW/UTIL/clang/13.3.0/lib:$PATH
# export PATH=/rwthfs/rz/SW/UTIL/clang/13.3.0/lib64:$PATH
# export PATH=/usr/local_rwth/sw/gcc/13.3.0/bin/gfortran:$PATH
# export PATH=/usr/local_rwth/sw/gcc/13.3.0/lib64:$PATH
# export PATH=/lib64:$PATH
# export LD_LIBRARY_PATH=/usr/local_rwth/sw/gcc/13.3.0/bin/gfortran:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/cvmfs/software.hpc.rwth.de/Linux/RH8/x86_64/intel/skylake_avx512/software/GCCcore/13.3.0/bin/gfortran:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/cvmfs/software.hpc.rwth.de/Linux/RH8/x86_64/intel/skylake_avx512/software/GCCcore/13.3.0/lib64:$LD_LIBRARY_PATH
# # 
# export LD_LIBRARY_PATH=/rwthfs/rz/SW/UTIL/clang/13.3.0/lib64:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/usr/local_rwth/sw/gcc/13.3.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/etc/OpenCL/vendors:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/mz071159/EMsoftSDK/Debug/fftw-3.3.8/lib64:$LD_LIBRARY_PATH
# echo $LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# which EMMCOpenCL
srun EMMCOpenCL __nml_file__
# srun EMOpenCLinfo

# echo; export; echo; nvidia-smi; echo

# $CUDA_ROOT/extras/demo_suite/deviceQuery -noprompt
