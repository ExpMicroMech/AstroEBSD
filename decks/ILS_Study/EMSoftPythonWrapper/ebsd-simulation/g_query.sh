#!/usr/local_rwth/bin/zsh

####Job Name
#SBATCH --job-name=openclinfo
####Outputfile
#SBATCH --output=/rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/NbNi_dt0.250_c26.5_a4.9/20kV/emmopcenl%j_.out

#SBATCH --account=rwth1308

###execution time
#SBATCH --time=0:00:45

##Memory need / node
#SBATCH --mem=2000M

#number of threads
#SBATCH --cpus-per-task=1 
#SBATCH --gres=gpu:1

##notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.berners@rwth-aachen.de


# module load GCCcore/13.3.0
# module load LLVM/18.1.8
# module load Clang/18.1.8
# module load GCC/13.3.0
# module load CUDA/12.6.3
# # module load OpenMPI/5.0.3
# module load CMake/3.31.7

# # module load Python/3.9.6
# module load HDF5/1.14.5
# module load tbb/2021.13.0
# module load intel-compilers/2024.2.0 
# # module load intel/2024a
# module load impi/2021.13.0
# module load imkl-FFTW/2024.2.0
## reg shell commands
#SBATCH --constraint Rocky9 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PATH=~/EMsoftBuild/Release/Bin:$PATH
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
srun EMOpenCLinfo
