#!/usr/local_rwth/bin/zsh

#SBATCH -o /rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/fe6ta7/20kV/mg_multi_%j.out
#SBATCH -e /rwthfs/rz/cluster/hpcwork/mz071159/EMsoft_work/fe6ta7/20kV/mg_multi_%j.err

sbatch __mc_sh__ -W
wait

sbatch __master_sh__
