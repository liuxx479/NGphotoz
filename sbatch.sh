#! /bin/bash
#SBATCH --nodes=2
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH -C haswell
#SBATCH -t 0:30:00
#SBATCH -J pzNG
#SBATCH -o pzNG.o%j
#SBATCH -e pzNG.e%j
#SBATCH --qos=debug ##regular
#SBATCH -A m1727
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jialiu@berkeley.edu

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
#export OMP_NUM_THREADS=2

srun python /global/u1/j/jialiu/NGphotoz/map_stats.py