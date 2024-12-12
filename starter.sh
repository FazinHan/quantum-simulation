#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=fauseweh-zhu
#SBATCH --output=qs.out
#SBATCH --time=4-00:00:00
#SBATCH --ntasks-per-node=2
#SBATCH -A physics_engg
#SBATCH --mail-user=fizaan.khan.phy21@iitbhu.ac.in

time mpiexec.hydra -n $SLURM_NTASKS python parallel_julia.py