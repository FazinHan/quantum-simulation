#!/bin/bash
#SBATCH --job-name=jeff
#SBATCH --ntasks=2
#SBATCH --output=qs.out
#SBATCH --time=4-00:00:00

source ~/.bashrc
time mpiexec.hydra -n $SLURM_NTASKS python parallel_julia.py
#python parallel_julia.py 4
#python parallel_julia.py 6
#python parallel_julia.py 8
#push
