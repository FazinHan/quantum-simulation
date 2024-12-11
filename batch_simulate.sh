#SBATCH -ntasks 3

parallel python parallel_julia.py ::: 4 6
#python parallel_julia.py 4
#python parallel_julia.py 6
#python parallel_julia.py 8
push
