#!/bin/bash

#SBATCH --job-name=kmeans-MPI
#SBATCH --nodes=2 
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=4GB 
#SBATCH --time=00:30:00 
#SBATCH --verbose

module purge 
module load  amber/openmpi/intel/20.06 

mpirun -np 50 --oversubscribe ./kmeansMPI 10 42000 784 train.csv
