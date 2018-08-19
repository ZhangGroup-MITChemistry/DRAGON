#!/bin/bash

#SBATCH --job-name=lmps
#SBATCH -N 1
#SBATCH -n 14
#SBATCH --partition=mit
#SBATCH --no-requeue
#SBATCH --time=48:00:00
#SBATCH --export=ALL                                                           

module load gcc
module add mvapich2/gcc

lammpsdir="../LAMMPS-PreGenome/src/"
mpirun -np 14 $lammpsdir/lmp_openmpi -in in.chromosome
