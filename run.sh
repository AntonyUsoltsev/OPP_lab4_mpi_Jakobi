#!/bin/bash

#PBS -l select=1:ncpus=1:mpiprocs=1:mem=10000m,place=scatter
#PBS -l walltime=00:01:00
#PBS -m n
#PBS -o out.txt
#PBS -e err.txt

MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

cd $PBS_O_WORKDIR

echo "Running $MPI_NP MPI processes"

mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./main