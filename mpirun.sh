#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:10:00
#SBATCH --job-name=Connected-Components-Detection-OpenMPI+OpenMP
#SBATCH --output=becnhmarks/rome-8-16-%j.json
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

module load gcc openmpi

mkdir -p benchmarks

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false

export RETRIES=10
export FILEPATH=$SCRATCH/data/bin/com-Friendster.bin

srun --cpu-bind=cores ./bin/connected_components_mpi -n $RETRIES $FILEPATH
