#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=1:00:00
#SBATCH --job-name=Connected-Components-Detection
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16


module load gcc openmpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./bin/connected_components_mpi
