#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:10:00
#SBATCH --job-name=Connected-Components-Detection
#SBATCH --output=becnhmarks/rome-8-16-%j.json
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

module load gcc openmpi

mkdir -p benchmarks

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export RETRIES=10
export FILEPATH=$SCRATCH/data/bin/com-Friendster.bin

srun ./bin/connected_components_mpi -n $RETRIES $FILEPATH
