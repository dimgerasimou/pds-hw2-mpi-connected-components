#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:30:00
#SBATCH --job-name=Connected-Components-Detection-OpenMPI+OpenMP
#SBATCH --output=benchmarks/rome-4-16-%j.json
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --mem-bind=local

module load gcc openmpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OMP_STACKSIZE=256M

export RETRIES=10
export FILEPATH=$SCRATCH/data/bin/com-Friendster.bin

srun --cpu-bind=cores ./bin/connected_components_mpi -n $RETRIES $FILEPATH
