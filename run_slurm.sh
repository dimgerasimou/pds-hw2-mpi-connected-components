#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:30:00
#SBATCH --job-name=Connected-Components-Detection-OpenMPI+OpenMP
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --mem-bind=local
#SBATCH --output=benchmarks/rome-4-16-%j.out
#SBATCH --error=err/slurm-%j.err

set -euo pipefail

module load gcc/9.2.0 openmpi/3.1.4

make clean 1>/dev/null
make -j "${SLURM_CPUS_PER_TASK}" 1>/dev/null

# OpenMP runtime
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OMP_STACKSIZE=256M

# App args
export RETRIES=10
export MAWI="$SCRATCH/data/bin/mawi_201512020330.bin"
export FRIENDSTER="${SCRATCH}/data/bin/com-Friendster.bin"

srun --mpi=pmix --cpu-bind=cores ./bin/connected_components_mpi \
  -n "${RETRIES}" "${FRIENDSTER}"

