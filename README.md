# Connected Components — Distributed-Memory (MPI) Implementation

<p align="center">
  <img src="https://img.shields.io/badge/parallelism-MPI%20%2B%20OpenMP-blueviolet" alt="MPI + OpenMP">
</p>

Assignment #2 for the Parallel and Distributed Systems course.

A **distributed-memory connected components** implementation for massive sparse graphs, designed for modern **HPC systems** using **MPI + OpenMP**.

The focus is on scalability beyond a single node, efficient memory usage, and realistic benchmarking on large real-world graphs.

## Overview

This project implements a high-performance connected components algorithm for distributed-memory systems.

Key goals:
- Scale to graphs with billions of edges
- Minimize memory usage per MPI rank
- Combine MPI-based distribution with OpenMP intra-node parallelism

## Features

- Hybrid **MPI + OpenMP** parallelism
- Distributed-memory connected components algorithm
- Column-block graph partitioning across MPI ranks
- Custom binary **CSC graph format** for fast I/O
- Automated benchmarking and JSON output
- Designed for SLURM-managed HPC clusters

## Algorithm

The implementation follows a **label-propagation-style** connected components algorithm adapted for distributed memory:

1. Each MPI rank owns a contiguous block of columns (CSC format)
2. Labels are initialized locally
3. Iterative propagation:
   - Local relaxation using OpenMP
   - Boundary updates exchanged via MPI
4. Convergence detected using collective reductions
5. Final component count computed collectively

This approach limits communication to boundary data and minimizes per-rank memory usage.

## Build

### Requirements
- MPI implementation (OpenMPI / MPICH / system MPI)
- OpenMP
- C compiler (`gcc` or `clang`)
- `make`

### Compile
```bash
make
```

Produces:
```
bin/connected_components_mpi
```

## Input Format

Graphs are stored in a custom binary **Compressed Sparse Column (CSC)** format optimized for fast loading and low memory overhead.

Matrix Market files must be converted before execution.

```bash
make converter
./bin/mtx_to_bin input.mtx output.bin
```

## Usage

Execution is intended for HPC environments using SLURM.

```bash
sbatch run_slurm.sh
```

The executable is launched via `srun` inside the SLURM job script.  
The SLURM script controls:
- Number of nodes
- MPI ranks per node
- OpenMP threads per rank

## Performance Results

Detailed benchmark tables and scaling results are available in [docs/performance.md](docs/performance.md).

## Notes on Scalability

- Single-node performance is primarily memory-bandwidth bound
- At scale, MPI communication and synchronization dominate
- Best results achieved with:
  - Few MPI ranks per node
  - Many OpenMP threads per rank
- Optimal configurations are hardware- and workload-dependent

---

<p align="center"><sub>December 2025 • Aristotle University of Thessaloniki</sub></p>
