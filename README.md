<h1 align="center">Connected Components: Distributed-Memory Implementation</h1>
<h3 align="center">Parallel and Distributed Systems</h3>

<p align="center">
  <em>Department of Electrical and Computer Engineering</em><br>
  <em>Aristotle University of Thessaloniki</em><br>
  <strong>Homework #2</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/language-C-blue.svg" alt="Language">
  <img src="https://img.shields.io/badge/parallelism-MPI%20%2B%20OpenMP-success" alt="MPI + OpenMP">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20HPC-lightgrey" alt="Platform">
</p>

## Overview

This project implements a **high-performance distributed-memory connected components algorithm** for **massive sparse graphs**, targeting modern HPC systems.

The focus is on **scalability beyond a single node**, combining:
- **MPI** for inter-node communication and data distribution
- **OpenMP** for intra-node shared-memory parallelism

The implementation is designed to process graphs with **billions of edges**, using a **custom binary CSC format** to minimize I/O overhead and memory footprint.

## Features

- **Distributed-memory connected components**
  - Graph partitioned by column blocks across MPI ranks
  - Each rank loads *only its local partition*
- **Hybrid parallelism**
  - MPI between nodes
  - OpenMP within each node
- **Binary CSC graph format**
  - Fast loading
  - Compact storage
  - Avoids Matrix Market parsing overhead at scale
- **Scales to billions of edges**
  - Tested on real-world graphs (Friendster, MAWI)
- **Automated benchmarking**
  - Multiple trials
  - Mean, median, std dev, min/max
  - Throughput in edges/sec
- **NUMA- and cluster-aware**
  - Explicit control over ranks, threads, and memory distribution

## Algorithm

The project implements a **parallel label propagation–style connected components algorithm** adapted for distributed memory:

1. Each MPI rank owns a contiguous block of columns (CSC format)
2. Labels are initialized locally
3. Iterative propagation:
   - Local relaxation with OpenMP
   - Boundary updates exchanged via MPI
4. Convergence detection via global reductions
5. Final component count computed collectively

This approach minimizes memory usage per rank and limits communication to boundary data.

## Dependencies


### Required

| Dependency | Purpose |
|----------|---------|
| `MPI` (OpenMPI / MPICH / system MPI) | Distributed execution |
| `OpenMP` | Intra-node parallelism |
| `gcc` or `clang` | Compilation |
| `make` | Build system |

### Environment

The code is intended to run on:
- Linux-based HPC clusters
- SLURM-managed systems

## Build

Compile main target:

```bash
make
```

This produces: `bin/connected_components_mpi`

Clean build:

```bash
make clean
make
```

Use `make converter` to compile the matrix market to binary converted needed.


## Input Format

### Binary CSC Format

Graphs are stored in a custom binary **Compressed Sparse Column CSC** format:

```c
uint32_t nrows
uint32_t ncols
size_t   nnz
uint32_t col_ptr[ncols + 1]
uint32_t row_idx[nnz]
```

### Convert from Matrix Market

```bash
./bin/mtx_to_bin input.mtx output.bin
```

This conversion step is **mandatory** before execution.

## Usage

First you need to convert the matrix market array to a binary file for the program to read. This can be done with:
```bash
make converter
./bin/mtx_to_bin input.mtx output.bin
```

Then edit `run_slurm.sh` to match your configuration and run it with:
```bash
sbatch run_slurm.sh
```

## Output

For each run, the program reports:
- Number of connected components
- Execution time statistics
- Throughput in edges per second

Results are printed in a machine-readable .json format suitable for later analysis.

## Performance Results
**Example: Friendster Social Network**
 - **Nodes**: 65,608,366
 - **Edges**: 3,612,134,270
 - **System**: AMD EPYC 7662
 - **Memory**: ~256 GB per node

| MPI Ranks	| Threads / Rank |	Total Cores	| Mean Time (s)	| Throughput (M edges/s) |
|----------:|---------------:|-------------:|--------------:|-----------------------:|
| 1         | 16             | 16           | 11.10         | 325                    |
| 1         | 64             | 64           | 4.84          | 736                    |
| 2         | 64             | 128          | 29.35         | 123                    |
| 4         | 64             | 256          | 16.44         | 220                    |
| 16        | 64             | 1024         | 6.40          | 565                    |

Results show:
- Communication overhead dominating at higher rank counts.
- Stable convergence behavior across configurations.

## Notes on Scalability
- The implementation is **memory-bandwidth bound** on a single node.
- At scale, **MPI communication and synchronization dominate**.

- Best performance achieved with:
  - Few ranks per node
  - Many OpenMP threads per rank

## Troubleshooting

**Out of Memory**
- Ensure the binary format is used.
- Reduce ranks per node.
- Increase `OMP_STACKSIZE`.

**Poor Scaling**
- Avoid oversubscribing cores.
- Prefer `--distribution=block:block`.
- Disable OpenMP dynamic teams:
```bash
export OMP_DYNAMIC=false
```

---

<p align="center"> <sub>December 2025 • Aristotle University of Thessaloniki</sub> </p>
