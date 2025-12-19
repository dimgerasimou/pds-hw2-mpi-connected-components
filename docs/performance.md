# Performance Evaluation

This document summarizes performance results for the distributed-memory connected components implementation.

---

## Experimental Setup

- **Dataset:** Friendster social network
- **Nodes:** 65,608,366
- **Edges:** 3,612,134,270
- **System:** AMD EPYC 7662
- **Memory:** ~256 GB per node

---

## Scaling Results

| MPI Ranks | Threads / Rank | Total Cores | Mean Time (s) | Throughput (M edges/s) |
|-----------|----------------|-------------|---------------|------------------------|
| 1         | 16             | 16          | 11.10         | 325                    |
| 1         | 64             | 64          | 4.84          | 736                    |
| 2         | 64             | 128         | 29.35         | 123                    |
| 4         | 64             | 256         | 16.44         | 220                    |
| 16        | 64             | 1024        | 6.40          | 565                    |

---

## Observations

- Communication overhead increases with MPI rank count
- Convergence behavior remains stable across configurations
- Hybrid parallelism mitigates intra-node overhead
