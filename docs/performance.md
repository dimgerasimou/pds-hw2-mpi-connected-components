# Performance

This report summarizes the benchmark runs produced by the current implementation.

## How to read the output filenames

Each benchmark output uses the format:

`rome-nnodes-ntasks_per_node-ncores_per_task-matrix_name-job_no.out`

Where:
- `nnodes`: number of nodes
- `ntasks_per_node`: MPI ranks per node
- `ncores_per_task`: OpenMP threads per rank
- `matrix_name`: input graph/matrix
- `job_no`: Slurm job id

## Benchmark environment

All benchmarking has been done on the “Aristotelis” HPC cluster, of the Aristotle University of Thessaloniki.

- **System tag:** `rome`
- **CPU:** AMD EPYC 7662 64-Core Processor
- **Memory (per node):** 256 GiB
- **Run date:** 2025-12-19

## Datasets

| Name                  | Nodes       | Edges         |
|----------------------:|------------:|--------------:|
| **com-Friendster**    | 65,608,366  | 3,612,134,270 |
| **mawi_201512020330** | 226,196,185 | 480,047,894   |

---

# com-Friendster

## 1 node: OpenMP-only scaling (1 rank/node)

| Threads per rank | Total cores | Mean time (s) | Throughput (Medges/s) |
|-----------------:|------------:|--------------:|----------------------:|
|                8 |           8 |        96.116 |                  37.6 |
|               16 |          16 |        60.748 |                  59.5 |
|               32 |          32 |        54.324 |                  66.5 |
|               64 |          64 |       104.843 |                  34.5 |
|              128 |         128 |       100.488 |                  35.9 |

## 1 node @ 128 cores: effect of rank splitting

| Ranks per node | Threads per rank | Mean time (s) | Throughput (Medges/s) | Load time (s) |
|---------------:|-----------------:|--------------:|----------------------:|--------------:|
|              1 |              128 |       100.488 |                  35.9 |         5.305 |
|              2 |               64 |        29.354 |                 123.1 |         1.662 |
|              4 |               32 |        28.766 |                 125.6 |         0.896 |
|              8 |               16 |        29.27  |                 123.4 |         0.55  |
|             16 |                8 |        29.204 |                 123.7 |         0.466 |

## Strong scaling with the recommended per-node configuration

| Number of nodes | Total Cores | Mean time (s) | Throughput (Medges/s) | Speedup vs 1 node (x) | Efficency vs 1 node (x) |
|----------------:|------------:|--------------:|----------------------:|----------------------:|------------------------:|
|               1 |         128 |        29.354 |                 123.1 |                  1    |                   1     |
|               2 |         256 |        16.438 |                 219.7 |                  1.79 |                   0.893 |
|               4 |         512 |         9.483 |                 380.9 |                  3.1  |                   0.774 |
|               8 |        1024 |         6.397 |                 564.7 |                  4.59 |                   0.574 |

<details>
<summary>Full friendster results</summary>

| Number of nodes | Ranks per node | Threads per rank | Total cores | Mean time (s) | Throughput (Medges/s) | Load time (s) |
|----------------:|---------------:|-----------------:|------------:|--------------:|----------------------:|--------------:|
|               1 |              1 |                8 |           8 |        96.116 |                  37.6 |         2.687 |
|               1 |              1 |               16 |          16 |        60.748 |                  59.5 |         2.394 |
|               1 |              1 |               32 |          32 |        54.324 |                  66.5 |         2.396 |
|               1 |              1 |               64 |          64 |       104.843 |                  34.5 |         4.829 |
|               1 |              1 |              128 |         128 |       100.488 |                  35.9 |         5.305 |
|               1 |              2 |               64 |         128 |        29.354 |                 123.1 |         1.662 |
|               1 |              4 |               32 |         128 |        28.766 |                 125.6 |         0.896 |
|               1 |              8 |               16 |         128 |        29.27  |                 123.4 |         0.55  |
|               1 |             16 |                8 |         128 |        29.204 |                 123.7 |         0.466 |
|               2 |              1 |               64 |         128 |        37.865 |                  95.4 |         1.525 |
|               2 |              1 |              128 |         256 |        62.946 |                  57.4 |         3.152 |
|               2 |              2 |               64 |         256 |        16.438 |                 219.7 |         2.492 |
|               2 |              4 |               32 |         256 |        15.848 |                 227.9 |         0.472 |
|               2 |              8 |               16 |         256 |        15.957 |                 226.4 |         2.47  |
|               4 |              1 |               64 |         256 |        24.593 |                 146.9 |         0.815 |
|               4 |              2 |               64 |         512 |         9.483 |                 380.9 |         2.406 |
|               8 |              2 |               64 |        1024 |         6.397 |                 564.7 |         0.212 |

</details>

---

# mawi_201512020330

## 1 node: OpenMP-only scaling (1 rank/node)

| Threads per rank | Total Cores | Mean time (s) | Throughput (Medges/s) |
|-----------------:|------------:|--------------:|----------------------:|
|                8 |           8 |         7.106 |                  67.6 |
|               16 |          16 |         6.299 |                  76.2 |
|               32 |          32 |         6.139 |                  78.2 |
|               64 |          64 |         6.089 |                  78.8 |
|              128 |         128 |         8.277 |                  58   |

## 1 node @ 128 cores: effect of rank splitting

| Ranks per node | Threads per rank | Mean time (s) | Throughput (Medges/s) | Load time (s) |
|---------------:|-----------------:|--------------:|----------------------:|--------------:|
|              1 |              128 |         8.277 |                  58   |         1.208 |
|              2 |               64 |         5.832 |                  82.3 |         0.564 |
|              4 |               32 |         6.895 |                  69.6 |         0.426 |
|              8 |               16 |         7.526 |                  63.8 |         0.094 |
|             16 |                8 |         8.44  |                  56.9 |         0.068 |

## Strong scaling with the recommended per-node configuration

| Number of nodes | Total Cores | Mean time (s) | Throughput (Medges/s) | Speedup vs 1 node (x) | Efficency vs 1 node (x) |
|----------------:|------------:|--------------:|----------------------:|----------------------:|------------------------:|
|               1 |         128 |         5.832 |                  82.3 |                  1    |                   1     |
|               2 |         256 |         9.661 |                  49.7 |                  0.6  |                   0.302 |
|               4 |         512 |        10.235 |                  46.9 |                  0.57 |                   0.142 |
|               8 |        1024 |        10.542 |                  45.5 |                  0.55 |                   0.069 |

<details>
<summary>Full mawi results</summary>

| Number of nodes | Ranks per node | Threads per rank | Total cores | Mean time (s) | Throughput (Medges/s) | Load time (s) |
|----------------:|---------------:|-----------------:|------------:|--------------:|----------------------:|--------------:|
|               1 |              1 |                8 |           8 |         7.106 |                  67.6 |         0.853 |
|               1 |              1 |               16 |          16 |         6.299 |                  76.2 |         0.736 |
|               1 |              1 |               32 |          32 |         6.139 |                  78.2 |         0.714 |
|               1 |              1 |               64 |          64 |         6.089 |                  78.8 |         0.739 |
|               1 |              1 |              128 |         128 |         8.277 |                  58   |         1.208 |
|               1 |              2 |               64 |         128 |         5.832 |                  82.3 |         0.564 |
|               1 |              4 |               32 |         128 |         6.895 |                  69.6 |         0.426 |
|               1 |              8 |               16 |         128 |         7.526 |                  63.8 |         0.094 |
|               1 |             16 |                8 |         128 |         8.44  |                  56.9 |         0.068 |
|               2 |              1 |               64 |         128 |         8.244 |                  58.2 |         0.521 |
|               2 |              2 |               64 |         256 |         9.661 |                  49.7 |         0.145 |
|               2 |              4 |               32 |         256 |        10.076 |                  47.6 |         0.296 |
|               2 |              8 |               16 |         256 |        10.352 |                  46.4 |         0.044 |
|               4 |              1 |               64 |         256 |        10.697 |                  44.9 |         0.17  |
|               4 |              2 |               64 |         512 |        10.235 |                  46.9 |         0.199 |
|               8 |              2 |               64 |        1024 |        10.542 |                  45.5 |         0.166 |

</details>

---

## Observations

- On **friendster**, a single MPI rank per node with many threads performs poorly at high thread counts. Splitting the node into **multiple MPI ranks** is dramatically better.
- On **mawi**, the best result is still on **one node**; adding nodes increases overhead and makes runtime worse for this this implementation.
