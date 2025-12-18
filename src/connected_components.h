/**
 * @file connected_components.h
 * @brief Connected components counting algorithm for sparse binary matrices.
 *
 * Provides both single-process (OpenMP) and distributed (MPI+OpenMP)
 * implementations.
 */

#ifndef CONNECTED_COMPONENTS_H
#define CONNECTED_COMPONENTS_H

#include "matrix.h"

/**
 * @brief Computes connected components using a parallel union-find algorithm.
 *
 * Single-process version using only OpenMP.
 *
 * @param matrix Sparse binary matrix in CSC format
 * @return Number of connected components, or -1 on error
 */
int connected_components(const CSCBinaryMatrix *matrix);

/**
 * @brief Computes connected components using distributed parallel union-find.
 *
 * Multi-process version using MPI + OpenMP hybrid parallelization.
 *
 * @param matrix Local CSC partition (this rank's columns only)
 * @param mpi_rank Current MPI rank
 * @param mpi_size Total number of MPI ranks
 * @return Number of connected components, or -1 on error
 */
int connected_components_mpi(const CSCBinaryMatrix *matrix, int mpi_rank, int mpi_size);

#endif /* CONNECTED_COMPONENTS_H */
