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
 * @brief Computes connected components using distributed parallel union-find.
 *
 * Multi-process version using MPI + OpenMP hybrid parallelization.
 *
 * @param[in] matrix   Local CSC partition (this rank's columns only).
 * @param[in] mpi_rank Current MPI rank.
 * @param[in] mpi_size Total number of MPI ranks.
 *
 * @return Number of connected components, or -1 on error.
 */
int connected_components(const CSCBinaryMatrix *matrix, int mpi_rank, int mpi_size);

#endif /* CONNECTED_COMPONENTS_H */
