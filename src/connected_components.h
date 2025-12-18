/**
 * @file connected_components.h
 * @brief Connected components counting algorithm for sparse binary matrices.
 *
 * Provides a parallel implementation using OpenMP.
 */

#ifndef CONNECTED_COMPONENTS_H
#define CONNECTED_COMPONENTS_H

#include "matrix.h"

/**
 * @brief Computes connected components using a parallel union-find algorithm.
 *
 * @param matrix Sparse binary matrix in CSC format
 * @return Number of connected components, or -1 on error
 */
int connected_components(const CSCBinaryMatrix *matrix);

#endif /* CONNECTED_COMPONENTS_H */
