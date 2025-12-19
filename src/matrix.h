/**
 * @file matrix.h
 * @brief CSC (Compressed Sparse Column) binary matrix structure and API.
 *
 * Provides functionality to load and free sparse binary matrices stored in
 * CSC format. Non-zero values are represented implicitly as 1.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include <stdint.h>

/**
 * @struct CSCBinaryMatrix
 * @brief Compressed Sparse Column (CSC) representation of a binary matrix.
 *
 * Non-zero entries are implicitly 1. Stores only row indices and column pointers.
 */
typedef struct {
	size_t nrows;        /**< Number of rows in the matrix */
	size_t ncols;        /**< Number of columns in the local MPI matrix */
	size_t ncols_global; /**< Number of columns in the matrix */
	size_t nnz;          /**< Number of non-zero (1) entries */
	uint32_t *row_idx;   /**< Row indices of non-zero elements (length nnz) */
	uint32_t *col_ptr;   /**< Column pointers (length ncols + 1) */
} CSCBinaryMatrix;

/**
 * @brief Load a distributed partition of a sparse binary matrix (.bin only).
 *
 * Each MPI rank loads only its partition of columns. This significantly reduces
 * memory usage per node and speeds up loading.
 *
 * Partition strategy: Block distribution of columns across ranks.
 *
 * @param path Path to the matrix file (.bin format required).
 * @param mpi_rank Current MPI rank.
 * @param mpi_size Total number of MPI ranks.
 * @return Newly allocated CSCBinaryMatrix (local partition), or NULL on failure.
 *
 * @note The returned matrix must be freed using csc_free_matrix().
 */
CSCBinaryMatrix *csc_load_matrix(const char *path, int mpi_rank, int mpi_size);

/**
 * @brief Free a CSCBinaryMatrix and its associated memory.
 *
 * Safe to call with NULL.
 *
 * @param m CSC matrix to free.
 */
void csc_free_matrix(CSCBinaryMatrix *m);

#endif /* MATRIX_H */
