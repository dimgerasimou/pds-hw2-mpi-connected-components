/**
 * @file matrix.h
 * @brief CSC (Compressed Sparse Column) binary matrix structure and API.
 *
 * Provides functionality to load, free, and print sparse binary matrices
 * stored in CSC format. Non-zero values are represented implicitly as 1.
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
	size_t nrows;       /**< Number of rows in the matrix */
	size_t ncols;       /**< Number of columns in the matrix */
	size_t nnz;         /**< Number of non-zero (1) entries */
	uint32_t *row_idx;  /**< Row indices of non-zero elements (length nnz) */
	uint32_t *col_ptr;  /**< Column pointers (length ncols + 1) */
} CSCBinaryMatrix;

/** @brief Load a sparse binary matrix from a .mat or .mtx file.
 *
 * Dispatches automatically based on file extension.
 *
 * @param path Path to the matrix file.
 * @return Newly allocated CSCBinaryMatrix, or NULL on failure.
 *
 * @note The returned matrix must be freed using csc_free_matrix().
 */
CSCBinaryMatrix *csc_load_matrix(const char *path);

/**
 * @brief Free a CSCBinaryMatrix and its associated memory.
 *
 * Safe to call with NULL.
 *
 * @param m CSC matrix to free.
 */
void csc_free_matrix(CSCBinaryMatrix *m);

/**
 * @brief Print a sparse binary matrix in coordinate format.
 *
 * @param m Pointer to the CSCBinaryMatrix to print.
 *
 * @note Prints as (row, col) pairs. Indices are 1-based.
 */
void csc_print_matrix(CSCBinaryMatrix *m);

#endif /* MATRIX_H */
