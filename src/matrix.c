/**
 * @file matrix.c
 * @brief CSC (Compressed Sparse Column) binary matrix utilities.
 *
 * This module implements loading, storing, and printing of binary sparse
 * matrices in CSC format. It supports only matrix market (.mtx) format.
 *
 * Only binary matrices are represented. Any non-zero numeric values in
 * the input are treated as 1.
 */
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "error.h"

/* ------------------------------------------------------------------------- */
/*                            Static Helper Functions                        */
/* ------------------------------------------------------------------------- */

/**
 * @brief Skip comments and blank lines in a Matrix Market file.
 *
 * Advances the file cursor until the next non-comment, non-empty line.
 *
 * @param f Open file pointer.
 * @return Always 0.
 */
static int
mm_skip_comments(FILE *f)
{
	long pos;
	int c;

	while (1) {
		pos = ftell(f);
		c = fgetc(f);
		if (c == '%') {          /* comment line */
			while ((c = fgetc(f)) != '\n' && c != EOF);
			continue;
		} else if (isspace(c)) { /* blank line */
			continue;
		}
		fseek(f, pos, SEEK_SET);
		return 0;
	}
}

/**
 * @brief Load a CSC matrix from a Matrix Market (.mtx) file.
 *
 * Supports the following formats:
 *
 * - coordinate or array
 * - pattern or real-valued
 * - general, symmetric, skew-symmetric, hermitian
 *
 * Only non-zero entries are stored (binary interpretation).
 *
 * @param filename Path to the .mtx file.
 * @return Newly allocated CSCBinaryMatrix on success, NULL on error.
 */
CSCBinaryMatrix*
csc_load_matrix(const char *filename)
{
	FILE *f = fopen(filename, "r");
	if (!f) {
		print_error(__func__, "failed to open .mtx file", errno);
		return NULL;
	}

	/* --- Header -------------------------------------------------------- */
	char format[64], field[64], symmetry[64];

	if (fscanf(f, "%%%%MatrixMarket matrix %63s %63s %63s",
				format, field, symmetry) != 3)
	{
		print_error(__func__, "invalid MatrixMarket header", 0);
		fclose(f);
		return NULL;
	}

	int is_coordinate = (strcmp(format, "coordinate") == 0);
	//int is_array      = (strcmp(format, "array") == 0);

	int is_pattern = (strcmp(field, "pattern") == 0);

	int symmetric   = (strcmp(symmetry, "symmetric") == 0);
	int skew        = (strcmp(symmetry, "skew-symmetric") == 0);
	int hermitian   = (strcmp(symmetry, "hermitian") == 0);
	int general     = (strcmp(symmetry, "general") == 0);

	if (!general && !symmetric && !skew && !hermitian) {
		print_error(__func__, "unsupported symmetry", 0);
		fclose(f);
		return NULL;
	}

	/* --- Sizes --------------------------------------------------------- */
	mm_skip_comments(f);

	size_t nrows, ncols, nnz;
	if (is_coordinate) {
		if (fscanf(f, "%zu %zu %zu", &nrows, &ncols, &nnz) != 3) {
			print_error(__func__, "invalid size line", 0);
			fclose(f);
			return NULL;
		}
	} else { /* array */
		if (fscanf(f, "%zu %zu", &nrows, &ncols) != 2) {
			print_error(__func__, "invalid array size line", 0);
			fclose(f);
			return NULL;
		}
		nnz = nrows * ncols; /* will filter zeroes later */
	}

	/* Allocate temporary COO arrays */
	size_t max_nnz = nnz * (symmetric ? 2 : 1) + 5;
	uint32_t *coo_i = malloc(max_nnz * sizeof(uint32_t));
	uint32_t *coo_j = malloc(max_nnz * sizeof(uint32_t));
	if (!coo_i || !coo_j) {
		print_error(__func__, "malloc failed", errno);
		fclose(f);
		free(coo_i); free(coo_j);
		return NULL;
	}

	size_t count = 0;

	/* --- Read entries -------------------------------------------------- */
	if (is_coordinate) {
		/* i j [value] */
		for (size_t k = 0; k < nnz; k++) {
			size_t i, j;
			double val = 1.0;

			if (is_pattern) {
				if (fscanf(f, "%zu %zu", &i, &j) != 2) {
					print_error(__func__, "bad coordinate entry", 0);
					goto fail;
				}
			} else {
				if (fscanf(f, "%zu %zu %lf", &i, &j, &val) != 3) {
					print_error(__func__, "bad coordinate entry", 0);
					goto fail;
				}
			}

			if (val != 0.0) {
				coo_i[count] = i - 1;
				coo_j[count] = j - 1;
				count++;

				if (symmetric && i != j) {
					coo_i[count] = j - 1;
					coo_j[count] = i - 1;
					count++;
				}
			}
		}
	}
	else {
		/* array: dense stored column-major */
		for (size_t j = 0; j < ncols; j++) {
			for (size_t i = 0; i < nrows; i++) {
				double val;
				if (fscanf(f, "%lf", &val) != 1)
					goto fail;

				if (val != 0.0) {
					coo_i[count] = i;
					coo_j[count] = j;
					count++;
				}
			}
		}
	}

	fclose(f);

	/* --- Convert COO → CSC binary -------------------------------------- */
	CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
	if (!m) goto fail;

	m->nrows = nrows;
	m->ncols = ncols;
	m->nnz   = count;

	m->row_idx = malloc(count * sizeof(uint32_t));
	m->col_ptr = malloc((ncols + 1) * sizeof(uint32_t));

	if (!m->row_idx || !m->col_ptr) {
		print_error(__func__, "malloc failed", errno);
		csc_free_matrix(m);
		goto fail;
	}

	/* count entries per column */
	memset(m->col_ptr, 0, (ncols + 1) * sizeof(uint32_t));

	for (size_t k = 0; k < count; k++)
		m->col_ptr[coo_j[k] + 1]++;

	for (size_t j = 0; j < ncols; j++)
		m->col_ptr[j+1] += m->col_ptr[j];

	/* fill rows */
	uint32_t *col_fill = calloc(ncols, sizeof(uint32_t));
	if (!col_fill) goto fail2;

	for (size_t k = 0; k < count; k++) {
		uint32_t j = coo_j[k];
		uint32_t dest = m->col_ptr[j] + col_fill[j];
		m->row_idx[dest] = coo_i[k];
		col_fill[j]++;
	}

	free(col_fill);
	free(coo_i);
	free(coo_j);

	return m;

fail2:
	csc_free_matrix(m);

fail:
	fclose(f);
	free(coo_i);
	free(coo_j);
	return NULL;
}

/**
 * @brief Free a CSCBinaryMatrix and its associated memory.
 *
 * Safe to call with NULL.
 *
 * @param m CSC matrix to free.
 */
void
csc_free_matrix(CSCBinaryMatrix *m)
{
	if (!m)
		return;

	if(m->row_idx){
		free(m->row_idx);
		m->row_idx = NULL;
	}

	if(m->col_ptr){
		free(m->col_ptr);
		m->col_ptr = NULL;
	}

	free(m);
	m = NULL;
}

