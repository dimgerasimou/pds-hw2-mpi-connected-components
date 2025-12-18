/**
 * @file matrix.c
 * @brief CSC (Compressed Sparse Column) binary matrix utilities.
 *
 * This module implements loading, storing, and printing of binary sparse
 * matrices in CSC format. It supports both matrix market (.mtx) format and
 * binary (.bin) format.
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

#ifdef _OPENMP
#include <omp.h>
#endif

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
 * @brief Get file extension from filename.
 *
 * @param filename Path to file.
 * @return Pointer to extension (without dot), or NULL if no extension.
 */
static const char*
get_file_extension(const char *filename)
{
	const char *dot = strrchr(filename, '.');
	if (!dot || dot == filename)
		return NULL;
	return dot + 1;
}

/**
 * @brief Load a CSC matrix from a binary (.bin) file with chunked I/O.
 *
 * Binary format:
 *   - uint32_t nrows
 *   - uint32_t ncols
 *   - size_t   nnz
 *   - uint32_t col_ptr[ncols + 1]
 *   - uint32_t row_idx[nnz]
 *
 * @param filename Path to the .bin file.
 * @return Newly allocated CSCBinaryMatrix on success, NULL on error.
 */
static CSCBinaryMatrix*
csc_load_binary(const char *filename)
{
	FILE *f = fopen(filename, "rb");
	if (!f) {
		print_error(__func__, "failed to open .bin file", errno);
		return NULL;
	}

	CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
	if (!m) {
		print_error(__func__, "malloc failed", errno);
		fclose(f);
		return NULL;
	}

	/* read header */
	uint32_t nrows_u32, ncols_u32;
	if (fread(&nrows_u32, sizeof(uint32_t), 1, f) != 1 ||
	    fread(&ncols_u32, sizeof(uint32_t), 1, f) != 1 ||
	    fread(&m->nnz, sizeof(size_t), 1, f) != 1)
	{
		print_error(__func__, "failed to read header", errno);
		free(m);
		fclose(f);
		return NULL;
	}

	m->nrows = nrows_u32;
	m->ncols = ncols_u32;

	/* allocate arrays */
	m->col_ptr = malloc((m->ncols + 1) * sizeof(uint32_t));
	m->row_idx = malloc(m->nnz * sizeof(uint32_t));

	if (!m->col_ptr || !m->row_idx) {
		print_error(__func__, "malloc failed", errno);
		csc_free_matrix(m);
		fclose(f);
		return NULL;
	}

	/* read col_ptr */
	if (fread(m->col_ptr, sizeof(uint32_t), m->ncols + 1, f) != m->ncols + 1) {
		print_error(__func__, "failed to read col_ptr", errno);
		csc_free_matrix(m);
		fclose(f);
		return NULL;
	}

	/* read row_idx in chunks for large matrices */
	size_t chunk_size = 100000000; /* 100M entries */
	size_t read_so_far = 0;

	while (read_so_far < m->nnz) {
		size_t to_read = (m->nnz - read_so_far < chunk_size) ?
		                 (m->nnz - read_so_far) : chunk_size;

		size_t actual = fread(&m->row_idx[read_so_far],
		                      sizeof(uint32_t), to_read, f);

		if (actual != to_read) {
			print_error(__func__, "failed to read row_idx chunk", errno);
			csc_free_matrix(m);
			fclose(f);
			return NULL;
		}

		read_so_far += actual;
	}

	fclose(f);
	return m;
}

/**
 * @brief Load a CSC matrix from a Matrix Market (.mtx) file with parallel processing.
 *
 * Supports the following formats:
 *
 * - coordinate or array
 * - pattern or real-valued
 * - general, symmetric, skew-symmetric, hermitian
 *
 * Only non-zero entries are stored (binary interpretation).
 * Uses OpenMP for parallel processing of large matrices.
 *
 * @param filename Path to the .mtx file.
 * @return Newly allocated CSCBinaryMatrix on success, NULL on error.
 */
static CSCBinaryMatrix*
csc_load_mtx(const char *filename)
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
	int is_pattern = (strcmp(field, "pattern") == 0);
	int symmetric = (strcmp(symmetry, "symmetric") == 0);
	int skew = (strcmp(symmetry, "skew-symmetric") == 0);
	int hermitian = (strcmp(symmetry, "hermitian") == 0);
	int general = (strcmp(symmetry, "general") == 0);

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

	/* --- Load all data into memory (sequential I/O) -------------------- */
	size_t *raw_i = malloc(nnz * sizeof(size_t));
	size_t *raw_j = malloc(nnz * sizeof(size_t));
	double *raw_v = is_pattern ? NULL : malloc(nnz * sizeof(double));

	if (!raw_i || !raw_j || (!is_pattern && !raw_v)) {
		print_error(__func__, "malloc failed for raw data", errno);
		fclose(f);
		free(raw_i); free(raw_j); free(raw_v);
		return NULL;
	}

	/* read entries */
	if (is_coordinate) {
		for (size_t k = 0; k < nnz; k++) {
			if (is_pattern) {
				if (fscanf(f, "%zu %zu", &raw_i[k], &raw_j[k]) != 2) {
					print_error(__func__, "bad coordinate entry", 0);
					goto fail_raw;
				}
			} else {
				if (fscanf(f, "%zu %zu %lf",
				           &raw_i[k], &raw_j[k], &raw_v[k]) != 3)
				{
					print_error(__func__, "bad coordinate entry", 0);
					goto fail_raw;
				}
			}
			/* convert to 0-indexed */
			raw_i[k]--;
			raw_j[k]--;
		}
	} else {
		/* array: dense stored column-major */
		size_t k = 0;
		for (size_t j = 0; j < ncols; j++) {
			for (size_t i = 0; i < nrows; i++) {
				if (fscanf(f, "%lf", &raw_v[k]) != 1) {
					print_error(__func__, "bad array entry", 0);
					goto fail_raw;
				}
				raw_i[k] = i;
				raw_j[k] = j;
				k++;
			}
		}
	}

	fclose(f);

	/* --- Parallel: Filter non-zeros and expand symmetric --------------- */
	size_t max_nnz = nnz * (symmetric ? 2 : 1);
	uint32_t *coo_i = malloc(max_nnz * sizeof(uint32_t));
	uint32_t *coo_j = malloc(max_nnz * sizeof(uint32_t));

	if (!coo_i || !coo_j) {
		print_error(__func__, "malloc failed for COO arrays", errno);
		goto fail_raw;
	}

#ifdef _OPENMP
	int nthreads = omp_get_max_threads();
	size_t *thread_counts = calloc(nthreads, sizeof(size_t));
	size_t *thread_offsets = calloc(nthreads, sizeof(size_t));

	if (!thread_counts || !thread_offsets) {
		print_error(__func__, "malloc failed for thread arrays", errno);
		free(coo_i); free(coo_j);
		goto fail_raw;
	}

	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		size_t local_count = 0;

		/* phase 1: count non-zeros per thread */
		#pragma omp for schedule(static) nowait
		for (size_t k = 0; k < nnz; k++) {
			double val = is_pattern ? 1.0 : raw_v[k];
			if (val != 0.0) {
				local_count++;
				if (symmetric && raw_i[k] != raw_j[k])
					local_count++;
			}
		}

		thread_counts[tid] = local_count;

		/* phase 2: compute offsets */
		#pragma omp barrier
		#pragma omp single
		{
			size_t offset = 0;
			for (int t = 0; t < nthreads; t++) {
				thread_offsets[t] = offset;
				offset += thread_counts[t];
			}
		}

		/* phase 3: write to COO arrays */
		size_t write_pos = thread_offsets[tid];
		#pragma omp for schedule(static)
		for (size_t k = 0; k < nnz; k++) {
			double val = is_pattern ? 1.0 : raw_v[k];
			if (val != 0.0) {
				coo_i[write_pos] = raw_i[k];
				coo_j[write_pos] = raw_j[k];
				write_pos++;

				if (symmetric && raw_i[k] != raw_j[k]) {
					coo_i[write_pos] = raw_j[k];
					coo_j[write_pos] = raw_i[k];
					write_pos++;
				}
			}
		}
	}

	size_t count = thread_offsets[nthreads - 1] + thread_counts[nthreads - 1];
	free(thread_counts);
	free(thread_offsets);
#else
	/* serial fallback */
	size_t count = 0;
	for (size_t k = 0; k < nnz; k++) {
		double val = is_pattern ? 1.0 : raw_v[k];
		if (val != 0.0) {
			coo_i[count] = raw_i[k];
			coo_j[count] = raw_j[k];
			count++;

			if (symmetric && raw_i[k] != raw_j[k]) {
				coo_i[count] = raw_j[k];
				coo_j[count] = raw_i[k];
				count++;
			}
		}
	}
#endif

	free(raw_i); free(raw_j); free(raw_v);

	/* --- Convert COO → CSC binary -------------------------------------- */
	CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
	if (!m) {
		print_error(__func__, "malloc failed for matrix", errno);
		goto fail_coo;
	}

	m->nrows = nrows;
	m->ncols = ncols;
	m->nnz = count;

	m->row_idx = malloc(count * sizeof(uint32_t));
	m->col_ptr = malloc((ncols + 1) * sizeof(uint32_t));

	if (!m->row_idx || !m->col_ptr) {
		print_error(__func__, "malloc failed for CSC arrays", errno);
		csc_free_matrix(m);
		goto fail_coo;
	}

	/* parallel: count entries per column */
	memset(m->col_ptr, 0, (ncols + 1) * sizeof(uint32_t));

#ifdef _OPENMP
	#pragma omp parallel
	{
		/* thread-local column counts */
		uint32_t *local_col_count = calloc(ncols + 1, sizeof(uint32_t));
		if (!local_col_count) {
			print_error(__func__, "malloc failed in parallel", errno);
		} else {
			#pragma omp for schedule(static) nowait
			for (size_t k = 0; k < count; k++)
				local_col_count[coo_j[k] + 1]++;

			/* reduce into global col_ptr */
			#pragma omp critical
			{
				for (size_t j = 0; j <= ncols; j++)
					m->col_ptr[j] += local_col_count[j];
			}

			free(local_col_count);
		}
	}
#else
	/* serial fallback */
	for (size_t k = 0; k < count; k++)
		m->col_ptr[coo_j[k] + 1]++;
#endif

	/* prefix sum */
	for (size_t j = 0; j < ncols; j++)
		m->col_ptr[j + 1] += m->col_ptr[j];

	/* fill row indices */
	uint32_t *col_fill = calloc(ncols, sizeof(uint32_t));
	if (!col_fill) {
		print_error(__func__, "malloc failed for col_fill", errno);
		csc_free_matrix(m);
		goto fail_coo;
	}

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

fail_coo:
	free(coo_i);
	free(coo_j);
	return NULL;

fail_raw:
	fclose(f);
	free(raw_i);
	free(raw_j);
	free(raw_v);
	return NULL;
}

/* ------------------------------------------------------------------------- */
/*                            Public API Functions                           */
/* ------------------------------------------------------------------------- */

/**
 * @brief Load a CSC matrix from file (auto-detects format by extension).
 *
 * Supports:
 * - .bin files: fast binary format
 * - .mtx files: Matrix Market format (uses parallel loading if OpenMP enabled)
 *
 * @param filename Path to the matrix file.
 * @return Newly allocated CSCBinaryMatrix on success, NULL on error.
 */
CSCBinaryMatrix*
csc_load_matrix(const char *filename)
{
	const char *ext = get_file_extension(filename);
	if (!ext) {
		print_error(__func__, "no file extension found", 0);
		return NULL;
	}

	if (strcmp(ext, "bin") == 0)
		return csc_load_binary(filename);
	else if (strcmp(ext, "mtx") == 0)
		return csc_load_mtx(filename);
	else {
		print_error(__func__, "unsupported file extension", 0);
		return NULL;
	}
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

	if (m->row_idx) {
		free(m->row_idx);
		m->row_idx = NULL;
	}

	if (m->col_ptr) {
		free(m->col_ptr);
		m->col_ptr = NULL;
	}

	free(m);
	m = NULL;
}
