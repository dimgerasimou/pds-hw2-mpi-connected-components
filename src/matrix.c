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
 *
 * MPI extension: Supports distributed loading where each rank loads only
 * its partition of columns to minimize memory usage.
 */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "matrix.h"

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

static CSCBinaryMatrix*
csc_load_bin(const char *filename)
{
	/* --- Open file with O_DIRECT hint for better performance ------- */
	int fd = open(filename, O_RDONLY);
	if (fd < 0) {
		print_error(__func__, "failed to open .bin file", errno);
		return NULL;
	}

	CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
	if (!m) {
		print_error(__func__, "malloc failed", errno);
		close(fd);
		return NULL;
	}

	off_t offset = 0;

	/* --- Read header (sequential, small) --------------------------- */
	uint32_t nrows_u32, ncols_u32;
	
	if (read(fd, &nrows_u32, sizeof(uint32_t)) != sizeof(uint32_t) ||
	    read(fd, &ncols_u32, sizeof(uint32_t)) != sizeof(uint32_t) ||
	    read(fd, &m->nnz, sizeof(size_t)) != sizeof(size_t))
	{
		print_error(__func__, "failed to read header", errno);
		free(m);
		close(fd);
		return NULL;
	}

	m->nrows = nrows_u32;
	m->ncols = ncols_u32;
	offset = sizeof(uint32_t) * 2 + sizeof(size_t);

	/* --- Read col_ptr (sequential, relatively small) --------------- */
	m->col_ptr = malloc((m->ncols + 1) * sizeof(uint32_t));
	if (!m->col_ptr) {
		print_error(__func__, "malloc failed for col_ptr", errno);
		free(m);
		close(fd);
		return NULL;
	}

	size_t col_ptr_bytes = (m->ncols + 1) * sizeof(uint32_t);
	if (read(fd, m->col_ptr, col_ptr_bytes) != (ssize_t)col_ptr_bytes) {
		print_error(__func__, "failed to read col_ptr", errno);
		csc_free_matrix(m);
		close(fd);
		return NULL;
	}
	offset += col_ptr_bytes;

	/* --- Allocate row_idx ------------------------------------------ */
	m->row_idx = malloc(m->nnz * sizeof(uint32_t));
	if (!m->row_idx) {
		print_error(__func__, "malloc failed for row_idx", errno);
		csc_free_matrix(m);
		close(fd);
		return NULL;
	}

	/* --- Parallel read of row_idx using pread ---------------------- */
	int n_threads = omp_get_max_threads();
	size_t chunk_size = (m->nnz + n_threads - 1) / n_threads;
	
	int read_error = 0;

	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		size_t start = tid * chunk_size;
		size_t end = start + chunk_size;
		if (end > m->nnz)
			end = m->nnz;

		if (start < end) {
			size_t count = end - start;
			off_t file_offset = offset + start * sizeof(uint32_t);
			size_t bytes_to_read = count * sizeof(uint32_t);

			/* pread is thread-safe: each thread reads from different offset */
			ssize_t bytes_read = pread(fd, &m->row_idx[start],
			                           bytes_to_read, file_offset);

			if (bytes_read != (ssize_t)bytes_to_read) {
				#pragma omp atomic write
				read_error = 1;
			}
		}
	}

	close(fd);

	if (read_error) {
		print_error(__func__, "parallel read failed", errno);
		csc_free_matrix(m);
		return NULL;
	}

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
 * @brief Load a CSC matrix from file (auto-detects format and uses best method).
 *
 * For .bin files, automatically uses parallel loading:
 *   - mmap with parallel copy (usually fastest)
 *   - Falls back to pread if mmap fails
 *
 * For .mtx files, uses parallel MTX loading if OpenMP is enabled.
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

	if (strcmp(ext, "bin") == 0) {
		return csc_load_bin(filename);
	} else if (strcmp(ext, "mtx") == 0) {
		return csc_load_mtx(filename);
	} else {
		print_error(__func__, "unsupported file extension", 0);
		return NULL;
	}
}

/**
 * @brief Load a distributed partition of a sparse binary matrix (.bin format only).
 *
 * Each MPI rank loads only its assigned partition of columns to minimize
 * memory usage and I/O bottlenecks. Uses block distribution strategy.
 *
 * Algorithm:
 * 1. All ranks read header independently (small, fast)
 * 2. Each rank calculates its column partition [start_col, end_col)
 * 3. Each rank reads full col_ptr to determine edge ranges
 * 4. Each rank reads only its portion of row_idx
 *
 * @param filename Path to the .bin file.
 * @param mpi_rank Current MPI rank.
 * @param mpi_size Total number of MPI ranks.
 * @return Newly allocated CSCBinaryMatrix (local partition), or NULL on error.
 */
CSCBinaryMatrix*
csc_load_matrix_distributed(const char *filename, int mpi_rank, int mpi_size)
{
	const char *ext = get_file_extension(filename);
	if (!ext) {
		if (mpi_rank == 0)
			print_error(__func__, "no file extension found", 0);
		return NULL;
	}

	/* For .bin format: efficient distributed loading */
	if (strcmp(ext, "bin") == 0) {
		/* Open file on all ranks */
		int fd = open(filename, O_RDONLY);
		if (fd < 0) {
			if (mpi_rank == 0)
				print_error(__func__, "failed to open .bin file", errno);
			return NULL;
		}

		/* Read header - all ranks read independently (small, fast) */
		uint32_t nrows_u32, ncols_u32;
		size_t nnz_total;
		
		if (read(fd, &nrows_u32, sizeof(uint32_t)) != sizeof(uint32_t) ||
		    read(fd, &ncols_u32, sizeof(uint32_t)) != sizeof(uint32_t) ||
		    read(fd, &nnz_total, sizeof(size_t)) != sizeof(size_t))
		{
			if (mpi_rank == 0)
				print_error(__func__, "failed to read header", errno);
			close(fd);
			return NULL;
		}

		size_t nrows = nrows_u32;
		size_t ncols = ncols_u32;

		/* Calculate this rank's column partition */
		size_t cols_per_rank = ncols / mpi_size;
		size_t remainder = ncols % mpi_size;
		
		size_t start_col = mpi_rank * cols_per_rank + (mpi_rank < (int)remainder ? mpi_rank : remainder);
		size_t end_col = start_col + cols_per_rank + (mpi_rank < (int)remainder ? 1 : 0);
		size_t local_ncols = end_col - start_col;

		/* Read full col_ptr array (needed to determine edge ranges) */
		uint32_t *full_col_ptr = malloc((ncols + 1) * sizeof(uint32_t));
		if (!full_col_ptr) {
			if (mpi_rank == 0)
				print_error(__func__, "malloc failed for col_ptr", errno);
			close(fd);
			return NULL;
		}

		off_t col_ptr_offset = sizeof(uint32_t) * 2 + sizeof(size_t);
		ssize_t col_ptr_bytes = (ncols + 1) * sizeof(uint32_t);
		
		if (pread(fd, full_col_ptr, col_ptr_bytes, col_ptr_offset) != col_ptr_bytes) {
			if (mpi_rank == 0)
				print_error(__func__, "failed to read col_ptr", errno);
			free(full_col_ptr);
			close(fd);
			return NULL;
		}

		/* Determine edge range for this rank's columns */
		uint32_t edge_start = full_col_ptr[start_col];
		uint32_t edge_end = full_col_ptr[end_col];
		size_t local_nnz = edge_end - edge_start;

		/* Allocate local matrix */
		CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
		if (!m) {
			if (mpi_rank == 0)
				print_error(__func__, "malloc failed", errno);
			free(full_col_ptr);
			close(fd);
			return NULL;
		}

		m->nrows = nrows;        /* Keep global row count */
		m->ncols = local_ncols;  /* Local column count */
		m->nnz = local_nnz;      /* Local edge count */

		/* Allocate local col_ptr (adjusted to local indexing) */
		m->col_ptr = malloc((local_ncols + 1) * sizeof(uint32_t));
		if (!m->col_ptr) {
			if (mpi_rank == 0)
				print_error(__func__, "malloc failed for local col_ptr", errno);
			free(full_col_ptr);
			free(m);
			close(fd);
			return NULL;
		}

		/* Copy and adjust col_ptr to local indexing */
		for (size_t i = 0; i <= local_ncols; i++) {
			m->col_ptr[i] = full_col_ptr[start_col + i] - edge_start;
		}
		free(full_col_ptr);

		/* Allocate local row_idx */
		m->row_idx = malloc(local_nnz * sizeof(uint32_t));
		if (!m->row_idx) {
			if (mpi_rank == 0)
				print_error(__func__, "malloc failed for row_idx", errno);
			csc_free_matrix(m);
			close(fd);
			return NULL;
		}

		/* Read this rank's portion of row_idx */
		off_t row_idx_offset = col_ptr_offset + col_ptr_bytes + edge_start * sizeof(uint32_t);
		ssize_t row_idx_bytes = local_nnz * sizeof(uint32_t);

		if (local_nnz > 0) {
			if (pread(fd, m->row_idx, row_idx_bytes, row_idx_offset) != row_idx_bytes) {
				if (mpi_rank == 0)
					print_error(__func__, "failed to read row_idx partition", errno);
				csc_free_matrix(m);
				close(fd);
				return NULL;
			}
		}

		close(fd);
		
		/* Debug output from all ranks */
		fprintf(stderr, "[Rank %d] Loaded partition: columns [%zu, %zu) = %zu cols, nnz=%zu\n",
		        mpi_rank, start_col, end_col, local_ncols, local_nnz);
		
		return m;
	}
	/* For .mtx format: rank 0 loads, others receive via MPI */
	else if (strcmp(ext, "mtx") == 0) {
		/* This is a simplified fallback - ideally should distribute after loading */
		if (mpi_rank == 0) {
			fprintf(stderr, "Warning: .mtx distributed loading not optimized\n");
			fprintf(stderr, "For best performance, convert to .bin format first\n");
		}
		
		/* For now, all ranks load full matrix (inefficient but functional) */
		/* TODO: Implement proper distribution */
		return csc_load_matrix(filename);
	}
	else {
		if (mpi_rank == 0)
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
	if (m->row_idx)
		free(m->row_idx);
	if (m->col_ptr)
		free(m->col_ptr);

	free(m);
	m = NULL;
}
