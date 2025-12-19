/**
 * @file matrix.c
 * @brief CSC (Compressed Sparse Column) binary matrix utilities.
 *
 * This module implements loading and freeing of sparse binary matrices in
 * CSC format.
 *
 * Only binary matrices are represented. Any non-zero numeric values in the
 * input are treated as 1.
 *
 * MPI extension: Supports distributed loading where each rank loads only
 * its partition of columns to minimize memory usage.
 */

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "error.h"
#include "matrix.h"

/* ------------------------------------------------------------------------- */
/*                            Static Helper Functions                        */
/* ------------------------------------------------------------------------- */

/**
 * @brief Robust pread() wrapper that reads exactly @p nbytes bytes.
 *
 * Some parallel/distributed filesystems can transiently return EINTR, EAGAIN,
 * or short reads under load. This helper retries with exponential backoff
 * (bounded) and continues until all bytes are read or a hard error occurs.
 *
 * @param fd File descriptor.
 * @param buf Destination buffer.
 * @param nbytes Number of bytes to read.
 * @param offset File offset.
 * @return 0 on success, -1 on failure with errno set.
 */
static int
pread_full(int fd, void *buf, size_t nbytes, off_t offset)
{
	unsigned char *p = (unsigned char *)buf;
	size_t done = 0;
	int retries = 0;

	/* Tunables: keep small but sufficient for shared FS hiccups */
	const int MAX_RETRIES = 20;
	const long BASE_BACKOFF_MS = 5;   /* 5ms */
	const long MAX_BACKOFF_MS  = 500; /* 0.5s per retry cap */

	while (done < nbytes) {
		ssize_t r = pread(fd, p + done, nbytes - done, offset + (off_t)done);

		if (r > 0) {
			done += (size_t)r;
			continue;
		}

		if (r == 0) {
			/* Unexpected EOF: file truncated or wrong offsets */
			errno = EIO;
			return -1;
		}

		/* r < 0 */
		if (errno == EINTR) {
			continue;
		}

		if (errno == EAGAIN) {
			if (++retries > MAX_RETRIES)
				return -1;

			long backoff = BASE_BACKOFF_MS << (retries < 10 ? retries : 10);
			if (backoff > MAX_BACKOFF_MS)
				backoff = MAX_BACKOFF_MS;

			struct timespec ts;
			ts.tv_sec = backoff / 1000;
			ts.tv_nsec = (backoff % 1000) * 1000000L;
			nanosleep(&ts, NULL);
			continue;
		}

		return -1;
	}

	return 0;
}

static int
has_bin_extension(const char *filename)
{
	const char *dot = strrchr(filename, '.');
	if (!dot)
		return 0;
	return (strcmp(dot + 1, "bin") == 0);
}

/* ------------------------------------------------------------------------- */
/*                            Public API Functions                           */
/* ------------------------------------------------------------------------- */

/**
 * @brief Load a distributed partition of a sparse binary matrix (.bin only).
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
csc_load_matrix(const char *filename, int mpi_rank, int mpi_size)
{
	/* Only .bin supported */
	if (!has_bin_extension(filename)) {
		if (mpi_rank == 0)
			print_error(__func__, "only .bin format is supported", 0);
		return NULL;
	}

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

	size_t nrows = (size_t)nrows_u32;
	size_t ncols = (size_t)ncols_u32;

	/* Calculate this rank's column partition */
	size_t cols_per_rank = ncols / (size_t)mpi_size;
	size_t remainder = ncols % (size_t)mpi_size;

	size_t start_col = (size_t)mpi_rank * cols_per_rank +
	    ((size_t)mpi_rank < remainder ? (size_t)mpi_rank : remainder);
	size_t end_col = start_col + cols_per_rank +
	    ((size_t)mpi_rank < remainder ? 1 : 0);
	size_t local_ncols = end_col - start_col;

	/* Read full col_ptr array (needed to determine edge ranges) */
	uint32_t *full_col_ptr = malloc((ncols + 1) * sizeof(uint32_t));
	if (!full_col_ptr) {
		if (mpi_rank == 0)
			print_error(__func__, "malloc failed for col_ptr", errno);
		close(fd);
		return NULL;
	}

	off_t col_ptr_offset = (off_t)(sizeof(uint32_t) * 2 + sizeof(size_t));
	size_t col_ptr_bytes = (ncols + 1) * sizeof(uint32_t);

	if (pread_full(fd, full_col_ptr, col_ptr_bytes, col_ptr_offset) != 0) {
		if (mpi_rank == 0)
			print_error(__func__, "failed to read col_ptr", errno);
		free(full_col_ptr);
		close(fd);
		return NULL;
	}

	/* Basic sanity checks */
	if (full_col_ptr[0] != 0 || full_col_ptr[ncols] != (uint32_t)nnz_total) {
		if (mpi_rank == 0)
			print_error(__func__, "corrupt col_ptr header", 0);
		free(full_col_ptr);
		close(fd);
		return NULL;
	}

	for (size_t j = 0; j < ncols; j++) {
		if (full_col_ptr[j] > full_col_ptr[j + 1]) {
			if (mpi_rank == 0)
				print_error(__func__, "col_ptr not monotone", 0);
			free(full_col_ptr);
			close(fd);
			return NULL;
		}
	}

	/* Determine edge range for this rank's columns */
	uint32_t edge_start = full_col_ptr[start_col];
	uint32_t edge_end = full_col_ptr[end_col];
	if (edge_end < edge_start || edge_end > (uint32_t)nnz_total) {
		if (mpi_rank == 0)
			print_error(__func__, "bad edge range", 0);
		free(full_col_ptr);
		close(fd);
		return NULL;
	}

	size_t local_nnz = (size_t)(edge_end - edge_start);

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
	m->ncols_global = ncols;

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
	m->row_idx = NULL;
	if (local_nnz > 0) {
		m->row_idx = malloc(local_nnz * sizeof(uint32_t));
		if (!m->row_idx) {
			if (mpi_rank == 0)
				print_error(__func__, "malloc failed for row_idx", errno);
			csc_free_matrix(m);
			close(fd);
			return NULL;
		}
	}

	/* Read this rank's portion of row_idx */
	off_t row_idx_offset = col_ptr_offset + (off_t)col_ptr_bytes +
	    (off_t)edge_start * (off_t)sizeof(uint32_t);
	size_t row_idx_bytes = local_nnz * sizeof(uint32_t);

	if (local_nnz > 0) {
		if (pread_full(fd, m->row_idx, row_idx_bytes, row_idx_offset) != 0) {
			if (mpi_rank == 0)
				print_error(__func__, "failed to read row_idx partition", errno);
			csc_free_matrix(m);
			close(fd);
			return NULL;
		}
	}

	close(fd);
	return m;
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

