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

static int
pread_full(int fd, void *buf, size_t nbytes, off_t offset)
{
	unsigned char *p = (unsigned char *)buf;
	size_t done = 0;
	int retries = 0;

	const int MAX_RETRIES = 20;
	const long BASE_BACKOFF_MS = 5;
	const long MAX_BACKOFF_MS  = 500;

	while (done < nbytes) {
		ssize_t r = pread(fd, p + done, nbytes - done, offset + (off_t)done);

		if (r > 0) {
			done += (size_t)r;
			continue;
		}

		if (r == 0) {
			errno = EIO;
			return -1;
		}

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

CSCBinaryMatrix*
csc_load_matrix(const char *filename, int mpi_rank, int mpi_size)
{
	if (!has_bin_extension(filename)) {
		if (mpi_rank == 0)
			print_error(__func__, "only .bin format is supported", 0);
		return NULL;
	}

	int fd = open(filename, O_RDONLY);
	if (fd < 0) {
		if (mpi_rank == 0)
			print_error(__func__, "failed to open .bin file", errno);
		return NULL;
	}

	/* Header */
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

	/* Block partition of columns */
	size_t cols_per_rank = ncols / (size_t)mpi_size;
	size_t remainder = ncols % (size_t)mpi_size;

	size_t start_col = (size_t)mpi_rank * cols_per_rank +
	    ((size_t)mpi_rank < remainder ? (size_t)mpi_rank : remainder);
	size_t end_col = start_col + cols_per_rank +
	    ((size_t)mpi_rank < remainder ? 1 : 0);
	size_t local_ncols = end_col - start_col;

	off_t col_ptr_offset = (off_t)(sizeof(uint32_t) * 2 + sizeof(size_t));

	/* ------------------------------------------------------------------ */
	/* Read ONLY needed col_ptr window: col_ptr[start_col .. end_col]     */
	/* ------------------------------------------------------------------ */
	uint32_t *win_col_ptr = malloc((local_ncols + 1) * sizeof(uint32_t));
	if (!win_col_ptr) {
		if (mpi_rank == 0)
			print_error(__func__, "malloc failed for col_ptr window", errno);
		close(fd);
		return NULL;
	}

	off_t win_off = col_ptr_offset + (off_t)start_col * (off_t)sizeof(uint32_t);
	size_t win_bytes = (local_ncols + 1) * sizeof(uint32_t);

	if (pread_full(fd, win_col_ptr, win_bytes, win_off) != 0) {
		if (mpi_rank == 0)
			print_error(__func__, "failed to read col_ptr window", errno);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	/* Minimal sanity: read col_ptr[0] and col_ptr[ncols] only */
	uint32_t first_ptr, last_ptr;
	if (pread_full(fd, &first_ptr, sizeof(uint32_t), col_ptr_offset) != 0) {
		if (mpi_rank == 0)
			print_error(__func__, "failed to read col_ptr[0]", errno);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	off_t last_off = col_ptr_offset + (off_t)ncols * (off_t)sizeof(uint32_t);
	if (pread_full(fd, &last_ptr, sizeof(uint32_t), last_off) != 0) {
		if (mpi_rank == 0)
			print_error(__func__, "failed to read col_ptr[ncols]", errno);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	if (first_ptr != 0 || last_ptr != (uint32_t)nnz_total) {
		if (mpi_rank == 0)
			print_error(__func__, "corrupt col_ptr header", 0);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	/* Window sanity: monotone in our local window */
	for (size_t j = 0; j < local_ncols; j++) {
		if (win_col_ptr[j] > win_col_ptr[j + 1]) {
			if (mpi_rank == 0)
				print_error(__func__, "col_ptr not monotone", 0);
			free(win_col_ptr);
			close(fd);
			return NULL;
		}
	}

	uint32_t edge_start = win_col_ptr[0];
	uint32_t edge_end   = win_col_ptr[local_ncols];

	if (edge_end < edge_start || edge_end > (uint32_t)nnz_total) {
		if (mpi_rank == 0)
			print_error(__func__, "bad edge range", 0);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	size_t local_nnz = (size_t)(edge_end - edge_start);

	CSCBinaryMatrix *m = malloc(sizeof(CSCBinaryMatrix));
	if (!m) {
		if (mpi_rank == 0)
			print_error(__func__, "malloc failed", errno);
		free(win_col_ptr);
		close(fd);
		return NULL;
	}

	m->nrows = nrows;
	m->ncols = local_ncols;
	m->ncols_global = ncols;
	m->nnz = local_nnz;

	m->col_ptr = malloc((local_ncols + 1) * sizeof(uint32_t));
	if (!m->col_ptr) {
		if (mpi_rank == 0)
			print_error(__func__, "malloc failed for local col_ptr", errno);
		free(win_col_ptr);
		free(m);
		close(fd);
		return NULL;
	}

	for (size_t i = 0; i <= local_ncols; i++)
		m->col_ptr[i] = win_col_ptr[i] - edge_start;

	free(win_col_ptr);

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
	size_t col_ptr_bytes_total = (ncols + 1) * sizeof(uint32_t);
	off_t row_idx_offset = col_ptr_offset + (off_t)col_ptr_bytes_total +
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

