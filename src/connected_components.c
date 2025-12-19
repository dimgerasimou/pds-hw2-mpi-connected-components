/**
 * @file connected_components.c
 * @brief Hybrid MPI+OpenMP implementations for computing connected components.
 *
 * Single MPI algorithm (speed/complexity sweet spot):
 *   - Deterministic synchronous min-label propagation (LP)
 *   - Sparse ghost exchange via MPI_Alltoallv (no global MPI_Allgatherv(label))
 *
 * Single-process OpenMP implementation remains unchanged.
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <omp.h>
#include <mpi.h>

#include "connected_components.h"

/* ========================================================================== */
/*                    ORIGINAL OPENMP-ONLY VERSION (UNCHANGED)                */
/* ========================================================================== */

static inline uint32_t
find_compress(uint32_t *label, uint32_t x)
{
	uint32_t root = x;

	while (label[root] != root)
		root = label[root];

	while (x != root) {
		uint32_t next = label[x];
		label[x] = root;
		x = next;
	}

	return root;
}

static inline void
union_rem(uint32_t *label, uint32_t a, uint32_t b)
{
	const int MAX_RETRIES = 10;

	for (int retry = 0; retry < MAX_RETRIES; retry++) {
		a = find_compress(label, a);
		b = find_compress(label, b);

		if (a == b)
			return;

		if (a > b) {
			uint32_t temp = a;
			a = b;
			b = temp;
		}

		uint32_t expected = b;
		if (__atomic_compare_exchange_n(&label[b], &expected, a,
		                                0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
			return;
		}

		b = expected;
	}

	a = find_compress(label, a);
	b = find_compress(label, b);
	if (a != b) {
		if (a > b) {
			uint32_t temp = a;
			a = b;
			b = temp;
		}
		__atomic_store_n(&label[b], a, __ATOMIC_RELEASE);
	}
}

static int
connected_components_omp(const CSCBinaryMatrix *matrix)
{
	if (!matrix || matrix->nrows == 0)
		return 0;

	const uint32_t n = (uint32_t)matrix->nrows;
	uint32_t *label = malloc(n * sizeof(uint32_t));
	if (!label)
		return -1;

	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n; i++)
		label[i] = i;

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, 128) nowait
		for (uint32_t col = 0; col < matrix->ncols; col++) {
			uint32_t start = matrix->col_ptr[col];
			uint32_t end = matrix->col_ptr[col + 1];

			for (uint32_t j = start; j < end; j++) {
				uint32_t row = matrix->row_idx[j];
				if (row < n)
					union_rem(label, row, col);
			}
		}
	}

	#pragma omp parallel for schedule(static, 2048)
	for (uint32_t i = 0; i < n; i++)
		find_compress(label, i);

	uint32_t count = 0;
	#pragma omp parallel for reduction(+:count) schedule(static, 2048)
	for (uint32_t i = 0; i < n; i++)
		if (label[i] == i)
			count++;

	free(label);
	return (int)count;
}

/* ========================================================================== */
/*                                MPI VERSION                                 */
/* ========================================================================== */

static inline uint32_t
min_u32(uint32_t a, uint32_t b)
{
	return (a < b) ? a : b;
}

static inline uint32_t
max_u32(uint32_t a, uint32_t b)
{
	return (a > b) ? a : b;
}

static inline void
compute_partition_u32(uint32_t n, int mpi_size, int mpi_rank,
                      uint32_t *global_offset, uint32_t *local_n)
{
	uint32_t base = n / (uint32_t)mpi_size;
	uint32_t rem  = n % (uint32_t)mpi_size;

	*global_offset = (uint32_t)mpi_rank * base +
	                 (mpi_rank < (int)rem ? (uint32_t)mpi_rank : rem);
	*local_n = base + (mpi_rank < (int)rem ? 1u : 0u);
}

static inline int
is_local_u32(uint32_t v, uint32_t global_offset, uint32_t local_n)
{
	return (v >= global_offset) && (v < global_offset + local_n);
}

/* Same owner() you already effectively use */
static inline int
owner_rank_u32(uint32_t v, uint32_t n, int mpi_size)
{
	uint32_t base = n / (uint32_t)mpi_size;
	uint32_t rem  = n % (uint32_t)mpi_size;
	uint32_t big  = base + 1;

	if (v < rem * big)
		return (int)(v / big);

	return (int)(rem + (v - rem * big) / base);
}

/* A small bounded "find root" that is correct and fast once compressed */
static inline uint32_t
find_root_bounded(uint32_t *parent, uint32_t x)
{
	/* After a few pointer jumps per round, depth is small */
	for (int k = 0; k < 8; k++) {
		uint32_t px = parent[x];
		if (px == x)
			return x;
		uint32_t ppx = parent[px];
		if (ppx == px)
			return px;
		x = ppx;
	}
	/* Fallback: finish */
	while (parent[x] != x)
		x = parent[x];
	return x;
}

int
connected_components(const CSCBinaryMatrix *matrix, int mpi_rank, int mpi_size)
{
	if (!matrix || matrix->nrows == 0)
		return 0;

	if (mpi_size == 1)
		return connected_components_omp(matrix);

	const uint32_t n = (uint32_t)matrix->nrows;

	uint32_t global_offset, local_n;
	compute_partition_u32(n, mpi_size, mpi_rank, &global_offset, &local_n);

	/* Allgatherv layout */
	int *recvcounts = malloc((size_t)mpi_size * sizeof(int));
	int *displs     = malloc((size_t)mpi_size * sizeof(int));
	if (!recvcounts || !displs)
		return -1;

	for (int r = 0; r < mpi_size; r++) {
		uint32_t off_r, n_r;
		compute_partition_u32(n, mpi_size, r, &off_r, &n_r);
		recvcounts[r] = (int)n_r;
		displs[r] = (int)off_r;
	}

	uint32_t *parent = malloc((size_t)n * sizeof(uint32_t));
	uint32_t *next_local = malloc((size_t)local_n * sizeof(uint32_t));
	if (!parent || !next_local) {
		free(recvcounts);
		free(displs);
		free(parent);
		free(next_local);
		return -1;
	}

	/* init local */
	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < local_n; i++)
		parent[global_offset + i] = global_offset + i;

	MPI_Allgatherv(parent + global_offset, (int)local_n, MPI_UINT32_T,
	               parent, recvcounts, displs, MPI_UINT32_T, MPI_COMM_WORLD);

	const int MAX_ROUNDS = 80;
	int converged = 0;

	for (int round = 0; round < MAX_ROUNDS && !converged; round++) {
		/* ------------------------------------------------------------ */
		/* 1) Pointer jumping (a few passes)                             */
		/* ------------------------------------------------------------ */
		#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < local_n; i++) {
			uint32_t u = global_offset + i;
			/* 3 jumps typically enough */
			parent[u] = parent[parent[u]];
			parent[u] = parent[parent[u]];
			parent[u] = parent[parent[u]];
		}

		/* Sync after compression so "find_root" is consistent */
		MPI_Allgatherv(parent + global_offset, (int)local_n, MPI_UINT32_T,
		               parent, recvcounts, displs, MPI_UINT32_T, MPI_COMM_WORLD);

		/* ------------------------------------------------------------ */
		/* 2) Owner-hooking: only owner of the HIGHER root updates it     */
		/*    next_local stores proposed new parent for owned vertices    */
		/* ------------------------------------------------------------ */
		#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < local_n; i++) {
			uint32_t u = global_offset + i;
			next_local[i] = parent[u];
		}

		int local_changed = 0;

		#pragma omp parallel for schedule(guided) reduction(|:local_changed)
		for (uint32_t col = 0; col < local_n; col++) {
			uint32_t u = global_offset + col;
			uint32_t ru = find_root_bounded(parent, u);

			uint32_t start = matrix->col_ptr[col];
			uint32_t end   = matrix->col_ptr[col + 1];

			for (uint32_t j = start; j < end; j++) {
				uint32_t v = matrix->row_idx[j];
				if (v >= n)
					continue;

				uint32_t rv = find_root_bounded(parent, v);
				if (ru == rv)
					continue;

				uint32_t low  = min_u32(ru, rv);
				uint32_t high = max_u32(ru, rv);

				/* Only the owner of 'high' is allowed to change parent[high] */
				if (owner_rank_u32(high, n, mpi_size) != mpi_rank)
					continue;

				uint32_t idx = high - global_offset;
				/* high must be local if we are its owner */
				if (idx < local_n) {
					uint32_t cur = next_local[idx];
					uint32_t newp = min_u32(cur, low);
					if (newp != cur) {
						next_local[idx] = newp;
						local_changed = 1;
					}
				}
			}
		}

		/* Commit owned changes */
		#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < local_n; i++)
			parent[global_offset + i] = next_local[i];

		/* Sync parents */
		MPI_Allgatherv(parent + global_offset, (int)local_n, MPI_UINT32_T,
		               parent, recvcounts, displs, MPI_UINT32_T, MPI_COMM_WORLD);

		int global_changed = 0;
		MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
		if (!global_changed)
			converged = 1;
	}

	/* Final strong compression */
	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < local_n; i++) {
		uint32_t u = global_offset + i;
		for (int k = 0; k < 8; k++)
			parent[u] = parent[parent[u]];
	}

	MPI_Allgatherv(parent + global_offset, (int)local_n, MPI_UINT32_T,
	               parent, recvcounts, displs, MPI_UINT32_T, MPI_COMM_WORLD);

	/* Count roots */
	uint32_t local_count = 0;
	#pragma omp parallel for schedule(static) reduction(+:local_count)
	for (uint32_t i = 0; i < local_n; i++) {
		uint32_t u = global_offset + i;
		if (parent[u] == u)
			local_count++;
	}

	uint32_t global_count = 0;
	MPI_Reduce(&local_count, &global_count, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global_count, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	free(parent);
	free(next_local);
	free(recvcounts);
	free(displs);

	return (int)global_count;
}


