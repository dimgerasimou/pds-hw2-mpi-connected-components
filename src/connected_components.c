/**
 * @file connected_components.c
 * @brief Hybrid MPI+OpenMP connected components.
 *
 * OpenMP:
 *   - Union-Find with lock-free unions (unchanged).
 *
 * MPI+OpenMP :
 *   - Deterministic synchronous min-label propagation (LP)
 *   - Uses MPI_Allgatherv to get a full label snapshot each round
 *   - Light pointer-jumping to reduce iterations
 *
 * REQUIREMENT:
 *   Input .bin should store explicit undirected adjacency (both directions),
 *   so that vertex u's neighbors appear in column u.
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <omp.h>
#include <mpi.h>

#include "connected_components.h"

static inline uint32_t
min_u32(uint32_t a, uint32_t b)
{
	return (a < b) ? a : b;
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

int
connected_components(const CSCBinaryMatrix *matrix, int mpi_rank, int mpi_size)
{
	if (!matrix || matrix->nrows == 0)
		return 0;
	
	const uint32_t n = (uint32_t)matrix->ncols_global;

	uint32_t global_offset, local_n;
	compute_partition_u32(n, mpi_size, mpi_rank, &global_offset, &local_n);

	if ((uint32_t)matrix->ncols != local_n)
		return -1;

	int *recvcounts = malloc((size_t)mpi_size * sizeof(int));
	int *displs     = malloc((size_t)mpi_size * sizeof(int));
	if (!recvcounts || !displs) {
		free(recvcounts);
		free(displs);
		return -1;
	}

	for (int r = 0; r < mpi_size; r++) {
		uint32_t off_r, n_r;
		compute_partition_u32(n, mpi_size, r, &off_r, &n_r);
		recvcounts[r] = (int)n_r;
		displs[r] = (int)off_r;
	}

	uint32_t *label_local = malloc((size_t)local_n * sizeof(uint32_t));
	uint32_t *next_local  = malloc((size_t)local_n * sizeof(uint32_t));
	uint32_t *label_global = malloc((size_t)n * sizeof(uint32_t));
	if (!label_local || !next_local || !label_global) {
		free(recvcounts);
		free(displs);
		free(label_local);
		free(next_local);
		free(label_global);
		return -1;
	}

	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < local_n; i++)
		label_local[i] = global_offset + i;

	/* Initial global snapshot */
	MPI_Allgatherv(label_local, (int)local_n, MPI_UINT32_T,
	               label_global, recvcounts, displs, MPI_UINT32_T,
	               MPI_COMM_WORLD);

	const int MAX_ITER = 512;
	int converged = 0;

	for (int iter = 0; iter < MAX_ITER && !converged; iter++) {
		int local_changed = 0;

		#pragma omp parallel for schedule(guided) reduction(|:local_changed)
		for (uint32_t li = 0; li < local_n; li++) {
			uint32_t u = global_offset + li;
			uint32_t best = label_global[u];

			uint32_t start = matrix->col_ptr[li];
			uint32_t end   = matrix->col_ptr[li + 1];

			for (uint32_t j = start; j < end; j++) {
				uint32_t v = matrix->row_idx[j];
				if (v < n)
					best = min_u32(best, label_global[v]);
			}

			best = min_u32(best, label_global[best]);

			next_local[li] = best;
			if (best != label_local[li])
				local_changed = 1;
		}

		uint32_t *tmp = label_local;
		label_local = next_local;
		next_local = tmp;

		/* publish updated local labels to global snapshot */
		MPI_Allgatherv(label_local, (int)local_n, MPI_UINT32_T,
		               label_global, recvcounts, displs, MPI_UINT32_T,
		               MPI_COMM_WORLD);

		int global_changed = 0;
		MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
		if (!global_changed)
			converged = 1;
	}

	/* Final compression so counting roots is accurate */
	#pragma omp parallel for schedule(static)
	for (uint32_t li = 0; li < local_n; li++) {
		uint32_t x = label_local[li];
		for (int k = 0; k < 8; k++)
			x = label_global[x];
		label_local[li] = x;
	}

	MPI_Allgatherv(label_local, (int)local_n, MPI_UINT32_T,
	               label_global, recvcounts, displs, MPI_UINT32_T,
	               MPI_COMM_WORLD);

	uint32_t local_count = 0;
	#pragma omp parallel for schedule(static) reduction(+:local_count)
	for (uint32_t li = 0; li < local_n; li++) {
		uint32_t u = global_offset + li;
		if (label_global[u] == u)
			local_count++;
	}

	uint32_t global_count = 0;
	MPI_Reduce(&local_count, &global_count, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global_count, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	free(recvcounts);
	free(displs);
	free(label_local);
	free(next_local);
	free(label_global);

	return (int)global_count;
}

