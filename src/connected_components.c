/**
 * @file connected_components.c
 * @brief Hybrid MPI+OpenMP implementation for computing connected components.
 *
 * This module implements a distributed parallel algorithm for finding connected
 * components in an undirected graph using MPI for inter-node communication
 * and OpenMP for intra-node parallelization.
 *
 * Key optimizations:
 * - Each rank loads only its partition of the graph
 * - Only changed labels are communicated between ranks
 * - Uses existing Union-Find algorithm with minimal modifications
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "connected_components.h"

/* ========================================================================== */
/*                           UNION-FIND UTILITIES                             */
/* ========================================================================== */

/**
 * @brief Finds the root of a node with path compression.
 *
 * Traverses parent pointers until reaching the root, and compresses the
 * path by making all nodes along the path point directly to the root.
 * The early-exit optimization avoids redundant writes if the path is
 * already compressed.
 *
 * @param label Array of parent pointers representing disjoint sets
 * @param x Node index to find the root for
 * @return Root of the set containing x
 */
static inline uint32_t
find_compress(uint32_t *label, uint32_t x)
{
	uint32_t root = x;
	
	/* Find the root */
	while (label[root] != root)
		root = label[root];
	
	/* Compress the path */
	while (x != root) {
		uint32_t next = label[x];
		if (label[x] == next)
			break;  /* Already compressed */
		label[x] = root;
		x = next;
	}
	
	return root;
}

/**
 * @brief Unites two disjoint sets using Rem's algorithm.
 *
 * Implements lock-free parallel union-find using compare-and-swap (CAS)
 * operations. Retries up to MAX_RETRIES times before falling back to a
 * simpler atomic store. Canonical ordering (smaller index as root) ensures
 * deterministic results in parallel execution.
 *
 * @param label Array of parent pointers representing disjoint sets
 * @param a First node
 * @param b Second node
 */
static inline void
union_rem(uint32_t *label, uint32_t a, uint32_t b)
{
	const int MAX_RETRIES = 10;
	
	/* Retry loop with CAS operations */
	for (int retry = 0; retry < MAX_RETRIES; retry++) {
		a = find_compress(label, a);
		b = find_compress(label, b);
		
		if (a == b)
			return;
		
		/* Canonical ordering: smaller index as root */
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
	
	/* Fallback after maximum retries */
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

/* ========================================================================== */
/*                    DISTRIBUTED CONNECTED COMPONENTS                        */
/* ========================================================================== */

/**
 * @brief Pack changed labels for MPI communication.
 *
 * Identifies labels that have changed and packs them into send buffer.
 * This optimization reduces communication volume significantly.
 *
 * @param label Current label array
 * @param prev_label Previous label array (from last iteration)
 * @param n Number of local vertices
 * @param global_offset Starting vertex ID for this rank
 * @param send_buf Output buffer for (vertex_id, label) pairs
 * @return Number of changed labels
 */
static size_t
pack_changed_labels(const uint32_t *label, const uint32_t *prev_label,
                    uint32_t n, uint32_t global_offset, uint32_t *send_buf)
{
	size_t count = 0;
	
	/* Sequential packing - check ONLY this rank's partition */
	for (uint32_t i = 0; i < n; i++) {
		uint32_t global_id = global_offset + i;  /* Convert to global vertex ID */
		if (label[global_id] != prev_label[global_id]) {
			send_buf[count * 2] = global_id;         /* vertex ID (global) */
			send_buf[count * 2 + 1] = label[global_id];  /* new label */
			count++;
		}
	}
	
	return count;
}

/**
 * @brief Computes connected components using distributed parallel union-find.
 *
 * Algorithm phases (per iteration):
 * 1. Local union-find on partition edges (OpenMP parallel)
 * 2. Pack only changed labels for communication
 * 3. Exchange changed labels via MPI_Allgatherv
 * 4. Apply received label updates
 * 5. Check global convergence
 *
 * @param matrix Local CSC partition (only this rank's columns)
 * @param mpi_rank Current MPI rank
 * @param mpi_size Total number of MPI ranks
 * @return Number of connected components, or -1 on error
 */
int
connected_components_mpi(const CSCBinaryMatrix *matrix, int mpi_rank, int mpi_size)
{
	if (!matrix || matrix->nrows == 0)
		return 0;
	
	/* Fast path: single rank = use original OpenMP-only algorithm */
	if (mpi_size == 1) {
		return connected_components(matrix);
	}
	
	/* For distributed case:
	 * - Each rank owns a partition of COLUMNS (vertices)
	 * - matrix->ncols = number of local columns (vertices) this rank owns
	 * - matrix->nrows = total number of rows (may equal total vertices for square matrices)
	 * 
	 * CRITICAL: We must calculate global_offset the SAME WAY as the loader!
	 * The loader partitions based on the ORIGINAL ncols from the file.
	 * After loading, matrix->nrows tells us the original total.
	 */
	const uint32_t n = (uint32_t)matrix->nrows;  /* total vertices (global) */
	const uint32_t local_n = (uint32_t)matrix->ncols;  /* local vertices owned by this rank */
	
	/* Calculate global_offset using SAME logic as csc_load_matrix_distributed */
	/* It partitions based on total columns (which = nrows for square adjacency matrix) */
	uint32_t cols_per_rank = n / mpi_size;
	uint32_t remainder = n % mpi_size;
	uint32_t global_offset = mpi_rank * cols_per_rank + 
	                         (mpi_rank < (int)remainder ? mpi_rank : (int)remainder);
	
	/* Allocate label arrays */
	uint32_t *label = malloc(n * sizeof(uint32_t));
	uint32_t *prev_label = malloc(n * sizeof(uint32_t));
	
	if (!label || !prev_label) {
		free(label);
		free(prev_label);
		return -1;
	}
	
	/* Initialize: each node as its own parent */
	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n; i++) {
		label[i] = i;
	}
	
	/* Initialize prev_label to match initial state */
	memcpy(prev_label, label, n * sizeof(uint32_t));
	
	/* Allocate buffers for changed labels communication */
	uint32_t max_changes = local_n * 2;  /* (vertex_id, label) pairs */
	uint32_t *send_buf = malloc(max_changes * sizeof(uint32_t));
	
	if (!send_buf) {
		free(label);
		free(prev_label);
		free(send_buf);
		return -1;
	}
	
	/* Iterative convergence loop */
	const int MAX_ITERATIONS = 1000;
	int converged = 0;
	
	for (int iter = 0; iter < MAX_ITERATIONS && !converged; iter++) {
		/* Phase 1: Local union-find on partition edges */
		#pragma omp parallel
		{
			#pragma omp for schedule(dynamic, 128) nowait
			for (uint32_t col = 0; col < matrix->ncols; col++) {
				uint32_t start = matrix->col_ptr[col];
				uint32_t end = matrix->col_ptr[col + 1];
				uint32_t global_col = global_offset + col;
				
				for (uint32_t j = start; j < end; j++) {
					uint32_t row = matrix->row_idx[j];
					if (row < n)
						union_rem(label, row, global_col);
				}
			}
		}
		
		/* Phase 2: Compress paths locally */
		#pragma omp parallel for schedule(static, 2048)
		for (uint32_t i = 0; i < n; i++)
			find_compress(label, i);
		
		/* Phase 3: Pack changed labels only */
		size_t num_changed = pack_changed_labels(label, prev_label, local_n,
		                                          global_offset, send_buf);
		
		/* Phase 4: Exchange changed labels via MPI_Allgatherv */
		int *recvcounts = malloc(mpi_size * sizeof(int));
		int *displs = malloc(mpi_size * sizeof(int));
		
		if (!recvcounts || !displs) {
			free(label);
			free(prev_label);
			free(send_buf);
			free(recvcounts);
			free(displs);
			return -1;
		}
		
		int send_count = (int)num_changed * 2;  /* 2 ints per change */
		MPI_Allgather(&send_count, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
		
		/* Calculate displacements and total receive size */
		int total_recv = 0;
		displs[0] = 0;
		for (int i = 0; i < mpi_size; i++) {
			if (i > 0)
				displs[i] = displs[i-1] + recvcounts[i-1];
			total_recv += recvcounts[i];
		}
		
		uint32_t *recv_buf = malloc(total_recv * sizeof(uint32_t));
		if (!recv_buf) {
			free(label);
			free(prev_label);
			free(send_buf);
			free(recvcounts);
			free(displs);
			return -1;
		}
		
		MPI_Allgatherv(send_buf, send_count, MPI_UINT32_T,
		               recv_buf, recvcounts, displs, MPI_UINT32_T,
		               MPI_COMM_WORLD);
		
		/* Phase 5: Apply received label updates */
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < total_recv / 2; i++) {
			uint32_t vertex = recv_buf[i * 2];
			uint32_t new_label = recv_buf[i * 2 + 1];
			
			if (vertex < n) {
				/* Apply update using union operation */
				union_rem(label, vertex, new_label);
			}
		}
		
		free(recv_buf);
		free(recvcounts);
		free(displs);
		
		/* Phase 6: Check convergence (any changes globally?) */
		int local_changed = (num_changed > 0) ? 1 : 0;
		int global_changed = 0;
		MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
		
		if (!global_changed) {
			converged = 1;
		}
		
		/* Update previous labels for next iteration */
		memcpy(prev_label, label, n * sizeof(uint32_t));
		
	}
	
	/* Final compression pass on ALL vertices (needed for correctness) */
	#pragma omp parallel for schedule(static, 2048)
	for (uint32_t i = 0; i < n; i++)
		find_compress(label, i);
	
	/* Count roots ONLY in this rank's partition */
	uint32_t local_count = 0;
	#pragma omp parallel for reduction(+:local_count) schedule(static, 2048)
	for (uint32_t i = 0; i < local_n; i++) {
		uint32_t global_id = global_offset + i;
		if (label[global_id] == global_id)
			local_count++;
	}
	
	/* Sum component counts across all ranks */
	uint32_t global_count = 0;
	MPI_Reduce(&local_count, &global_count, 1, MPI_UINT32_T, MPI_SUM,
	           0, MPI_COMM_WORLD);
	
	/* Broadcast result to all ranks */
	MPI_Bcast(&global_count, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
	
	free(label);
	free(prev_label);
	free(send_buf);
	
	return (int)global_count;
}

/**
 * @brief Wrapper function for backwards compatibility.
 *
 * For non-MPI builds or single-process execution.
 */
int
connected_components(const CSCBinaryMatrix *matrix)
{
	/* Single-process version - original algorithm */
	if (!matrix || matrix->nrows == 0)
		return 0;
	
	const uint32_t n = (uint32_t)matrix->nrows;
	uint32_t *label = malloc(n * sizeof(uint32_t));
	if (!label)
		return -1;
	
	/* Initialize: each node as its own parent */
	#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n; i++)
		label[i] = i;
	
	/* Process all edges: union connected nodes */
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
	
	/* Final compression pass: flatten all paths */
	#pragma omp parallel for schedule(static, 2048)
	for (uint32_t i = 0; i < n; i++)
		find_compress(label, i);
	
	/* Count roots (each root represents one component) */
	uint32_t count = 0;
	#pragma omp parallel for reduction(+:count) schedule(static, 2048)
	for (uint32_t i = 0; i < n; i++)
		if (label[i] == i)
			count++;
	
	free(label);
	return (int)count;
}
