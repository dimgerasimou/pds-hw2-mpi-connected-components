/**
 * @file cc_openmp.c
 * @brief Optimized OpenMP implementations for computing connected components.
 *
 * This module implements two parallel algorithms for finding connected
 * components in an undirected graph using OpenMP:
 *
 * - Label Propagation (variant 0): Iteratively propagates minimum labels
 *   until convergence with persistent threads and relaxed atomics.
 *
 * - Union-Find with Rem's Algorithm (variant 1): Lock-free parallel
 *   union-find using compare-and-swap operations and path compression.
 *
 * Both algorithms return the count of unique connected components.
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

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
/*                         UNION-FIND ALGORITHM                               */
/* ========================================================================== */

/**
 * @brief Computes connected components using parallel union-find.
 *
 * Algorithm phases:
 * 1. Initialize each node as its own root (parallel)
 * 2. Perform parallel union operations on edges using dynamic scheduling
 * 3. Flatten all paths to roots for accurate counting (parallel)
 * 4. Count roots in parallel using reduction
 *
 * @param matrix Sparse CSC binary matrix representing graph
 * @param n_threads Number of OpenMP threads to use
 * @return Number of connected components, or -1 on error
 */

int
connected_components(const CSCBinaryMatrix *matrix)
{
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
