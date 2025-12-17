/**
 * @file main.c
 * @brief Entry point for the connected components benchmark program.
 *
 * This implementation was developed for the purposes of the class:
 * Parallel and Distributed Systems,
 * Department of Electrical and Computer Engineering,
 * Aristotle University of Thessaloniki.
 *
 * Loads a sparse binary matrix in CSC format, based on the selected
 * connected components implementation (sequential or parallel),
 * runs a benchmark, and prints the statistics. The implementations
 * are selected through a series of definitions (through compiler flags).
 *
 * Supported implementations:
 * - USE_SEQUENTIAL
 * - USE_OPENMP
 * - USE_PTHREADS
 * - USE_CILK
 *
 * Usage: ./connected_components [-t n_threads] [-n n_trials] ./data_filepath
 */

#include "connected_components.h"
#include "matrix.h"
#include "error.h"
#include "args.h"

#include <stdio.h>

const char *program_name = "connected_components_mpi";

int
main(int argc, char *argv[])
{
	CSCBinaryMatrix *matrix;
	char *filepath;
	unsigned int n_trials;
	int ret = 0;

	/* Initialize program name for error reporting */
	set_program_name(argv[0]);

	/* Parse command line arguments */
	if (parseargs(argc, argv, &n_trials, &filepath)) {
		return 1;
	}
	
	/* Load the sparse matrix */
	matrix = csc_load_matrix(filepath);
	if (!matrix)
		return 1;

	ret = connected_components(matrix);

	if (ret > 0) {
		printf("CC:%d\n", ret);
		ret = 0;
	} else {
		ret = 1;
	}

	/* Cleanup */
	csc_free_matrix(matrix);
	
	return ret;
}
