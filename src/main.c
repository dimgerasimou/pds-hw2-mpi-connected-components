/**
 * @file main.c
 * @brief Entry point for the connected components benchmark program.
 *
 * This implementation was developed for the purposes of the class:
 * Parallel and Distributed Systems,
 * Department of Electrical and Computer Engineering,
 * Aristotle University of Thessaloniki.
 *
 * Loads a sparse binary matrix in CSC format, runs a benchmark,
 * and prints the statistics.
 *
 * Usage: ./connected_components_mpi [-n n_trials] ./data_filepath
 */

#include "connected_components.h"
#include "matrix.h"
#include "error.h"
#include "args.h"
#include "benchmark.h"

const char *program_name = "connected_components_mpi";

int
main(int argc, char *argv[])
{
	CSCBinaryMatrix *matrix;
	Benchmark *benchmark = NULL;
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
	double t_load_start = now_sec();
	matrix = csc_load_matrix(filepath);
	double t_load_end = now_sec();
	if (!matrix)
		return 1;

	/* Initialize benchmarking structure */
	benchmark = benchmark_init(filepath, n_trials, matrix);
	if (!benchmark) {
		csc_free_matrix(matrix);
		return 1;
	}
	benchmark->matrix_info.load_time_s = t_load_end - t_load_start;
	/* Actually run the benchmark */
	ret = benchmark_cc(matrix, benchmark);

	benchmark_print(benchmark);

	/* Cleanup */
	csc_free_matrix(matrix);
	
	return ret;
}
