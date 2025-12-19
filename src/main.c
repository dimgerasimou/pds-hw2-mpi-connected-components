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
 * Usage: mpirun -np N ./connected_components_mpi [-n n_trials] ./data_filepath
 */

#include <mpi.h>
#include "connected_components.h"
#include "matrix.h"
#include "error.h"
#include "args.h"
#include "benchmark.h"

#include <stdio.h>

const char *program_name = "connected_components_mpi";

int
main(int argc, char *argv[])
{
	CSCBinaryMatrix *matrix;
	Benchmark *benchmark = NULL;
	char *filepath;
	unsigned int n_trials;
	int ret = 0;
	int mpi_rank, mpi_size;
	int provided;

	/* Initialize MPI with thread support */
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED) {
		print_error(__func__, "MPI does not support required threading level\n", 0);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	/* Initialize program name for error reporting */
	set_program_name(argv[0]);

	/* Parse command line arguments */
	if (parseargs(argc, argv, &n_trials, &filepath)) {
		MPI_Finalize();
		return 1;
	}
	
	/* Load the sparse matrix - distributed across ranks */
	double t_load_start = now_sec();
	matrix = csc_load_matrix(filepath, mpi_rank, mpi_size);
	double t_load_end = now_sec();
	
	if (!matrix) {
		if (mpi_rank == 0) {
			print_error(__func__, "Failed to load matrix", 0);
		}
		MPI_Finalize();
		return 1;
	}

	/* Initialize benchmarking structure (only rank 0 prints) */
	benchmark = benchmark_init(filepath, n_trials, matrix);
	if (!benchmark) {
		csc_free_matrix(matrix);
		MPI_Finalize();
		return 1;
	}
	benchmark->matrix_info.load_time_s = t_load_end - t_load_start;
	benchmark->benchmark_info.mpi_ranks = mpi_size;
	benchmark->benchmark_info.mpi_rank = mpi_rank;
	
	/* For multi-rank: gather global matrix dimensions */
	if (mpi_size > 1) {
		/* Gather total columns across all ranks */
		unsigned int local_cols = matrix->ncols;
		unsigned int total_cols = 0;
		MPI_Reduce(&local_cols, &total_cols, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
		
		/* Gather total nnz across all ranks */
		unsigned int local_nnz = matrix->nnz;
		unsigned int total_nnz = 0;
		MPI_Reduce(&local_nnz, &total_nnz, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (mpi_rank == 0) {
			benchmark->matrix_info.cols = total_cols;
			benchmark->matrix_info.nnz = total_nnz;
			/* rows stays as is (total vertices) */
		}
	}

	/* Actually run the benchmark */
	ret = benchmark_cc_mpi(matrix, benchmark, mpi_rank, mpi_size);

	/* Only rank 0 prints results */
	if (mpi_rank == 0) {
		benchmark_print(benchmark);
	}

	/* Cleanup */
	benchmark_free(benchmark);
	csc_free_matrix(matrix);
	
	MPI_Finalize();
	return ret;
}
