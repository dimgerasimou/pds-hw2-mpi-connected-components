#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "matrix.h"

/**
 * @struct Statistics
 * @brief Statistical summary of benchmark timing results
 * 
 * Provides comprehensive timing statistics computed across multiple
 * benchmark trials.
 */
typedef struct {
	double mean_time_s;    /**< Mean execution time in seconds */
	double std_dev_s;      /**< Standard deviation of execution time in seconds */
	double median_time_s;  /**< Median execution time in seconds */
	double min_time_s;     /**< Minimum execution time in seconds */
	double max_time_s;     /**< Maximum execution time in seconds */
} Statistics;

/**
 * @struct Result
 * @brief Complete benchmark result for a single algorithm
 * 
 * Contains all measured and computed metrics for one algorithm implementation.
 */
typedef struct {
	unsigned int connected_components;   /**< Number of connected components found */
	Statistics stats;                    /**< Timing statistics */
	double throughput_edges_per_sec;     /**< Processing throughput in edges per second */
} Result;

/**
 * @struct SystemInfo
 * @brief System information captured during benchmark execution
 * 
 * Contains details about the hardware and system configuration where
 * the benchmark was executed.
 */
typedef struct {
	char timestamp[32];    /**< ISO 8601 timestamp of benchmark execution */
	char cpu_info[128];    /**< CPU model and specifications */
	double ram_mb;         /**< Total RAM in megabytes */
	double swap_mb;        /**< Total swap space in megabytes */
} SystemInfo;

/**
 * @struct MatrixInfo
 * @brief Information about the input matrix/graph
 * 
 * Describes the sparse matrix used as input for the connected
 * components algorithm, including its dimensions and sparsity.
 */
typedef struct {
	char path[256];     /**< File path to the matrix */
	unsigned int rows;  /**< Number of rows in the matrix */
	unsigned int cols;  /**< Number of columns in the matrix */
	unsigned int nnz;   /**< Number of non-zero elements (edges in graph) */
} MatrixInfo;

/**
 * @struct BenchmarkInfo
 * @brief Benchmark execution parameters
 * 
 * Contains the configuration parameters used for running the benchmark.
 */
typedef struct {
	unsigned int threads;  /**< Number of threads used for parallel execution */
	// unsigned int mpi_nodes;
	// unsigned int mpi_rank_per_node;
	unsigned int trials;   /**< Number of benchmark trials performed */
} BenchmarkInfo;

/**
 * @brief Holds benchmark results and metadata.
 */
typedef struct {
	double *times;                /**< Array of trial execution times in seconds. */
	SystemInfo sys_info;          /**< System information */
	MatrixInfo matrix_info;       /**< Matrix/graph information */
	BenchmarkInfo benchmark_info; /**< Benchmark parameters */
	Result result;                /**< Algorithm result */
} Benchmark;

/**
 * @brief Initializes a benchmark structure.
 *
 * Allocates and populates a new Benchmark instance for the specified
 * algorithm and dataset. Also allocates memory for timing results.
 *
 * @param name Name of the algorithm being benchmarked.
 * @param filepath Path to the dataset file.
 * @param n_trials Number of trials to run.
 * @param mat Pointer to the CSCBinaryMatrix used as input.
 *
 * @return Pointer to a newly allocated Benchmark structure, or `NULL` on failure.
 */
Benchmark* benchmark_init(const char *filepath, const unsigned int n_trials, const CSCBinaryMatrix *mat);

/**
 * @brief Frees a Benchmark structure and all allocated resources.
 *
 * @param b Pointer to the Benchmark structure to free. Safe to call with NULL.
 */
void benchmark_free(Benchmark *b);

/**
 * @brief Runs a connected components benchmark.
 *
 * Executes the provided connected components function multiple times,
 * measuring execution time per trial and verifying consistency of results.
 *
 * @param m Input CSCBinaryMatrix.
 * @param b Benchmark object containing configuration and result storage.
 *
 * @return
 * - `0` on success,
 * - `1` on algorithm failure or invalid data,
 * - `2` if results differ between trials.
 */
int benchmark_cc(const CSCBinaryMatrix *m, Benchmark *b);

/**
 * @brief Prints benchmark results in structured JSON format.
 *
 * Outputs benchmark metadata, timing statistics, system information,
 * and matrix properties in JSON form for easy parsing or logging.
 *
 * @param b Pointer to the Benchmark structure with populated data.
 */
void benchmark_print(Benchmark *b);

#endif /* BENCHMARK_H */
