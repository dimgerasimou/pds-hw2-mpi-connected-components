/**
 * @file json.c
 * @brief Minimal JSON printer for benchmark output
 */

#include <stdio.h>
#include <string.h>

#include "json.h"

/**
 * @brief Print system information as formatted JSON.
 */
void
print_sys_info(const SystemInfo *info, int indent_level)
{
	printf("%*s\"sys_info\": {\n", indent_level, "");
	printf("%*s\"timestamp\": \"%s\",\n", indent_level + 2, "", info->timestamp);
	printf("%*s\"cpu_info\": \"%s\",\n", indent_level + 2, "", info->cpu_info);
	printf("%*s\"ram_mb\": %.2f,\n", indent_level + 2, "", info->ram_mb);
	printf("%*s\"swap_mb\": %.2f\n", indent_level + 2, "", info->swap_mb);
	printf("%*s}", indent_level, "");
}

/**
 * @brief Print matrix information as formatted JSON.
 */
void
print_matrix_info(const MatrixInfo *info, int indent_level)
{
	printf("%*s\"matrix_info\": {\n", indent_level, "");
	printf("%*s\"path\": \"%s\",\n", indent_level + 2, "", info->path);
	printf("%*s\"rows\": %u,\n", indent_level + 2, "", info->rows);
	printf("%*s\"cols\": %u,\n", indent_level + 2, "", info->cols);
	printf("%*s\"nnz\": %u\n", indent_level + 2, "", info->nnz);
	printf("%*s}", indent_level, "");
}

/**
 * @brief Print benchmark parameters as formatted JSON.
 */
void
print_benchmark_info(const BenchmarkInfo *info, int indent_level)
{
	printf("%*s\"benchmark_info\": {\n", indent_level, "");
	printf("%*s\"threads\": %u,\n", indent_level + 2, "", info->threads);
	printf("%*s\"trials\": %u\n", indent_level + 2, "", info->trials);
	printf("%*s}", indent_level, "");
}

/**
 * @brief Print algorithm result as formatted JSON.
 */
void
print_result(const Result *result, int indent_level)
{
	printf("%*s{\n", indent_level, "");
	printf("%*s\"connected_components\": %u,\n", indent_level + 2, "", result->connected_components);
	printf("%*s\"statistics\": {\n", indent_level + 2, "");
	printf("%*s\"mean_time_s\": %.6f,\n", indent_level + 4, "", result->stats.mean_time_s);
	printf("%*s\"std_dev_s\": %.6f,\n", indent_level + 4, "", result->stats.std_dev_s);
	printf("%*s\"median_time_s\": %.6f,\n", indent_level + 4, "", result->stats.median_time_s);
	printf("%*s\"min_time_s\": %.6f,\n", indent_level + 4, "", result->stats.min_time_s);
	printf("%*s\"max_time_s\": %.6f\n", indent_level + 4, "", result->stats.max_time_s);
	printf("%*s},\n", indent_level + 2, "");
	printf("%*s\"throughput_edges_per_sec\": %.2f,\n", indent_level + 2, "", result->throughput_edges_per_sec);
	printf("%*s}", indent_level, "");
}
