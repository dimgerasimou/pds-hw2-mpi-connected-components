/**
 * @file json.h
 * @brief Minimal JSON printer for benchmark output.
 */

#ifndef JSON_H
#define JSON_H

#include "benchmark.h"

/**
 * @brief Print system information as formatted JSON
 * 
 * @param info Pointer to SystemInfo structure to print
 * @param indent_level Number of spaces to indent the output
 * 
 * @note Output is written to stdout
 */
void print_sys_info(const SystemInfo *info, int indent_level);

/**
 * @brief Print matrix information as formatted JSON
 * 
 * @param info Pointer to MatrixInfo structure to print
 * @param indent_level Number of spaces to indent the output
 * 
 * @note Output is written to stdout
 */
void print_matrix_info(const MatrixInfo *info, int indent_level);

/**
 * @brief Print benchmark parameters as formatted JSON
 * 
 * @param info Pointer to BenchmarkInfo structure to print
 * @param indent_level Number of spaces to indent the output
 * 
 * @note Output is written to stdout
 */
void print_benchmark_info(const BenchmarkInfo *info, int indent_level);

/**
 * @brief Print algorithm result as formatted JSON
 * 
 * @param result Pointer to Result structure to print
 * @param indent_level Number of spaces to indent the output
 * 
 * @note Output is written to stdout
 */
void print_result(const Result *result, int indent_level);

#endif /* JSON_H */
