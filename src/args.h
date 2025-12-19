/**
 * @file args.h
 * @brief Command-line argument parsing interface.
 *
 * This header declares the function used to parse command-line arguments
 * for configuring the program's execution parameters such as
 * number of trials and input file path.
 */

#ifndef ARGS_H
#define ARGS_H

/**
 * @brief Parses command-line arguments.
 *
 * Supported options:
 *   -n <trials>    Number of trials (default: 3)
 *   -h             Show usage and exit
 *
 * @param[in]  argc     Argument count.
 * @param[in]  argv     Argument vector.
 * @param[out] n_trials Output: number of trials.
 * @param[out] filepath Output: path to matrix file.
 *
 * @return 0 on success, -1 if help requested, 1 on error.
 */
int parseargs(int argc, char *argv[], unsigned int *n_trials, char **filepath);

#endif /* ARGS_H */
