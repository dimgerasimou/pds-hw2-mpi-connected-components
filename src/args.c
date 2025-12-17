/**
 * @file args.c
 * @brief Command-line argument parsing implementation.
 *
 * Provides functions to parse program arguments that specify
 * the number of threads, number of trials, and input file path.
 */

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "args.h"
#include "error.h"

extern const char *program_name;

/**
 * @brief Checks if a string represents an unsigned integer.
 *
 * @param[in] s String to check.
 * @return 1 if @p s is a valid unsigned integer, 0 otherwise.
 */
static int
isuint(const char *s)
{
	if (!s || s[0] == '\0')
		return 0;

	for (const char *ptr = s; *ptr != '\0'; ptr++) {
		if (!(*ptr >= '0' && *ptr <= '9'))
			return 0;
	}

	return 1;
}

/**
 * @brief Prints program usage instructions to stdout.
 *
 * Displays the valid command-line options and their expected arguments.
 */
static void
usage(void) {
	fprintf(stderr,
		"Usage: %s [OPTIONS] <matrix_file>\n\n"
		"Options:\n"
		"  -n <trials>        Number of benchmark trials (default: 3)\n"
		"  -h                 Show this help message and exit\n\n"
		"Arguments:\n"
		"  matrix_file Path to the input matrix file (Matrix Market format)\n\n"
		"Example:\n"
		"  %s -n 10 ./data/matrix.mat\n",
		program_name, program_name
	);
}

/**
 * @copydoc parseargs()
 */
int
parseargs(int argc, char *argv[],
          unsigned int *n_trials,
          char **filepath)
{
	*n_trials = 3;
	*filepath = NULL;

	opterr = 0;

	int opt;
	while ((opt = getopt(argc, argv, "+n:h")) != -1) {
		switch (opt) {
		case 'n': {
			if (!optarg || !isuint(optarg)) {
				char err[128];
				snprintf(err, sizeof(err), "invalid or missing argument for -%c", opt);
				print_error(__func__, err, 0);
				usage();
				return 1;
			}
			int val = atoi(optarg);
			if (!val) {
				char err[128];
				snprintf(err, sizeof(err), "%s must be > 0", (opt == 't') ? "threads" : "trials");
				print_error(__func__, err, 0);
				usage();
				return 1;
			}
			*n_trials = val;
			break;
		}
		case 'h':
			usage();
			return -1;
		
		case '?':
		default: {
			char err[128];
			if (optopt == 't' || optopt == 'n')
				snprintf(err, sizeof(err), "missing argument for -%c", optopt);
			else
				snprintf(err, sizeof(err), "unknown option '-%c'", optopt ? optopt : '?');
			print_error(__func__, err, 0);
			usage();
			return 1;
		}
		}
	}

	if (optind < argc) {
		*filepath = argv[optind];
		if (access(*filepath, R_OK) != 0) {
			char err[256];
			snprintf(err, sizeof(err), "cannot access file: \"%s\"", *filepath);
			print_error(__func__, err, errno);
			usage();
			return 1;
		}
	} else {
		print_error(__func__, "no input file specified", 0);
		usage();
		return 1;
	}

	return 0;
}
