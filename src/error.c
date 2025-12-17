/**
 * @file error.c
 * @brief Implementation of error handling utilities.
 */

#include <stdio.h>
#include <string.h>

#include "error.h"

/**
 * @brief Global variable storing the program name for error messages.
 *
 * This variable should be initialized by calling `set_program_name()`
 * before any calls to `print_error()`.
 */
extern const char *program_name;

/**
 * @copydoc set_program_name()
 */
void
set_program_name(const char *argv0)
{
	if (argv0) {
		const char *slash = strrchr(argv0, '/');
		program_name = slash ? slash + 1 : argv0;
	}
}

/**
 * @copydoc print_error()
 */
void
print_error(const char *func, const char *msg, int err)
{
	if (err)
		fprintf(stderr, "%s: %s: %s: %s\n", program_name, func, msg, strerror(err));
	else
		fprintf(stderr, "%s: %s: %s\n", program_name, func, msg);
}
