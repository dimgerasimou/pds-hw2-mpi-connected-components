/**
 * @file mtx_to_bin.c
 * @brief Convert Matrix Market (.mtx) files to fast binary CSC format.
 *
 * Usage: ./mtx2bin input.mtx output.bin
 *
 * Binary format:
 *   - uint32_t nrows
 *   - uint32_t ncols
 *   - size_t   nnz
 *   - uint32_t col_ptr[ncols + 1]
 *   - uint32_t row_idx[nnz]
 */
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

/* ------------------------------------------------------------------------- */
/*                            Helper Functions                               */
/* ------------------------------------------------------------------------- */

static int
mm_skip_comments(FILE *f)
{
	long pos;
	int c;

	while (1) {
		pos = ftell(f);
		c = fgetc(f);
		if (c == '%') {
			while ((c = fgetc(f)) != '\n' && c != EOF);
			continue;
		} else if (isspace(c)) {
			continue;
		}
		fseek(f, pos, SEEK_SET);
		return 0;
	}
}

static void
print_progress(size_t current, size_t total, time_t start_time)
{
	if (current % 100000000 == 0 || current == total) {
		time_t now = time(NULL);
		double elapsed = difftime(now, start_time);
		double rate = current / (elapsed + 0.001);
		double eta = (total - current) / (rate + 0.001);

		fprintf(stderr, "\rProcessed: %zu / %zu edges (%.1f%%) "
		        "Rate: %.2fM edges/s ETA: %.0fs   ",
		        current, total, 100.0 * current / total,
		        rate / 1e6, eta);
		fflush(stderr);
	}
}

/* ------------------------------------------------------------------------- */
/*                            Main Conversion                                */
/* ------------------------------------------------------------------------- */

int
main(int argc, char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s input.mtx output.bin\n", argv[0]);
		return 1;
	}

	const char *infile = argv[1];
	const char *outfile = argv[2];
	time_t start_time = time(NULL);

	/* --- Open input ---------------------------------------------------- */
	FILE *fin = fopen(infile, "r");
	if (!fin) {
		perror("fopen input");
		return 1;
	}

	/* --- Parse header -------------------------------------------------- */
	char format[64], field[64], symmetry[64];

	if (fscanf(fin, "%%%%MatrixMarket matrix %63s %63s %63s",
	           format, field, symmetry) != 3)
	{
		fprintf(stderr, "Error: invalid MatrixMarket header\n");
		fclose(fin);
		return 1;
	}

	int is_coordinate = (strcmp(format, "coordinate") == 0);
	int is_pattern = (strcmp(field, "pattern") == 0);
	int symmetric = (strcmp(symmetry, "symmetric") == 0);

	fprintf(stderr, "Format: %s, Field: %s, Symmetry: %s\n",
	        format, field, symmetry);

	/* --- Parse dimensions ---------------------------------------------- */
	mm_skip_comments(fin);

	size_t nrows, ncols, nnz;
	if (is_coordinate) {
		if (fscanf(fin, "%zu %zu %zu", &nrows, &ncols, &nnz) != 3) {
			fprintf(stderr, "Error: invalid size line\n");
			fclose(fin);
			return 1;
		}
	} else {
		if (fscanf(fin, "%zu %zu", &nrows, &ncols) != 2) {
			fprintf(stderr, "Error: invalid array size\n");
			fclose(fin);
			return 1;
		}
		nnz = nrows * ncols;
	}

	fprintf(stderr, "Matrix: %zu x %zu, NNZ: %zu\n", nrows, ncols, nnz);

	/* --- Allocate COO arrays ------------------------------------------- */
	size_t max_nnz = nnz * (symmetric ? 2 : 1);
	uint32_t *coo_i = malloc(max_nnz * sizeof(uint32_t));
	uint32_t *coo_j = malloc(max_nnz * sizeof(uint32_t));

	if (!coo_i || !coo_j) {
		fprintf(stderr, "Error: malloc failed\n");
		fclose(fin);
		return 1;
	}

	/* --- Read edges ---------------------------------------------------- */
	fprintf(stderr, "Reading edges...\n");
	size_t count = 0;

	if (is_coordinate) {
		for (size_t k = 0; k < nnz; k++) {
			size_t i, j;
			double val = 1.0;

			if (is_pattern) {
				if (fscanf(fin, "%zu %zu", &i, &j) != 2) {
					fprintf(stderr, "Error: bad entry at line %zu\n", k);
					goto cleanup;
				}
			} else {
				if (fscanf(fin, "%zu %zu %lf", &i, &j, &val) != 3) {
					fprintf(stderr, "Error: bad entry at line %zu\n", k);
					goto cleanup;
				}
			}

			if (val != 0.0) {
				coo_i[count] = i - 1;
				coo_j[count] = j - 1;
				count++;

				if (symmetric && i != j) {
					coo_i[count] = j - 1;
					coo_j[count] = i - 1;
					count++;
				}
			}

			print_progress(k + 1, nnz, start_time);
		}
	} else {
		/* array format */
		for (size_t j = 0; j < ncols; j++) {
			for (size_t i = 0; i < nrows; i++) {
				double val;
				if (fscanf(fin, "%lf", &val) != 1) {
					fprintf(stderr, "Error: bad array entry\n");
					goto cleanup;
				}

				if (val != 0.0) {
					coo_i[count] = i;
					coo_j[count] = j;
					count++;
				}

				print_progress(i + j * nrows + 1, nnz, start_time);
			}
		}
	}

	fclose(fin);
	fprintf(stderr, "\nRead %zu edges\n", count);

	/* --- Convert COO to CSC -------------------------------------------- */
	fprintf(stderr, "Converting to CSC format...\n");

	uint32_t *col_ptr = calloc(ncols + 1, sizeof(uint32_t));
	uint32_t *row_idx = malloc(count * sizeof(uint32_t));

	if (!col_ptr || !row_idx) {
		fprintf(stderr, "Error: malloc failed for CSC\n");
		goto cleanup;
	}

	/* count entries per column */
	for (size_t k = 0; k < count; k++)
		col_ptr[coo_j[k] + 1]++;

	/* prefix sum */
	for (size_t j = 0; j < ncols; j++)
		col_ptr[j + 1] += col_ptr[j];

	/* fill row indices */
	uint32_t *col_fill = calloc(ncols, sizeof(uint32_t));
	if (!col_fill) {
		fprintf(stderr, "Error: malloc failed\n");
		free(col_ptr);
		free(row_idx);
		goto cleanup;
	}

	for (size_t k = 0; k < count; k++) {
		uint32_t j = coo_j[k];
		uint32_t dest = col_ptr[j] + col_fill[j];
		row_idx[dest] = coo_i[k];
		col_fill[j]++;
	}

	free(col_fill);
	free(coo_i);
	free(coo_j);

	/* --- Write binary file --------------------------------------------- */
	fprintf(stderr, "Writing binary file...\n");

	FILE *fout = fopen(outfile, "wb");
	if (!fout) {
		perror("fopen output");
		free(col_ptr);
		free(row_idx);
		return 1;
	}

	uint32_t nrows_u32 = nrows;
	uint32_t ncols_u32 = ncols;

	/* write header */
	if (fwrite(&nrows_u32, sizeof(uint32_t), 1, fout) != 1 ||
	    fwrite(&ncols_u32, sizeof(uint32_t), 1, fout) != 1 ||
	    fwrite(&count, sizeof(size_t), 1, fout) != 1)
	{
		perror("fwrite header");
		fclose(fout);
		free(col_ptr);
		free(row_idx);
		return 1;
	}

	/* write col_ptr */
	if (fwrite(col_ptr, sizeof(uint32_t), ncols + 1, fout) != ncols + 1) {
		perror("fwrite col_ptr");
		fclose(fout);
		free(col_ptr);
		free(row_idx);
		return 1;
	}

	/* write row_idx in chunks for better performance */
	size_t chunk_size = 100000000; /* 100M entries at a time */
	size_t written = 0;
	while (written < count) {
		size_t to_write = (count - written < chunk_size) ?
		                  (count - written) : chunk_size;

		if (fwrite(&row_idx[written], sizeof(uint32_t), to_write, fout) != to_write) {
			perror("fwrite row_idx");
			fclose(fout);
			free(col_ptr);
			free(row_idx);
			return 1;
		}

		written += to_write;
		fprintf(stderr, "\rWritten: %zu / %zu edges (%.1f%%)",
		        written, count, 100.0 * written / count);
		fflush(stderr);
	}

	fclose(fout);
	free(col_ptr);
	free(row_idx);

	time_t end_time = time(NULL);
	double total_time = difftime(end_time, start_time);

	fprintf(stderr, "\n\nConversion complete!\n");
	fprintf(stderr, "Total time: %.0f seconds\n", total_time);
	fprintf(stderr, "Output file: %s\n", outfile);

	return 0;

cleanup:
	fclose(fin);
	free(coo_i);
	free(coo_j);
	return 1;
}
