#include <getopt.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include "mtwist.h"

#include "alphabet.h"
#include "array.h"
#include "config.h"
#include "matrix.h"
#include "motif.h"
#include "motif-in.h"


VERBOSE_T verbosity = QUIET_VERBOSE;

typedef struct {
  mt_state prng;
  uint32_t seed;
  ARRAYLST_T *in_files;
  char *out_file;
} OPTIONS_T;

/**************************************************************************
 * Prints a usage message and exits. 
 * If given an error message it prints that first and will exit with
 * return code of EXIT_FAILURE.
 **************************************************************************/
static void usage(char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n" 
    "    motif-shuffle-columns [options] <motif db>+\n\n"
    "Options:\n"
    "    -o <file name>   output file name\n"
    "    -seed <seed>     pseudo-random number generator seed\n"
    "    -help            print this usage message\n\n"
    "Description:\n"
    "    Creates motifs from the shuffled columns of the input motifs.\n";

  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fputs(usage, stderr);
    fflush(stderr);
  } else {
    puts(usage);
  }
  if (format) exit(EXIT_FAILURE);
  exit(EXIT_SUCCESS);
}

/**************************************************************************
 * Process the command line arguments and put the results in the
 * options structure.
 **************************************************************************/
void process_arguments(int argc, char **argv, OPTIONS_T *options) {
  bool bad_argument = false;

  struct option msc_options[] = {
    {"seed", required_argument, NULL, 's'},
    {"o", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set defaults
  options->seed = mts_seed(&(options->prng));
  options->out_file = NULL;
  options->in_files = arraylst_create();

  // parse optional arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", msc_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case 's': // pseudo-random number generator seed
        if (sscanf(optarg, "%" SCNu32, &(options->seed)) == 1) {
          mts_seed32new(&(options->prng), options->seed);
        } else {
          usage("Seed \"%s\" could not be interpreted as an unsigned 32bit integer", optarg);
        }
        break;
      case 'o': //-o <file name>
        options->out_file = optarg;
        break;
      case 'h': //-help
        usage(NULL);
        break;
      case '?': //unrecognised or ambiguous argument
        bad_argument = TRUE;
    }
  }
  if (bad_argument) usage("One or more unknown or ambiguous options were supplied.");

  // get the motif file
  if (optind >= argc) usage("No motif db specified.");
  while (optind < argc) arraylst_add(argv[optind++], options->in_files);
}

/**************************************************************************
 * Create a number in the range 0 to n (inclusive).
 **************************************************************************/
static inline uint32_t rand_uint32(mt_state* prng, const uint32_t n) {
  uint64_t d, threshold;
  uint32_t r;
  assert(n > 0);
  d = (uint64_t)n + 1;
  threshold = (uint64_t)UINT32_MAX + 1;
  threshold = threshold - (threshold % d);
  assert(threshold > 0);
  do {
    r = mts_lrand(prng);
  } while (r >= threshold);
  return r % d;
}

/**************************************************************************
 * Run the program
 **************************************************************************/
int main(int argc, char** argv) {
  OPTIONS_T options;
  ALPH_T *alph;
  ARRAY_T *bg;
  int i, j, k, len, offset, cols_alloc, cols_used, lens_alloc, lens_used;
  int *lens;
  double **cols;
  FILE *out;
  // setup memory
  alph = NULL;
  bg = NULL;
  cols_alloc = 512;
  cols_used = 0;
  cols = mm_calloc(cols_alloc, sizeof(double*));
  lens_alloc = 64;
  lens_used = 0;
  lens = mm_calloc(lens_alloc, sizeof(int));
  // process args
  process_arguments(argc, argv, &options);
  // open output file
  if (options.out_file) {
    if ((out = fopen(options.out_file, "w")) == NULL) {
      die("Failed to open \"%s\" for writing", options.out_file);
    }
  } else {
    out = stdout;
  }
  // read in all motif columns
  for (i = 0; i < arraylst_size(options.in_files); i++) {
    MREAD_T *mread;
    char *file;
    file = (char*)arraylst_get(i, options.in_files);
    mread = mread_create(file, OPEN_MFILE);
    if (alph == NULL) {
      alph = alph_hold(mread_get_alphabet(mread));
    } else if (!alph_equal(alph, mread_get_alphabet(mread))) {
      die("alphabet of %s is not %s", file, alph_name(alph));
    }
    if (bg == NULL) bg = mread_get_background(mread);
    while (mread_has_motif(mread)) {
      MOTIF_T *motif;
      MATRIX_T *freqs;
      motif = mread_next_motif(mread);
      len = get_motif_length(motif);
      if (lens_used >= lens_alloc) {
        lens_alloc *= 2;
        lens = mm_realloc(lens, lens_alloc * sizeof(int));
        memset(lens+lens_used, 0, (lens_alloc - lens_used) * sizeof(int));
      }
      lens[lens_used++] = len;
      if ((cols_used + len) >= cols_alloc) {
        do {
          cols_alloc *= 2;
        } while ((cols_used + len) >= cols_alloc);
        cols = mm_realloc(cols, cols_alloc * sizeof(double*));
        memset(cols+cols_used, 0, (cols_alloc - cols_used) * sizeof(double*));
      }
      freqs = get_motif_freqs(motif);
      for (j = 0; j < len; j++) {
        cols[cols_used] = mm_calloc(alph_size_core(alph), sizeof(double));
        for (k = 0; k < alph_size_core(alph); k++) {
          cols[cols_used][k] = get_matrix_cell(j, k, freqs);
        }
        cols_used++;
      }
    }
  }
  // shuffle all columns
  for (i = cols_used - 1; i >= 1; i--) {
    double *tmp;
    // generate number 0 <= j <= i.
    j = rand_uint32(&(options.prng), i);
    // swap columns
    tmp = cols[i];
    cols[i] = cols[j];
    cols[j] = tmp;
  }
  // write motifs
  if (cols_used > 0) {
    fprintf(out, "# PRNG seed %" PRIu32 "\n", options.seed);
    fprintf(out, "MEME version %s\n\n", VERSION);
    alph_print(alph, true, out);
    fputs("END ALPHABET\n\n", out);
    fputs("Background letter frequencies\n", out);
    for (i = 0; i < alph_size_core(alph); i++) {
      if (i != 0) fputs(" ", out);
      fprintf(out, "%c %.3f", alph_char(alph, i), get_array_item(i, bg)); 
    }
    fputs("\n\n", out);
    offset = 0;
    for (i = 0; i < lens_used; i++) {
      len = lens[i];
      fprintf(out, "MOTIF %d random\n", i + 1);
      fprintf(out, "letter-probability matrix: alength= %d w= %d\n", alph_size_core(alph), lens[i]);
      for (j = 0; j < len; j++) {
        double *col;
        col = cols[offset + j];
        for (k = 0; k < alph_size_core(alph); k++) {
          if (k != 0) fputs(" ", out);
          fprintf(out, "%.6f", col[k]);
        }
        fputs("\n", out);
      }
      fputs("\n", out);
      offset += len;
    }
  }
  // clean up
  for (i = 0; i < cols_used; i++) free(cols[i]);
  free(cols);
  free(lens);
  return EXIT_SUCCESS;
}

