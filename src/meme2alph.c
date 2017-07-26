
#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>

#include "alphabet.h"
#include "motif-in.h"
#include "utils.h"

VERBOSE_T verbosity = QUIET_VERBOSE;

typedef struct motif_get_alph_options {
  char *motifs_file;
  char *alphabet_file;
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
    "    meme2alph [options] <motifs file> [<alphabet file>]\n"
    "Options:\n"
    "    -help            print this usage message\n"
    "Description:\n"
    "    Extracts the motif alphabet and prints it to standard output or the <alphabet file>.\n";

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
  BOOLEAN_T bad_argument = FALSE;

  struct option meme2alph_options[] = {
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set defaults
  options->motifs_file = NULL;
  options->alphabet_file = "-";

  // parse optional arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", meme2alph_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case 'h':           //-help
        usage(NULL);
        break;
      case '?':           //unrecognised or ambiguous argument
        bad_argument = TRUE;
    }
  }
  if (bad_argument) usage("One or more unknown or ambiguous options were supplied.");

  // get the motif file
  if (optind >= argc) usage("No motif file!");
  options->motifs_file = argv[optind++];
  if (!file_exists(options->motifs_file)) 
    usage("Motif file \"%s\" does not exist!", options->motifs_file);

  // get the output directory
  if (optind < argc) {
    options->alphabet_file = argv[optind++];
  }
}

/**************************************************************************
 * Load the alphabet from the motif file
 **************************************************************************/
ALPH_T* load_alphabet(OPTIONS_T *options) {
  MREAD_T *mread;
  ALPH_T *alph;
  mread = mread_create(options->motifs_file, OPEN_MFILE);
  alph = alph_hold(mread_get_alphabet(mread));
  mread_destroy(mread);
  return alph;
}

/**************************************************************************
 * Print the alphabet to the output file
 **************************************************************************/
void output_alphabet(OPTIONS_T *options, ALPH_T *alph) {
  FILE *fh;
  if (strcmp("-", options->alphabet_file) == 0) {
    fh = stdout;
  } else {
    fh = fopen(options->alphabet_file, "w");
    if (fh == NULL) {
      die("Unable to open alphabet file \"%s\" for writing", options->alphabet_file);
    }
  }
  alph_print(alph, true, fh);
  if (fh != stdout) fclose(fh);
}

/**************************************************************************
 * Run the program
 **************************************************************************/
int main(int argc, char** argv) {
  OPTIONS_T options;
  ALPH_T *alph;
  process_arguments(argc, argv, &options);
  alph = load_alphabet(&options);
  if (alph != NULL) {
    output_alphabet(&options, alph);
    alph_release(alph);
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
