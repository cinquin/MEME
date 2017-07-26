#include <getopt.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "config.h"
#include "motif-in.h"
#include "red-black-tree.h"
#include "string-builder.h"


#define PAGEWIDTH 80

VERBOSE_T verbosity = QUIET_VERBOSE;

/**************************************************************************
 * Holds the program options.
 **************************************************************************/
typedef struct meme_get_motif_options {
  RBTREE_T *ids;
  bool match_all;
  bool match_alternate;
  bool match_id_or_alternate;
  bool reverse_complement;
  char* motif_file;
} OPTIONS_T;

/**************************************************************************
 * Prints a usage message and exits. 
 * If given an error message it prints that first and will exit with
 * return code of EXIT_FAILURE.
 **************************************************************************/
static void usage(char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage: meme-get-motif [options] [<MEME file>]\n"
    "Options:\n"
    "    -id <id>       id of motif to extract from the MEME file\n"
    "    -a             match alternate id instead of id\n"
    "    -ia            match either id or alternate id\n"
    "    -rc            reverse complement motifs (assuming alphabet allows)\n"
    "    -all           get all motifs in the MEME file\n\n"
    "Description:\n"
    "    Extract motif(s) from a MEME formatted file and writes out to standard\n"
    "    output.\n";


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

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {OPT_HELP, OPT_ID, OPT_MATCH_ALT, OPT_MATCH_ID_OR_ALT, OPT_REV_COMP, OPT_ALL};
/**************************************************************************
 * Process the command line arguments and put the results in the
 * options structure.
 **************************************************************************/
void process_arguments(int argc, char **argv, OPTIONS_T *options) {
  struct option arg_opts[] = {
    {"help", no_argument, NULL, OPT_HELP},
    {"id", required_argument, NULL, OPT_ID},
    {"a", no_argument, NULL, OPT_MATCH_ALT},
    {"ia", no_argument, NULL, OPT_MATCH_ID_OR_ALT},
    {"rc", no_argument, NULL, OPT_REV_COMP},
    {"all", no_argument, NULL, OPT_ALL},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set defaults
  memset(options, 0, sizeof(OPTIONS_T));
  options->ids = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
  options->match_all = false;
  options->match_alternate = false;
  options->match_id_or_alternate = false;
  options->reverse_complement = false;
  options->motif_file = "-";

  // parse optional arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", arg_opts, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_HELP:           //-help
        usage(NULL);
        break;
      case OPT_ID:
        rbtree_make(options->ids, optarg, NULL);
        break;
      case OPT_MATCH_ALT:
        options->match_alternate = true;
        break;
      case OPT_MATCH_ID_OR_ALT:
        options->match_id_or_alternate = true;
        break;
      case OPT_ALL:
        options->match_all = true;
        break;
      case OPT_REV_COMP:
        options->reverse_complement = true;
        break;
      case '?': //unrecognised or ambiguous argument
        usage("Unknown or ambiguous option");
        break;
      default: // just in case we forget to handle a option
        die("Unhandled option %d", opt);
    }
  }
  // get the motif file
  if (optind < argc) {
    options->motif_file = argv[optind];
    if (!file_exists(options->motif_file)) {
      usage("Motif file \"%s\" does not exist!", options->motif_file);
    }
  }
}
/**************************************************************************
 * Generate logos for all motifs in a file
 **************************************************************************/
static void generate_minimal_meme(OPTIONS_T *options, FILE *out) {
  int i, j, k, strands, pcol;
  ALPH_T *alph;
  ARRAY_T *bg;
  MREAD_T *mread;
  STR_T *buf;
  buf = str_create(10);
  mread = mread_create(options->motif_file, OPEN_MFILE | SKIP_POST_PROCESSING);
  alph = alph_hold(mread_get_alphabet(mread));
  if (alph == NULL) die("Could not load alphabet from \"%s\"!", options->motif_file);
  strands = mread_get_strands(mread);
  bg = mread_get_background(mread);
  if (bg == NULL) die("Could not load background from \"%s\"!", options->motif_file);
  // output the motif file header
  fputs("MEME version " VERSION "\n\n", out);
  alph_print(alph, true, out);
  fputs("END ALPHABET\n\n", out);
  if (strands == 2) {
    fputs("strands: + -\n\n", out);
  } else if (strands == 1) {
    fputs("strands: +\n\n", out);
  }
  fprintf(out, "Background letter frequencies (from unknown source):\n");
  for (i = 0, pcol = 0; i < alph_size_core(alph); i++) {
    pcol += 8; // start of next printed thing 
    if (pcol >= PAGEWIDTH) {
      pcol = 8;
      fputc('\n', out);
    } else {
      fputc(' ', out);
    }
    fprintf(out, "%c %5.3f", alph_char(alph, i), get_array_item(i, bg));
  }
  fputs("\n\n", out);
  // output the motifs
  while (mread_has_motif(mread)) {
    MOTIF_T *motif;
    MATRIX_T *freqs, *scores;
    char *alt;
    bool print;
    // get the next motif
    motif = mread_next_motif(mread);
    // skip motifs if we don't have the IDs in our set
    print = true;
    if (!options->match_all) {
      if (options->match_id_or_alternate) {
        if (rbtree_find(options->ids, get_motif_id(motif)) == NULL &&
            rbtree_find(options->ids, get_motif_id2(motif)) == NULL) print = false;
      } else if (options->match_alternate) {
        if (rbtree_find(options->ids, get_motif_id2(motif)) == NULL) print = false;
      } else {
        if (rbtree_find(options->ids, get_motif_id(motif)) == NULL) print = false;
      }
    }
    if (print) {
      fprintf(out, "\nMOTIF %s", get_motif_id(motif));
      alt = get_motif_id2(motif);
      if (alt != NULL && alt[0] != '\0') {
        fprintf(out, " %s", alt);
      }
      fputs("\n\n", out);
      if (alph_has_complement(alph) && options->reverse_complement) {
        reverse_complement_motif(motif);
      }
      scores = get_motif_scores(motif);
      // optionally output the log-odds
      if (scores != NULL) {
        fprintf(out, "log-odds matrix: alength= %d w= %d E= %s\n",
            alph_size_core(alph), get_motif_length(motif),
            str_evalue(buf, get_motif_log_evalue(motif), 1));
        for (j = 0; j < get_num_rows(scores); j++) {
          for (k = 0; k < get_num_cols(scores); k++) {
            if (k != 0) fputc(' ', out);
            fprintf(out, "%6.0f", get_matrix_cell(j, k, scores));
          }
          fputc('\n', out);
        }
        fputc('\n', out);
      }
      // output the probabilities
      freqs = get_motif_freqs(motif);
      if (freqs != NULL) {
        fprintf(out, "letter-probability matrix: alength= %d w= %d nsites= %.0f E= %s\n",
            alph_size_core(alph), get_motif_length(motif), get_motif_nsites(motif), 
            str_evalue(buf, get_motif_log_evalue(motif), 1));
        for (j = 0; j < get_num_rows(freqs); j++) {
          for (k = 0; k < get_num_cols(freqs); k++) {
            if (k != 0) fputc(' ', out);
            fprintf(out, "%9.6f", get_matrix_cell(j, k, freqs));
          }
          fputc('\n', out);
        }
      }
      // output the URL
      if (has_motif_url(motif)) {
        fprintf(out, "\nURL %s\n", get_motif_url(motif));
      }
    }
    destroy_motif(motif);
  }
  free_array(bg);
  alph_release(alph);
  mread_destroy(mread);
  str_destroy(buf, false);
}

/**************************************************************************
 * Run the program
 **************************************************************************/
int main(int argc, char** argv) {
  OPTIONS_T options;
  process_arguments(argc, argv, &options);
  generate_minimal_meme(&options, stdout);
  rbtree_destroy(options.ids);
  return EXIT_SUCCESS;
}
