#include <getopt.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "config.h"
#include "motif-in.h"
#include "red-black-tree.h"
#include "regex-utils.h"
#include "string-builder.h"


#define PAGEWIDTH 80

VERBOSE_T verbosity = QUIET_VERBOSE;

/**************************************************************************
 * Holds the program options.
 **************************************************************************/
typedef struct meme2meme_options {
  bool use_consensus;
  bool use_numbers;
  bool print_logodds;
  bool replace_url;
  bool xalph;
  char* bg_filename;
  char* url_pattern;
  char** motif_files;
  int motif_file_count;
} OPTIONS_T;

/**************************************************************************
 * Holds a list of MOTIF_T.
 **************************************************************************/
typedef struct motifs {
  MOTIF_T** list;
  int count;
} MOTIFS_T;

/**************************************************************************
 * Frees the memory used by a MOTIFS_T including all the referenced MOTIF_T.
 **************************************************************************/
static void destroy_motifs(void* value) {
  MOTIFS_T *motifs;
  int i;
  motifs = (MOTIFS_T*)value;
  for (i = 0; i < motifs->count; i++) {
    destroy_motif(motifs->list[i]);
  }
  free(motifs->list);
  free(motifs);
}

/**************************************************************************
 * Prints a usage message and exits. 
 * If given an error message it prints that first and will exit with
 * return code of EXIT_FAILURE.
 **************************************************************************/
static void usage(char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n" 
    "    meme2meme [options] <meme file>+\n"
    "Options:\n"
    "    -consensus     numeric names are swapped for an IUPAC\n"
    "                     consensus; default: use existing names\n"
    "    -numbers       use numbers instead of strings for motif names;\n"
    "                     default: use existing ID\n"
    "    -bg <file>     file with background frequencies of letters;\n"
    "                     default: use first file background\n"
    "    -logodds       print log-odds matrix as well as frequency matrix;\n"
    "                     default: frequency matrix only\n"
    "    -url <website> website for the motif if it doesn't have one\n"
    "                     already; The motif name is substituted for\n"
    "                     MOTIF_NAME; default: use existing url\n"
    "    -forceurl      Existing urls are ignored\n"
    "    -xalph         Convert all motifs to use the same alphabet as\n"
    "                     specified in the first file which must be a superset;\n"
    "                     default: all alphabets must be identical\n"
    "Description:\n"
    "    Takes meme motifs in many forms and writes out a single database in\n"
    "    minimal meme format to standard output.\n";


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
  int i;

  struct option arg_opts[] = {
    {"help", no_argument, NULL, 'h'},
    {"consensus", no_argument, NULL, 'c'},
    {"numbers", no_argument, NULL, 'n'},
    {"logodds", no_argument, NULL, 'l'},
    {"forceurl", no_argument, NULL, 'f'},
    {"bg", required_argument, NULL, 'b'},
    {"url", required_argument, NULL, 'u'},
    {"xalph", no_argument, NULL, 'x'},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set defaults
  memset(options, 0, sizeof(OPTIONS_T));
  options->use_consensus = false;
  options->use_numbers = false;
  options->print_logodds = false;
  options->replace_url = false;
  options->bg_filename = NULL;
  options->url_pattern = NULL;
  options->motif_files = NULL;
  options->motif_file_count = 0;

  // parse optional arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", arg_opts, NULL);
    if (opt == -1) break;
    switch (opt) {
      case 'h':           //-help
        usage(NULL);
        break;
      case 'c':
        options->use_consensus = true;
        break;
      case 'n':
        options->use_numbers = true;
        break;
      case 'l':
        options->print_logodds = true;
        break;
      case 'f':
        options->replace_url = true;
        break;
      case 'b':
        options->bg_filename = optarg;
        break;
      case 'u':
        options->url_pattern = optarg;
        break;
      case 'x':
        options->xalph = true;
        break;
      case '?':           //unrecognised or ambiguous argument
        bad_argument = TRUE;
    }
  }
  if (bad_argument) usage("One or more unknown or ambiguous options were supplied.");
  // check background exists
  if (options->bg_filename && !file_exists(options->bg_filename)) {
    usage("Background file \"%s\" does not exist!", options->bg_filename); 
  }
  if (options->replace_url && options->url_pattern == NULL) {
    usage("The URL can only be replaced when a default is provided.");
  }
  // get the motif files
  if (optind >= argc) usage("No motif file!");
  options->motif_files = argv+optind;
  options->motif_file_count = argc - optind;
  for (i = 0; i < options->motif_file_count; i++) {
    if (!file_exists(options->motif_files[i])) {
      usage("Motif file \"%s\" does not exist!", options->motif_files[i]);
    }
  }
}

/**************************************************************************
 * Print out the name of a symbol
 **************************************************************************/
static void print_name(FILE *out, const char *name) {
  const char *c;
  fputc('"', out);
  for (c = name; *c != '\0'; c++) {
    switch (*c) {
      case '"': fputs("\\\"", out); break;
      case '/': fputs("\\/", out); break;
      case '\\': fputs("\\\\", out); break;
      default: fputc(*c, out);
    }
  }
  fputc('"', out);
}

/**************************************************************************
 * Splits the URL pattern into parts surrounding the MOTIF_NAME token
 **************************************************************************/
static char** split_url_pattern(char *url_pattern) {
  char *start, *location, *part;
  char **parts;
  int count, len;
  count = 1;
  parts = mm_malloc(sizeof(char*));
  parts[0] = NULL;
  start = url_pattern;
  while (true) {
    location = strstr(start, "MOTIF_NAME");
    if (location != NULL) {
      len = location - start;
    } else {
      len = strlen(start);
    }
    part = mm_malloc(sizeof(char) * (len + 1));
    strncpy(part, start, len);
    part[len] = '\0';
    count++;
    parts = mm_realloc(parts, sizeof(char*) * count);
    parts[count - 2] = part;
    parts[count - 1] = NULL;
    if (location == NULL) return parts;
    start = location+10; // skip over MOTIF_NAME
  }
}

/**************************************************************************
 * Cluster together items that say the same thing in different cases
 * but don't treat them as the same.
 **************************************************************************/
static int id_cmp(const void* p1, const void* p2) {
  int diff;
  diff = strcasecmp((const char*)p1, (const char*)p2);
  if (diff != 0) return diff;
  return strcmp((const char*)p1, (const char*)p2);
}

/**************************************************************************
 * Generate logos for all motifs in a file
 **************************************************************************/
static void generate_minimal_meme(OPTIONS_T *options, FILE *out) {
  int i, j, k, idx, strands, pcol;
  ALPH_T *alph;
  ARRAY_T *bg;
  MREAD_T *mread;
  RBNODE_T *node;
  RBTREE_T *all_motifs, *all_names;
  STR_T *buf, *id_buf;
  regex_t int_re;
  char **url_parts;
  if (options->motif_file_count <= 0) return;
  url_parts = NULL;
  if (options->url_pattern != NULL) {
    url_parts = split_url_pattern(options->url_pattern);
  }
  regcomp_or_die("integer pattern", &int_re, "^[0-9]+$", REG_EXTENDED);
  buf = str_create(10);
  id_buf = str_create(10);
  all_motifs = rbtree_create(id_cmp, rbtree_strcpy, free, NULL, destroy_motifs);
  all_names = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, NULL);
  // stop the compiler from complaining about uninitilized values
  // (they're all set by the first loop iteration)
  alph = NULL; bg = NULL; strands = 0;
  // loop over all the motif files
  for (i = 0; i < options->motif_file_count; i++) {
    // open motif file
    mread = mread_create(options->motif_files[i], OPEN_MFILE);
    if (i == 0) {
      // get the alphabet and background
      alph = alph_hold(mread_get_alphabet(mread));
      if (alph == NULL) die("Could not load alphabet from \"%s\"!", options->motif_files[i]);
      strands = mread_get_strands(mread);
      if (options->bg_filename) {
        bg = get_file_frequencies(alph, options->bg_filename);
      } else {
        bg = mread_get_background(mread);
        if (bg == NULL) die("Could not load background from \"%s\"!", options->motif_files[i]);
      }
    } else {
      // check that the alphabet matches the first file
      if (mread_has_motif(mread) && !alph_equal(alph, mread_get_alphabet(mread))) {
        if (options->xalph) {
          switch(alph_core_subset(mread_get_alphabet(mread), alph)) {
            case 0:
              die("The motifs in \"%s\" have the %s alphabet which is not "
                  "a subset of the %s alphabet.", options->motif_files[i], 
                  alph_name(mread_get_alphabet(mread)), alph_name(alph));
              break;
            case -1:
              fprintf(stderr, "Warning: the alphabet expansion from %s to %s"
                  " requires changing complementing rules.\n", 
                  alph_name(mread_get_alphabet(mread)), alph_name(alph));
              break;
          }
          // since we're not applying pseudocounts the background doesn't matter
          // so we use the uniform background
          mread_set_conversion(mread, alph, NULL);
        } else {
          die("The motifs in \"%s\" have the %s alphabet which is not "
              "the same as the expected %s alphabet.", options->motif_files[i], 
              alph_name(mread_get_alphabet(mread)), alph_name(alph));
        }
      }
    }
    // load the motifs
    while (mread_has_motif(mread)) {
      MOTIF_T *motif;
      MOTIFS_T *motifs;
      char *motif_id;
      bool new_id;
      motif = mread_next_motif(mread);
      motif_id = get_motif_id(motif);
      // rename numbered motifs if the "-consensus" option is set
      if (options->use_consensus && regexec_or_die("test num id", &int_re, motif_id, 0, NULL, 0)) {
        str_clear(id_buf);
        motif2consensus(motif, id_buf, false);
        motif_id = str_internal(id_buf);
      }
      node = rbtree_lookup(all_motifs, motif_id, true, &new_id);
      if (new_id) {
        motifs = mm_malloc(sizeof(MOTIFS_T));
        motifs->list = NULL;
        motifs->count = 0;
        rbtree_set(all_motifs, node, motifs);
      } else {
        motifs = (MOTIFS_T*)rbtree_value(node);
      }
      motifs->count += 1;
      motifs->list = mm_realloc(motifs->list, sizeof(MOTIF_T*) * motifs->count);
      motifs->list[motifs->count - 1] = motif;
    }
    mread_destroy(mread);
  }
  // output the motif file header
  fputs("MEME version " VERSION "\n\n", out);
  if (alph_is_builtin_dna(alph)) {
    fputs("ALPHABET= ACGT\n\n", out);
  } else if (alph_is_builtin_rna(alph)) {
    fputs("ALPHABET= ACGU\n\n", out);
  } else if (alph_is_builtin_protein(alph)) {
    fputs("ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n", out);
  } else {
    alph_print(alph, true, out);
    fputs("END ALPHABET\n\n", out);
  }
  if (strands == 2) {
    fputs("strands: + -\n\n", out);
  } else if (strands == 1) {
    fputs("strands: +\n\n", out);
  }
  fprintf(out, "Background letter frequencies (from %s):\n",
      options->bg_filename ? options->bg_filename : "unknown source");
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
  for (idx = 1, node = rbtree_first(all_motifs); node != NULL; node = rbtree_next(node)) {
    MOTIFS_T *motifs;
    char *id;
    id = (char*)rbtree_key(node);
    motifs = (MOTIFS_T*)rbtree_value(node);
    for (i = 0; i < motifs->count; i++, idx++) {
      MOTIF_T *motif;
      MATRIX_T *freqs, *scores;
      char *alt;
      // separate multiple motifs
      if (idx > 1) {
        fputs("\n\n", out);
      }
      motif = motifs->list[i];
      freqs = get_motif_freqs(motif);
      scores = get_motif_scores(motif);
      if (options->use_numbers) {
        // the easy option, we rename everthing using numbers
        str_setf(id_buf, "%d", idx);
      } else if (motifs->count > 1) {
        // use the orignal name but append a number to make it unique
        bool name_unique;
        j = i;
        do {
          str_setf(id_buf, "%s.%d", id, ++j);
          // check we don't have a motif with that name
          if (rbtree_find(all_motifs, str_internal(id_buf))) continue;
          // confirm we haven't output a motif using this name
          rbtree_lookup(all_names, str_internal(id_buf), true, &name_unique);
        } while(!name_unique);
      } else {
        // use the original name
        str_setf(id_buf, "%s", id);
      }
      fprintf(out, "MOTIF %s", str_internal(id_buf));
      alt = get_motif_id2(motif);
      if (alt != NULL && alt[0] != '\0') {
        fprintf(out, " %s", alt);
      }
      fputs("\n\n", out);
      // optionally output the log-odds
      if (options->print_logodds) {
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
      // output the URL
      if (!options->replace_url && has_motif_url(motif)) {
        fprintf(out, "\nURL %s\n", get_motif_url(motif));
      } else if (options->url_pattern) {
        fputs("\nURL ", out);
        for (i = 0; url_parts[i]; i++) {
          if (i > 0) fputs(str_internal(id_buf), out);
          fputs(url_parts[i], out);
        }
        fputs("\n", out);
      }
    }
  }
}

/**************************************************************************
 * Run the program
 **************************************************************************/
int main(int argc, char** argv) {
  OPTIONS_T options;
  process_arguments(argc, argv, &options);
  generate_minimal_meme(&options, stdout);
  return EXIT_SUCCESS;
}
