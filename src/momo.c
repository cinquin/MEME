/********************************************************************
 * MOMO Portal
 ********************************************************************/

#define DEFINE_GLOBALS

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "display.h"
#include "display_globals.h"
#include "dir.h"
#include "fasta-io.h"
#include "momo.h"
#include "momo-output.h"
#include "io.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "simple-getopt.h"
#include "utils.h"
#include "binomial.h"
#include "meme.h"
#include "read_seq_file.h"
#include "momo-input.h"
#include "momo-algorithm.h"

char* program_name = "momo";

/* local variables */
#define DATA_HASH_SIZE 100003

/***********************************************************************
 Initialize options processing
 ***********************************************************************/
static MOMO_OPTIONS_T init_options(ALGORITHM_T algorithm) {
  /* Make sure various options are set to NULL or defaults. */
  MOMO_OPTIONS_T options;
  
  options.algorithm = algorithm;
  options.phospho_filenames = arraylst_create();
  
  options.allow_clobber = FALSE;
  options.eliminate_repeats = TRUE;
  options.filter = FALSE;
  options.hash_fasta = TRUE;
  options.remove_unknowns = (algorithm == motifx) ? TRUE : FALSE;
  options.single_motif_per_mass = FALSE;
  
  options.command_line = NULL;
  options.filter_field = NULL;
  options.html_path = NULL;

  options.output_dirname = "momo_out";
  options.protein_database_filename = NULL;
  options.sequence_column = "sequence";
  options.text_path = NULL;

  options.filter_threshold = 0.0;
  options.score_threshold = (algorithm == motifx) ? 0.000001 : 0.0;
  
  options.fg_filetype = psm;
  options.bg_filetype = fasta;
  
  options.count_threshold = (algorithm == motifx) ? 20 : 0;
  options.eliminate_repeat_width = 7;
  options.max_iterations = 50;
  options.max_motifs = 100;
  options.min_occurrences = (algorithm == motifx) ? 0 : 5;
  options.no_decrease_stop_iteration = 10;
  options.hash_fasta_width = 6;
  options.width = 7;
  return options;
}

/***********************************************************************
 Initialize summmary
 ***********************************************************************/
static SUMMARY_T init_summary() {
  SUMMARY_T summary;
  summary.num_mod = 0;
  summary.num_modtype = 0;
  summary.num_mod_passing = 0;
  summary.num_modtype_passing = 0;
  
  summary.hash_fasta_table = NULL;
  summary.hash_fasta_table_keys = NULL;
  HASH_TABLE mod_table = hash_create(DATA_HASH_SIZE);	// hash table of mods
  ARRAYLST_T * mod_table_keys = arraylst_create(); // hash table keys
  summary.mod_table = mod_table;
  summary.mod_table_keys = mod_table_keys;
  
  // Initalize the protein alphabet within summary
  ALPH_T* alph = alph_protein();
  summary.alph = alph;
  
  // Get list of amino acids
  STR_T* alph_letters_string = str_create(MAX_ALPH_SIZE);
  const char* alph_letters = alph_string(summary.alph, alph_letters_string);
  summary.alph_letters = alph_letters;
  
  // cleanup
  char* stored_string = str_destroy(alph_letters_string, 1);
  return summary;
}

/***********************************************************************
 Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options(MOMO_OPTIONS_T* options) {
  myfree(options->command_line);
  myfree(options->html_path);
  myfree(options->text_path);
  if (options->filter) {
    myfree(options->filter_field);
  }
  arraylst_destroy(NULL, options->phospho_filenames);
}

/***********************************************************************
 Free memory allocated in summary
 ***********************************************************************/
static void cleanup_summary(SUMMARY_T* summary) {
  int i;

  // Each SUMMARY_T contains: alph, bg_freqs, mod_table/mod_table_keys, hash_fasta_table/hash_fasta_table_keys

  myfree(summary->alph_letters);
  
  // Clean up alph & bg_freqs
  alph_release(summary->alph);
  free_array(summary->bg_freqs);
  
  // Clean up mod_table/mod_table_keys
  for (i = 0; i < arraylst_size(summary->mod_table_keys); ++i) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, summary->mod_table_keys);
    MOD_INFO_T * modinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    cleanup_modinfo(modinfo);
  }
  arraylst_destroy(NULL, summary->mod_table_keys);
  hash_destroy(summary->mod_table);
  
  // Clean up hash_fasta_table/hash_fasta_table_keys if used.
  if (summary->hash_fasta_table != NULL) {
    for (i = 0; i < arraylst_size(summary->hash_fasta_table_keys); ++i) {
      HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, summary->hash_fasta_table_keys);
      ARRAYLST_T * currlist = (ARRAYLST_T *) hash_get_entry_value(hash_entry);
      arraylst_destroy(free, currlist);
    }
    arraylst_destroy(NULL, summary->hash_fasta_table_keys);
    hash_destroy(summary->hash_fasta_table);
  }
}

static char* get_usage_message(ALGORITHM_T algorithm) {
  if (algorithm == simple) {
    return "Usage: momo simple [options] <psm file>+\n"
    "\n"
    "   Options:\n"
    "     --bg-filetype fasta|prealigned (default fasta)\n"
    "     --eliminate-repeats [positive odd integer or 0 for no elimination] (default 7)\n"
    "     --fg-filetype psm|fasta|prealigned (default psm)\n"
    "     --filter [field],lt|le|eq|ge|gt,[threshold] (default no filter)\n"
    "     --hash-fasta [positive integer. 0 for linear search] (default 6)\n"
    "     --min-occurrences [positive integer] (default 5)\n"
    "     --o <output dir> (default=momo_out)\n"
    "     --oc <output dir> (default=momo_out)\n"
    "     --protein-database <protein fasta file> (default none)\n"
    "     --remove-unknowns T|F (default F)\n"
    "     --sequence-column [column name] (default \"sequence\")\n"
    "     --single-motif-per-mass T|F (default F)\n"
    "     --verbosity 1|2|3|4 (default 2)\n"
    "     --version (print the version and exit)\n"
    "     --width [positive odd integer] (default 7)\n"
    "\n";
  } else if (algorithm == motifx) {
    return "Usage: momo motifx [options] <psm file>+\n"
    "\n"
    "   Options:\n"
    "     --bg-filetype fasta|prealigned (default fasta)\n"
    "     --count-threshold [positive integer] (default 20)\n"
    "     --eliminate-repeats [positive odd integer or 0 for no elimination] (default 7)\n"
    "     --fg-filetype psm|fasta|prealigned (default psm)\n"
    "     --filter [field],lt|le|eq|ge|gt,[threshold] (default no filter)\n"
    "     --hash-fasta [positive integer. 0 for linear search] (default 6)\n"
    "     --o <output dir> (default=momo_out)\n"
    "     --oc <output dir> (default=momo_out)\n"
    "     --protein-database <protein fasta file> (default none)\n"
    "     --remove-unknowns T|F (default T)\n"
    "     --score-threshold [positive value] (default 10E-6)\n"
    "     --sequence-column [column name] (default \"sequence\")\n"
    "     --single-motif-per-mass T|F (default F)\n"
    "     --verbosity 1|2|3|4 (default 2)\n"
    "     --version (print the version and exit)\n"
    "     --width [positive odd integer] (default 7)\n"
    "\n";
  } else if (algorithm == modl) {
    return "Usage: momo modl [options] <psm file>+\n"
    "\n"
    "   Options:\n"
    "     --bg-filetype fasta|prealigned (default fasta)\n"
    "     --eliminate-repeats [positive odd integer or 0 for no elimination] (default 7)\n"
    "     --fg-filetype psm|fasta|prealigned (default psm)\n"
    "     --filter [field],lt|le|eq|ge|gt,[threshold] (default no filter)\n"
    "     --hash-fasta [positive integer. 0 for linear search] (default 6)\n"
    "     --max-iterations [positive integer] (default 50)\n"
    "     --max-motifs [positive integer] (default 100)\n"
    "     --no-decrease-stop-iteration [positive integer] (default 10)\n"
    "     --o <output dir> (default=momo_out)\n"
    "     --oc <output dir> (default=momo_out)\n"
    "     --protein-database <protein fasta file> (default none)\n"
    "     --remove-unknowns T|F (default T)\n"
    "     --sequence-column [column name] (default \"sequence\")\n"
    "     --single-motif-per-mass T|F (default F)\n"
    "     --verbosity 1|2|3|4 (default 2)\n"
    "     --version (print the version and exit)\n"
    "     --width [positive odd integer] (default 7)\n"
    "\n";
  }
  
  return "Usage: momo <algorithm> [options] <arguments>\n"
    "\n"
    "momo supports the following algorithms:\n"
    "simple\n"
    "motifx\n"
    "modl\n"
    "\n"
    "Options and arguments are specific to each command.\n"
    "Type \'momo <command>\' for details.\n";
}

/***********************************************************************
 Process command line options
 ***********************************************************************/
static MOMO_OPTIONS_T process_momo_command_line(
                                              int argc,
                                              char* argv[]
                                              ) {
  
  ALGORITHM_T algorithm = none;
  if (argc >= 2) {
    if (strcmp(argv[1], "simple") == 0) {
      algorithm = simple;
    } else if (strcmp(argv[1], "motifx") == 0) {
      algorithm = motifx;
    } else if (strcmp(argv[1], "modl") == 0) {
      algorithm = modl;
    }
  }
  
  MOMO_OPTIONS_T options = init_options(algorithm);
  
  // Define command line options.
  cmdoption const momo_simple_options[] = {
    {"bg-filetype", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"fg-filetype", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"min-occurrences", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
  };
  cmdoption const momo_motifx_options[] = {
    {"bg-filetype", REQUIRED_VALUE},
    {"count-threshold", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"fg-filetype", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"score-threshold", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
  };
  cmdoption const momo_modl_options[] = {
    {"bg-filetype", REQUIRED_VALUE},
    {"eliminate-repeats", REQUIRED_VALUE},
    {"fg-filetype", REQUIRED_VALUE},
    {"filter", REQUIRED_VALUE},
    {"hash-fasta", REQUIRED_VALUE},
    {"max-iterations", REQUIRED_VALUE},
    {"max-motifs", REQUIRED_VALUE},
    {"no-decrease-stop-iteration", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"protein-database", REQUIRED_VALUE},
    {"remove-unknowns", REQUIRED_VALUE},
    {"sequence-column", REQUIRED_VALUE},
    {"single-motif-per-mass", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"width", REQUIRED_VALUE},
  };
  
  // Define the usage message.
  options.usage = get_usage_message(algorithm);
  
  // Set verbosity to default
  verbosity = 2;
  
  // Get the options for the corresponding algorithm
  int option_index = 0;
  
  if (algorithm == simple) {
    const int num_options = sizeof(momo_simple_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_simple_options);
  } else if (algorithm == modl) {
    const int num_options = sizeof(momo_modl_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_modl_options);
  } else { // algorithm == none or algorithm == motifx
    const int num_options = sizeof(momo_motifx_options) / sizeof(cmdoption);
    simple_setopt(argc - 1, argv + 1, num_options, momo_motifx_options);
  }
  
  // Parse the command line.
  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;
    
    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    }
    else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }
    if (strcmp(option_name, "bg-filetype") == 0) {
      if (strcmp(option_value, "fasta") != 0 && strcmp(option_value, "prealigned") != 0) {
        die("Bg-filetype must be fasta or prealigned");
      }
      if (strcmp(option_value, "prealigned") == 0) {
        options.bg_filetype = prealigned;
        options.hash_fasta = FALSE;
      }
    }
    else if (strcmp(option_name, "count-threshold") == 0) {
      options.count_threshold = atoi(option_value);
    }
    else if (strcmp(option_name, "eliminate-repeats") == 0) {
      options.eliminate_repeat_width = atoi(option_value);
      options.eliminate_repeats = (options.eliminate_repeat_width > 0);
      if (options.eliminate_repeats && options.eliminate_repeat_width % 2 != 1) {
        die("Eliminate-repeats must be odd or 0.");
      }
    }
    else if (strcmp(option_name, "fg-filetype") == 0) {
      if (strcmp(option_value, "psm") != 0 && strcmp(option_value, "fasta") != 0 && strcmp(option_value, "prealigned") != 0) {
        die("Fg-filetype must be psm or fasta or prealigned");
      }
      if (strcmp(option_value, "fasta") == 0) {
        options.fg_filetype = fasta;
        options.hash_fasta = FALSE;
      } else if (strcmp(option_value, "prealigned") == 0) {
        options.fg_filetype = prealigned;
        options.hash_fasta = FALSE;
      }
    }
    else if (strcmp(option_name, "filter") == 0){
      options.filter = TRUE;
      char * value = strdup(option_value);
      char * pch = strtok(value, ",");
      options.filter_field = pch;
      pch = strtok(NULL, ",");
      if (strcmp(pch, "le") == 0) {
        options.filter_type = le;
      } else if (strcmp(pch, "lt") == 0) {
        options.filter_type = lt;
      } else if (strcmp(pch, "eq") == 0) {
        options.filter_type = eq;
      } else if (strcmp(pch, "gt") == 0) {
        options.filter_type = gt;
      } else if (strcmp(pch, "ge") == 0) {
        options.filter_type = ge;
      } else {
        die ("Error reading filter argument");
      }
      pch = strtok(NULL, ",");
      options.filter_threshold = atof(pch);
    }
    else if (strcmp(option_name, "hash-fasta") == 0){
      options.hash_fasta_width = atoi(option_value);
      options.hash_fasta = (options.hash_fasta_width > 0) ? TRUE : FALSE;
    }
    else if (strcmp(option_name, "max-iterations") == 0){
      options.max_iterations = atoi(option_value);
    }
    else if (strcmp(option_name, "max-motifs") == 0){
      options.max_motifs = atoi(option_value);
    }
    else if (strcmp(option_name, "no-decrease-stop-iteration") == 0){
      options.no_decrease_stop_iteration = atoi(option_value);
    }
    else if (strcmp(option_name, "min-occurrences") == 0){
      options.min_occurrences = atoi(option_value);
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      options.output_dirname = option_value;
      options.allow_clobber = FALSE;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      options.output_dirname = option_value;
      options.allow_clobber = TRUE;
    }
    else if (strcmp(option_name, "protein-database") == 0){
      options.protein_database_filename = option_value;
    }
    else if (strcmp(option_name, "score-threshold") == 0){
      options.score_threshold = atof(option_value);
    }
    else if (strcmp(option_name, "sequence-column") == 0){
      options.sequence_column = option_value;
    }
    else if (strcmp(option_name, "single-motif-per-mass") == 0){
      if (!(strcmp(option_value, "T") == 0 || strcmp(option_value, "F") == 0)) {
        die("Error reading single-motif-per-mass");
      }
      options.single_motif_per_mass = (strcmp(option_value, "T") == 0);
    }
    else if (strcmp(option_name, "remove-unknowns") == 0){
      if (!(strcmp(option_value, "T") == 0 || strcmp(option_value, "F") == 0)) {
        die("Error reading remove-unknowns");
      }
      options.remove_unknowns = (strcmp(option_value, "T") == 0);
    }
    else if (strcmp(option_name, "verbosity") == 0){
      verbosity = atoi(option_value);
    }
    else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    }
    else if (strcmp(option_name, "width") == 0){
      options.width = atoi(option_value);
      if (options.width % 2 != 1) {
        die("Width must be odd");
      }
    }
  }
  option_index++;
  
  // Must have algorithm
  if (algorithm == none) {
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  }
  
  // Cannot filter on non-psm input
  if (options.fg_filetype != psm && options.filter) {
    die("Cannot filter on non-tsv input");
  }
  
  // Must have peptide-spectrum-match (psm) file name
  if (argc < option_index + 1) {
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  }
  
  // Motif X must have protein database
  if ((algorithm == motifx || algorithm == modl) && ! options.protein_database_filename) {
    die("Must provide protein-database if using motif-x or modl algorithm");
  }
  
  // Record the command line
  options.command_line = get_command_line(argc, argv);
  
  // Record the input file names
  while (option_index != argc) {
    arraylst_add(argv[option_index], options.phospho_filenames);
    option_index++;
  }
  
  // Set up path values for needed stylesheets and output files.
  options.HTML_FILENAME = "momo.html";
  options.TEXT_FILENAME = "momo.txt";
  
  options.html_path = make_path_to_file(options.output_dirname, options.HTML_FILENAME);
  options.text_path = make_path_to_file(options.output_dirname, options.TEXT_FILENAME);
  
  return options;
  
}

/**
 * Sets the background frequencies. If protein database is provided, will count the number of 
 * occurrences of each amino acid and normalize. Otherwise will use NRDB frequencies.
 */
static void initialize_background_frequencies(MOMO_OPTIONS_T* options,
                                            SUMMARY_T* summary,
                                            SEQ_T** all_sequences,
                                            int num_sequences) {
  int i;
  int j;
  
  ARRAY_T * bg_freqs = NULL;
  const char* alph_letters = summary->alph_letters;
  
  if (options->protein_database_filename) {
    bg_freqs = allocate_array(strlen(alph_letters));
    init_array(0.0, bg_freqs);
    
    for (i = 0; i < num_sequences; ++i) {
      char* currseq = get_raw_sequence(all_sequences[i]);
      for (j = 0; j < strlen(currseq); ++j) {
        if (strchr(alph_letters, currseq[j])) {
          int idx = strchr(alph_letters, currseq[j]) - alph_letters;
          if (idx >= 0) {
            set_array_item_defcheck(idx, get_array_item_defcheck(idx, bg_freqs) + 1, bg_freqs);
          }
        }
      }
    }
    normalize(0.0, bg_freqs);
  } else {
    bg_freqs = get_nrdb_frequencies(summary->alph, bg_freqs);
  }
  summary->bg_freqs = bg_freqs;
}

/**
 * Analyzes the Peptide-Spectrum Match (PSM) files, create motifs, and returns a summary
 * on the results
 */
static SUMMARY_T get_summary(MOMO_OPTIONS_T* options) {
  int i;
  
  // Initialize the summary object
  SUMMARY_T summary = init_summary();
  
  // If a protein database is provided and in fasta format, we will read it into
  // an array of sequences.
  SEQ_T** all_sequences = NULL;
  int num_sequences = 0;
  read_protein_database_sequences(options, &summary, &all_sequences, &num_sequences);
  
  // See momo-input.c
  // If a protein database is provided and the user specified to use
  // an O(1) lookup table to speed up the process of finding a peptide
  // within the protein database, we initialize the O(1) lookup table.
  create_hash_fasta_preprocess_table(options, &summary, all_sequences, num_sequences);
  
  // Initialize background frequencies
  initialize_background_frequencies(options, &summary, all_sequences, num_sequences);
  
  // See momo-input.c
  // For each psm file, we will add the information to a mod table. Each mod is hashed
  // to its mod entry.
  add_phospho_files_to_table(options, &summary, all_sequences, num_sequences);
  
  // See momo-input.c
  // if using motifx, we want to set background sequences for each mod.
  add_background_sequences_to_table(options, &summary, all_sequences, num_sequences);
  
  // See momo-algorithm.c
  // Using the mod table, we will generate frequencies and create motifs for each mod.
  create_motifs(options, &summary, all_sequences, num_sequences);
  
  // Free up the memory used by sequences
  for (i = 0; i < num_sequences; ++i) {
    free_seq(all_sequences[i]);
  }
  if (all_sequences) {
    myfree(all_sequences);
  }
  
  return summary;
  
}

/*************************************************************************
 * Entry point for momo
 *************************************************************************/
int main(int argc, char *argv[]) {
  // Start timing
  clock_t start = clock(), diff;
  
  // Get command line arguments
  MOMO_OPTIONS_T options = process_momo_command_line(argc, argv);
  
  // Create motifs and obtain summary
  SUMMARY_T summary = get_summary(&options);
  
  // Print results
  print_momo_results(options, summary);

  // Clean up.
  cleanup_options(&options);
  cleanup_summary(&summary);
  
  // Print timing
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
//  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
  
  return 0;
}
