/********************************************************************
 * FILE: ame.c
 * AUTHOR: Robert McLeay
 * CREATE DATE: 19/08/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * ame is a yet unpublished algorithm that seeks to assist in
 * determining whether a given transcription factor regulates a set
 * of genes that have been found by an experimental method that
 * provides ranks, such as microarray or ChIp-chip.
 *
 *
 ********************************************************************/

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <libgen.h>
#include "alphabet.h"
#include "ame.h"
#include "array.h"
#include "cisml.h"
#include "config.h"
#include "fasta-io.h"
#include "fisher_exact.h"
#include "io.h"
#include "macros.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif_regexp.h"
#include "projrel.h"
#include "ramen_scan.h"
#include "ranksum_test.h"
#include "regress.h"
#include "simple-getopt.h"
#include "spearman-rank-correlation.h"
#include "string-list.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

static const double range = 100; // used for creating PSSMs
static const int gcbins = 1; // used for creating PSSMs

/*
 * Global Variables
 */

static rsc_motifs_t motifs;
static MATRIX_T* pssms;

static rsc_arg_t default_args;
static rsc_arg_t args;

static time_t  t0, t1; /* measuring time */
static clock_t c0, c1; /* measuring cpu_time */

static double** results;
static STRING_LIST_T* seq_ids;
static ARRAY_T* seq_fscores;
static double* factorial_array;

// maximum line length for output (set buffer size bigger than needed)
// no real global use of this, global for convenience: keep it that way
#define MAXLINE 1000
static char buffer [MAXLINE];

#define DEFAULTDIR "ame_out"

static char *default_output_dirname = DEFAULTDIR;  // default output directory name
static const char *text_filename = "ame.txt";
static const char *html_filename = "ame.html";
static const char *template_filename = "ame_template.html";

const char *pvalue_method_names[] = {"Wilcoxon rank-sum test", "Fisher's exact test", "multihg", 
  "long_multihg", "linear regression", "Spearman's correlation coefficient"};
const char *scoring_method_names[] = {"avg_odds", "max_odds", "sum_odds", 
  "total_hits"};
const char *bg_format_names[] = {"uniform", "motif", "file"};

static rsc_result_t * init_result (rsc_result_t *new_result, // NULL => create a new one and return it
    int db_idx,
    MOTIF_T *motif,           // motif to copy into result
    int split,                // where to split between foreground, background
    double pleft,             // p-values for left, right and both tails
    double pright,
    double pboth,
    double u);

static void rsc_print_results (char *format, ...); // send output to the text file
static void final_print_results (); // finalise outputs and close files

static char * my_strdup(const char *s1); // returns "" string if passed a NULL pointer

/*************************************************************************
 * Initializations of parameters to defaults
 *************************************************************************/

void rsc_set_defaults() {
  default_args.alph_filename = NULL;
  default_args.bg_filename = NULL;
  default_args.outputdir = default_output_dirname;
  default_args.bg_format = MOTIF_BG;
  default_args.motif_filenames = NULL;
  default_args.number_motif_files = 0;
  default_args.sequence_filename = NULL;
  default_args.control_filename = NULL;
  default_args.commandline = NULL;
  default_args.pseudocount = 0.25;
  default_args.scoring = TOTAL_HITS ; 
  default_args.verbose = NORMAL_VERBOSE;
  default_args.rs_method = QUICK_RS;
  default_args.positive_list = POS_FL;
  default_args.pvalue_method = FISHER_METHOD; 
  default_args.fisher_pwm_threshold = 1;
  default_args.fisher_fasta_threshold = 1e-3;
  default_args.pvalue_threshold = 2e-4;
  default_args.pvalue_report_threshold = 0.05;
  default_args.length_correction = FALSE;
  default_args.log_pwmscores = FALSE;
  default_args.log_fscores = FALSE;
  default_args.linreg_normalise = FALSE;
  default_args.linreg_dump_dir = NULL;
  default_args.linreg_switchxy = TRUE;
  default_args.clobber = TRUE;
  default_args.silent = FALSE;
  default_args.fix_partition = -1; //i.e. disabled.
  default_args.text_output = NULL;
  default_args.html_output = NULL;
  default_args.json_output = NULL;
  //Now copy the defaults into the real args
  memcpy(&args, &default_args, sizeof(rsc_arg_t));
}


/*************************************************************************
 * Entry point for AME
 *************************************************************************/
int main(int argc, char *argv[]) {

  int i;

  char *version_message = "AME (Analysis of Motif Enrichment): Compiled on " __DATE__ " at " __TIME__ "\n"
      "------------------------------\n"
      "Copyright (c) Robert McLeay <r.mcleay@imb.uq.edu.au> & Timothy Bailey <t.bailey@imb.uq.edu.au>, 2009.\n\n";

  /* Record the execution start and end times */
  if (verbosity >= HIGH_VERBOSE) {
    t0 = time(NULL);
    c0 = clock();
  }

  rsc_set_defaults(); //Set default cmd line args
  rsc_getopt(argc, argv); //Get command line args

  // skip this stuff if testing: don't really care if paths change or compiled on a different date
  if (!args.silent) {
    rsc_print_results ("%s", version_message);
    rsc_print_results ("Command line\n%s\n\n", args.commandline);
  }

  rsc_load_motifs(); //Load the motifs from the file
  rsc_scan_sequences(); //Load each sequence, scanning with each motif to determine score

  if(args.control_filename) {
    rsc_print_results("Not in partition maximization mode. "
	"Fixing partition at the number of primary sequences (%d).\n\n", args.fix_partition);
  } else if (args.fix_partition > 0) {
    rsc_print_results("Not in partition maximization mode. "
        "Fixing partition at %i.\n\n", args.fix_partition);
  } else {
    rsc_print_results("In partition maximization mode.\n\n");
  }

  rsc_print_results("Threshold p-value for reporting results: %g\n", args.pvalue_report_threshold);

  if (args.pvalue_method != RANKSUM_METHOD) {
    fisher_exact_init(get_num_strings(seq_ids)); // Generate the table of log factorials
    factorial_array = fisher_exact_get_log_factorials();
  }

  rsc_get_scores(); //Do the significance tests to associate each motif with the set.

  //print_cismlq(stdout, output);

  rsc_terminate(0); //Successfully end.

  return 0; //shuts up a warning.
}

// send line unaltered to the text output
void rsc_print_results (char* format, ...) {
  va_list arg_ptrs;
  va_start(arg_ptrs, format);
  vfprintf(args.text_output, format, arg_ptrs);
  va_end(arg_ptrs);
}

static void record_result(rsc_result_t* result, double corrected_pvalue) {
  int i, j, len, alen;
  MATRIX_T *freqs;
  ARRAY_T *row;
  jsonwr_start_object_value(args.json_output);
  jsonwr_lng_prop(args.json_output, "db", result->db_idx);
  jsonwr_str_prop(args.json_output, "id", get_motif_id(result->motif));
  if (get_motif_id2(result->motif)[0] != '\0')
    jsonwr_str_prop(args.json_output, "alt", get_motif_id2(result->motif));
  jsonwr_lng_prop(args.json_output, "len", get_motif_length(result->motif));
  jsonwr_dbl_prop(args.json_output, "motif_evalue", get_motif_evalue(result->motif));
  jsonwr_dbl_prop(args.json_output, "motif_nsites", get_motif_nsites(result->motif));
  if (has_motif_url(result->motif)) 
    jsonwr_str_prop(args.json_output, "url", get_motif_url(result->motif));
  alen = alph_size_core(get_motif_alph(result->motif));
  len = get_motif_length(result->motif);
  freqs = get_motif_freqs(result->motif);
  jsonwr_property(args.json_output, "pwm");
  jsonwr_start_array_value(args.json_output);
  for (i = 0; i < len; i++) {
    row = get_matrix_row(i, freqs);
    jsonwr_start_array_value(args.json_output);
    for (j = 0; j < alen; j++) {
      jsonwr_dbl_value(args.json_output, get_array_item(j, row));
    }
    jsonwr_end_array_value(args.json_output);
  }
  jsonwr_end_array_value(args.json_output);
    
  if (args.fix_partition <= 0) 
    jsonwr_lng_prop(args.json_output, "split", result->split);
  if (args.pvalue_method == LINREG_METHOD) {
    jsonwr_dbl_prop(args.json_output, "mean_square_error", result->pboth);
    jsonwr_dbl_prop(args.json_output, "m", result->pleft);
    jsonwr_dbl_prop(args.json_output, "b", result->pright);
  } else if (args.pvalue_method == SPEARMAN_METHOD) {
    jsonwr_dbl_prop(args.json_output, "spearmans_rho", result->pboth);
  } else { // RANKSUM_METHOD, FISHER_METHOD, MULTIHG_METHOD, LONG_MULTIHG_METHOD
    jsonwr_dbl_prop(args.json_output, "pvalue", result->pright);
    jsonwr_dbl_prop(args.json_output, "corrected_pvalue", corrected_pvalue);
  }
  jsonwr_end_object_value(args.json_output);
}

// close the output files
static void final_print_results () {
  fclose(args.text_output);
  if (args.html_output) {
    if (htmlwr_output(args.html_output) != NULL) {
      die("Expected only one replacement variable in template.\n");
    }
    htmlwr_destroy(args.html_output);
    args.html_output = NULL;
  }
}


void rsc_getopt(int argc, char *argv[]) {
  const int num_options = 26; // change this if the number of options changes
  cmdoption const motif_scan_options[] = {
    {"xalph", REQUIRED_VALUE},
    {"control", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"bgformat", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"pseudocount", REQUIRED_VALUE},
    {"fasta-threshold", REQUIRED_VALUE},
    {"pwm-threshold", REQUIRED_VALUE},
    {"pvalue-threshold", REQUIRED_VALUE},
    {"pvalue-report-threshold", REQUIRED_VALUE},
    {"scoring", REQUIRED_VALUE},
    {"method", REQUIRED_VALUE},
    {"rsmethod", REQUIRED_VALUE},
    {"poslist", REQUIRED_VALUE},
    {"verbose", REQUIRED_VALUE},
    {"silent", NO_VALUE},
    {"length-correction", NO_VALUE},
    {"linreg-dumpdir", REQUIRED_VALUE},
    {"linreg-switchxy", NO_VALUE},
    {"log-fscores", NO_VALUE},
    {"log-pwmscores", NO_VALUE},
    {"normalise-affinity", NO_VALUE},
    {"fix-partition", REQUIRED_VALUE},
    {"help", NO_VALUE},
    {"version", NO_VALUE}
  };

  int option_index = 0;
  char* option_name = NULL;
  char* option_value = NULL;
  const char * message = NULL;
  BOOLEAN_T bad_options = FALSE;
  int i;

  if (simple_setopt(argc, argv, num_options, motif_scan_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }

  /*
   * Now parse the command line options
   */
  //simple_getopt will return 0 when there are no more options to parse
  while(simple_getopt(&option_name, &option_value, &option_index) > 0) {
    //fprintf(stderr, "Option name: %s || %s\n", option_name, option_value);
    if (args.verbose >= HIGHER_VERBOSE) {
      fprintf(stderr, "ame Option: %s - %s\n",option_name, option_value);
    }
    if (strcmp(option_name, "xalph") == 0) {
      args.alph_filename = option_value;
    } else if (strcmp(option_name, "control") == 0) {
      args.control_filename = option_value;
    } else if (strcmp(option_name, "bgfile") == 0) {
      args.bg_filename = option_value;
    } else if (strcmp(option_name, "o") == 0) {
      args.outputdir = option_value;
    } else if (strcmp(option_name, "oc") == 0) {
      args.outputdir = option_value;
      args.clobber = TRUE;
    } else if (strcmp(option_name, "bgformat") == 0) {
      if (atoi(option_value)==MOTIF_BG) {
        args.bg_format = MOTIF_BG;
      } else if (atoi(option_value)==FILE_BG) {
        args.bg_format = FILE_BG;
      } else if (atoi(option_value)==UNIFORM_BG) {
        args.bg_format = UNIFORM_BG;
      } else {
        rsc_usage();
        rsc_terminate(1);
      }

    } else if (strcmp(option_name, "rsmethod") == 0) {
      if (strcmp(option_value,"quick")==0) {
        args.rs_method = QUICK_RS;
        fprintf(stderr, "Using quick ranksum calculation method.\n");
      } else if (strcmp(option_value,"better") == 0) {
        args.rs_method = BETTER_RS;
        fprintf(stderr, "More accurate ranksum calculation method not yet supported.\n");
        exit(EXIT_FAILURE);
      } else {
        rsc_usage();
        rsc_terminate(1);
      }
    } else if (strcmp(option_name, "pseudocount") == 0) {
      args.pseudocount = atof(option_value);
    } else if (strcmp(option_name, "pwm-threshold") == 0) {
      args.fisher_pwm_threshold = atof(option_value);
    } else if (strcmp(option_name, "fasta-threshold") == 0) {
      args.fisher_fasta_threshold = atof(option_value);
    } else if (strcmp(option_name, "pvalue-threshold") == 0) {
      args.pvalue_threshold = atof(option_value);
    } else if (strcmp(option_name, "pvalue-report-threshold") == 0) {
      args.pvalue_report_threshold = atof(option_value);
    } else if (strcmp(option_name, "scoring") == 0) {
      if (strcmp(option_value,"avg")==0) {
        args.scoring = AVG_ODDS;
        fprintf(stderr, "Using average odds scoring.\n");
      } else if (strcmp(option_value,"sum") == 0) {
        args.scoring = SUM_ODDS;
        fprintf(stderr, "Using sum of odds scoring.\n");
      } else if (strcmp(option_value,"max") == 0) {
        args.scoring = MAX_ODDS;
        fprintf(stderr, "Using maximum odds scoring.\n");
      } else if (strcmp(option_value,"totalhits") == 0) {
        args.scoring = TOTAL_HITS;
        fprintf(stderr, "Using total hits scoring.\n");
      } else {
        rsc_usage();
        rsc_terminate(1);
      }
    } else if (strcmp(option_name, "poslist") == 0) {
      if (strcmp(option_value,"fasta")==0) {
        args.positive_list = POS_FL;
        fprintf(stderr, "Using FASTA scores as positive.\n");
      } else if (strcmp(option_value,"pwm") == 0) {
        args.positive_list = POS_PWM;
        fprintf(stderr, "Using PWM scores as positive.\n");
      } else {
        rsc_usage();
        rsc_terminate(1);
      }
    } else if (strcmp(option_name, "method") == 0) {
      if (strcmp(option_value,"ranksum")==0) {
        args.pvalue_method = RANKSUM_METHOD;
      } else if (strcmp(option_value,"fisher") == 0) {
        args.pvalue_method = FISHER_METHOD;
      } else if (strcmp(option_value,"mhg") == 0) {
        args.pvalue_method = MULTIHG_METHOD;
      } else if (strcmp(option_value,"4dmhg") == 0) {
        args.pvalue_method = LONG_MULTIHG_METHOD;
      } else if (strcmp(option_value,"linreg") == 0) {
        args.pvalue_method = LINREG_METHOD;
      } else if (strcmp(option_value,"spearman") == 0) {
        args.pvalue_method = SPEARMAN_METHOD;
      } else {
        rsc_usage();
        rsc_terminate(1);
      }
    } else if (strcmp(option_name, "verbose") == 0) {
      args.verbose = atoi(option_value);
      verbosity = args.verbose;
      if (args.verbose <= INVALID_VERBOSE || args.verbose > DUMP_VERBOSE) {
        rsc_usage();
        rsc_terminate(1);
      }
    } else if (strcmp(option_name, "silent") == 0) {
      args.silent = TRUE; // omit from usage message: for generating test output
    } else if (strcmp(option_name, "length-correction") == 0) {
      args.length_correction = TRUE;
    } else if (strcmp(option_name, "log-fscores") == 0) {
      args.log_fscores = TRUE;
    } else if (strcmp(option_name, "log-pwmscores") == 0) {
      args.log_pwmscores = TRUE;
    } else if (strcmp(option_name, "normalise-affinity") == 0) {
      args.linreg_normalise = TRUE;
      fprintf(stderr, "Normalising per-motif scores.\n");
    } else if (strcmp(option_name, "linreg-dumpdir") == 0) {
      args.linreg_dump_dir = option_value;    
    } else if (strcmp(option_name, "fix-partition") == 0) {
      args.fix_partition = atoi(option_value);
    } else if (strcmp(option_name, "linreg-switchxy") == 0) {
      args.linreg_switchxy = FALSE;
    } else if (strcmp(option_name, "help") == 0) {
      fprintf(stderr, "%s",rsc_get_usage()); //was not to stderr
      rsc_terminate(0);
    } else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    } else {
      fprintf(stderr, "Error: %s is not a recognised switch.\n", option_name);
      rsc_usage();
      rsc_terminate(1);
    }

    option_index++;
  }

  // Set --bgformat if --bgfile given.
  if (args.bg_filename) {
    args.bg_format = FILE_BG;
  }

  if (args.pvalue_method == LINREG_METHOD || args.pvalue_method == SPEARMAN_METHOD) {
    if (args.linreg_switchxy) {
      fprintf(stderr, "In LR/Spearman mode, x=PWM, y=FASTA score\n");
    } else {
      //Standard mode
      fprintf(stderr, "In LR/Spearman mode, x=FASTA score, y=PWM\n");
    }
  }

  // Must have sequence and motif files
  if (argc < option_index + 2) {
    fprintf(stderr, "Error: Must specify a sequence file and motif file(s).\n");
    rsc_usage();
    rsc_terminate(1);
  }
  args.sequence_filename = argv[option_index];
  args.motif_filenames = &argv[option_index+1]; // FIXME: must now iterate until argc getting all motif DB filenames
  args.number_motif_files = argc-option_index-1;
  if (args.bg_format == MOTIF_BG) { //bgfile is the same as the motif file.
    //TODO: fix ?? FIXME --  motif file name not used so no problem changing to motif DB of multiple files
    args.bg_filename = NULL;
  }

  /* Now validate the options.
   *
   * Illegal combinations are:
   *   - (control or spearman) and fix-partition
   *   - control and (spearman or linreg)
   *   - sequence bg and no sequence
   *   - no motif
   *   - no sequences
   *   - each file exists.
   */
 
  if (args.control_filename && args.fix_partition > 0) {
    die("Error: You may not specify '--fix-partition' and '--control'.\n");
  }
  if (args.fix_partition > 0 && args.pvalue_method == SPEARMAN_METHOD) {
    die("Error: You may not specify '--fix-partition' and '--method spearman'.\n");
  }

  if (args.control_filename && args.pvalue_method == SPEARMAN_METHOD) {
    die("Error: You may not specify '--control' and '--method spearman'.\n");
  }
  if (args.control_filename && args.pvalue_method == LINREG_METHOD) {
    die("Error: You may not specify '--control' and '--method linreg'.\n");
  }

  if (args.motif_filenames == NULL) {
    fprintf(stderr, "Error: Motif file not specified.\n");
    bad_options = TRUE;
  } else {
    int i;
    for (i = 0; i < args.number_motif_files; i++) {
      if (!file_exists(args.motif_filenames[i])) {
        fprintf(stderr, "Error: Specified motif '%s' file does not exist.\n", args.motif_filenames[i]);
        bad_options = TRUE;
      }
    }
  }
  if (args.sequence_filename == NULL) {
    fprintf(stderr, "Error: Sequence file not specified.\n");
    bad_options = TRUE;
  } else if (!file_exists(args.sequence_filename)) {
    fprintf(stderr, "Error: Specified sequence file does not exist.\n");
    bad_options = TRUE;
  }

  if (args.positive_list == POS_PWM && args.rs_method == BETTER_RS) {
    fprintf(stderr, "Error: You may not specify '--rsmethod better' with '--poslist pwm'.\n");
    bad_options = TRUE;
  }

  if (args.pvalue_method == MULTIHG_METHOD && args.scoring != TOTAL_HITS) {
    fprintf(stderr, "Error: You must use '--scoring totalhits' with '--method mhg'.\n");
    bad_options = TRUE;
  }

  if (args.length_correction && args.scoring != TOTAL_HITS) {
    fprintf(stderr, "Error: You must use '--scoring totalhits' with '--length-correction'.\n");
    bad_options = TRUE;
  } else if (args.length_correction) {
    fprintf(stderr, "Enabling length correction mode for total hits scoring.\n");
  }

  if (bad_options) {
    //rsc_usage();
    fprintf(stderr, "Type 'ame --help' to see the complete command usage message.\n");
    rsc_terminate(1);
  }
  // make enough space for all the command line args, with one space between each
  int line_length = 0;
  for (i = 0; i < argc; i++)
    line_length += strlen(i == 0 ? basename(argv[0]) : argv[i]);
  // add on argc to allow one char per word for separating space + terminal '\0'
  args.commandline = (char*) malloc(sizeof(char)*(line_length+argc));
  int nextpos = 0;
  for (i = 0; i < argc; i++) {
    // been here before? put in a space before adding the next word
    if (nextpos) {
      args.commandline[nextpos] = ' ';
      nextpos ++;
    }
    char * nextword = i == 0 ? basename(argv[0]) : argv[i];
    strcpy(&args.commandline[nextpos], nextword);
    nextpos += strlen (nextword);
  }
  if (create_output_directory(args.outputdir, args.clobber, 
        args.verbose > QUIET_VERBOSE)) { // only warn in higher verbose modes
    fprintf(stderr, "failed to create output directory `%s' or already exists\n", args.outputdir);
    exit(1);
  }
  char *path = make_path_to_file(args.outputdir, text_filename);
  args.text_output = fopen(path, "w"); //FIXME check for errors: MEME doesn't either and we at least know we have a good directory
  myfree(path);
  // setup the html output writer 
  args.html_output = htmlwr_create(get_meme_etc_dir(), template_filename, false);
  htmlwr_set_dest_name(args.html_output, args.outputdir, html_filename);
  htmlwr_replace(args.html_output, "ame_data.js", "data");
  args.json_output = htmlwr_output(args.html_output);
  // write out some information
  jsonwr_str_prop(args.json_output, "version", VERSION);
  jsonwr_str_prop(args.json_output, "revision", REVISION);
  jsonwr_str_prop(args.json_output, "release", ARCHIVE_DATE);
  jsonwr_str_array_prop(args.json_output, "cmd", argv, argc);
  // options
  jsonwr_property(args.json_output, "options");
  jsonwr_start_object_value(args.json_output);
  jsonwr_str_prop(args.json_output, "scoring", scoring_method_names[args.scoring]);
  jsonwr_str_prop(args.json_output, "rs_method", (args.rs_method ? "better" : "quick"));
  jsonwr_str_prop(args.json_output, "positive_list", (args.positive_list ? "MOTIF" : "FASTA"));
  jsonwr_str_prop(args.json_output, "pvalue_method", pvalue_method_names[args.pvalue_method]);
  jsonwr_dbl_prop(args.json_output, "pseudocount", args.pseudocount);
  jsonwr_dbl_prop(args.json_output, "fisher_pwm_threshold", args.fisher_pwm_threshold);
  jsonwr_dbl_prop(args.json_output, "fisher_fasta_threshold", args.fisher_fasta_threshold);
  jsonwr_dbl_prop(args.json_output, "pvalue_threshold", args.pvalue_threshold);
  jsonwr_dbl_prop(args.json_output, "pvalue_report_threshold", args.pvalue_report_threshold);
  jsonwr_bool_prop(args.json_output, "length_correction", args.length_correction);
  jsonwr_bool_prop(args.json_output, "log_fscores", args.log_fscores);
  jsonwr_bool_prop(args.json_output, "log_pwmscores", args.log_pwmscores);
  jsonwr_bool_prop(args.json_output, "linreg_normalise", args.linreg_normalise);
  jsonwr_bool_prop(args.json_output, "linreg_switchxy", args.linreg_switchxy);
  jsonwr_bool_prop(args.json_output, "fix_partition", 
    (args.control_filename || args.pvalue_method == SPEARMAN_METHOD || (args.fix_partition > 0)));
  if (args.control_filename) {
    jsonwr_str_prop(args.json_output, "partition", "number of primary sequences");
  } else if (args.pvalue_method == SPEARMAN_METHOD) {
    jsonwr_str_prop(args.json_output, "partition", "total number of sequences");
  } else if (args.fix_partition > 0) {
    jsonwr_lng_prop(args.json_output, "partition", args.fix_partition);
  }
  jsonwr_end_object_value(args.json_output);
}

const char* rsc_get_usage() {
  // Define the usage message.
  //TODO: sprintf in default values
  return
    "Usage: ame [options] <sequence file> <motif file>+\n"
    "	  <sequence file> sorted sequences in FASTA format\n"
    "	  <motif> 	  motifs in MEME format\n"
    "\n"
    "     --o  <output dir> (default=ame_out)\n"
    "     --oc <output dir> (default=ame_out)\n"
    "     --control <control file>	control sequences in FASTA format\n"
    "     --method  [fisher|mhg|4dmhg|ranksum|linreg|spearman] (default: fisher)\n"
    "     --scoring [avg|max|sum|totalhits] (default: totalhits)\n"
    "     --xalph <alph file> motifs will be converted to this custom alphabet\n"
    "     --bgformat [0|1|2] 0:uniform, 1: motif file freqs, 2: bgfile freqs(default: 1)\n"
    "     --bgfile <background>, overrides -bgformat (default: motif file freqs)\n"
    "     --length-correction (default: don't correct for sequence length\n"
    "     --pvalue-threshold <pvt> (default: 0.0002)\n"
    "     --pvalue-report-threshold <pv> (default: 0.05)\n"
    "     --pwm-threshold <pwmt> (default: 1.0) (Requires --poslist pwm)\n"
    "     --fasta-threshold  <ft> (default: 0.001) (Requires --poslist fasta)\n"
    "     --fix-partition <int> (default: unconstrained partition maximization)\n"
    "     --poslist       [fasta|pwm] (default: fasta)\n"
    "     --rsmethod [better|quick] (default: quick) (Requires ranksum)\n"
    "     --normalise-affinity (default: don't normalize motif scores (linreg or spearman) \n"
    "     --linreg-switchxy (default: don't switch (X=FASTA scores, Y=motif scores)\n"
    "     --log-fscores (default: regress on FASTA scores (linreg) or ranks (spearman))\n"
    "     --log-pwmscores (default: regress on the motif scores)\n"
    "     --pseudocount <float, default = 0.25> Pseudocount for motif affinity scan\n"
    "     --verbose     <1...5> Integer describing verbosity.\n"
    "     --help        Print this message and exit.\n"
    "     --version     Print the version and exit.\n"
    "\n";
}

void rsc_usage() {
  fprintf(stderr, "%s", rsc_get_usage());
}

/*
 * TODO: This code taken from ama.c - let's rewrite it properly.
 */
void rsc_load_motifs() {
  ARRAYLST_T* read_motifs;
  MREAD_T *mread;
  int num_motifs_before_rc = 0;
  int num_motifs_after_rc = 0;
  int i, j, before_count, db_size;
  int *db_sizes;
  char *bg_src;
  ALPH_T *alph;
  BOOLEAN_T use_rc;

  db_sizes = malloc(sizeof(int) * args.number_motif_files);
  before_count = 0;

  if (args.alph_filename != NULL) {
    alph = alph_load(args.alph_filename, true);
    if (alph == NULL) exit(EXIT_FAILURE);
    switch (args.bg_format) {
      case MOTIF_BG:
        fprintf(stderr, "Warning, can't use motif background when converting motif alphabet.\n");
        args.bg_format = UNIFORM_BG; // can't use motif bg with alphabet
      case UNIFORM_BG:
        motifs.bg_freqs = get_uniform_frequencies(alph, NULL);
        break;
      case FILE_BG:
        motifs.bg_freqs = get_file_frequencies(alph, args.bg_filename);
        break;
      default:
        die("Illegal bg_format value.\n");
    }
    bg_src = NULL; // make compiler happy
  } else {
    alph = NULL;
    switch (args.bg_format) {
      case UNIFORM_BG:
        bg_src = "--uniform--";
        break;
      case MOTIF_BG:
        bg_src = "--motif--";
        break;
      case FILE_BG:
        bg_src = args.bg_filename;
        break;
      default:
        die("Illegal bg_format value.\n");
        bg_src = NULL; // make compiler happy
    }
  }

  // write motif dbs to html data section
  jsonwr_property(args.json_output, "motif_dbs");
  jsonwr_start_array_value(args.json_output);
  // clear the motifs structure
  memset(&motifs, 0, sizeof(rsc_motifs_t));
  read_motifs = arraylst_create();
  for (i = 0; i < args.number_motif_files; i++) {
    mread = mread_create(args.motif_filenames[i], OPEN_MFILE);
    if (alph != NULL) {
      mread_set_conversion(mread, alph, motifs.bg_freqs);
    } else {
      if (i == 0) mread_set_bg_source(mread, bg_src);
      else mread_set_background(mread, motifs.bg_freqs);
    }
    mread_set_pseudocount(mread, args.pseudocount);

    mread_load(mread, read_motifs);
    if (arraylst_size(read_motifs) == 0) die("No motifs found.");
    if (motifs.bg_freqs == NULL) motifs.bg_freqs = mread_get_background(mread);

    mread_destroy(mread);
    // keep track of the motif counts
    db_size = arraylst_size(read_motifs) - before_count;
    db_sizes[i] = db_size * 2; // double for RC motifs
    before_count = arraylst_size(read_motifs);
    // write motif db to html data
    jsonwr_start_object_value(args.json_output);
    jsonwr_str_prop(args.json_output, "source", args.motif_filenames[i]);
    jsonwr_lng_prop(args.json_output, "count", db_size);
    jsonwr_end_object_value(args.json_output);
  }
  // finish writing motif dbs
  jsonwr_end_array_value(args.json_output);
  // write background to json data
  if (alph == NULL && arraylst_size(read_motifs) > 0) {
    alph = get_motif_alph((MOTIF_T*)arraylst_get(0, read_motifs));
  }
  jsonwr_property(args.json_output, "alphabet");
  alph_print_json(alph, args.json_output);
  jsonwr_property(args.json_output, "background");
  jsonwr_start_object_value(args.json_output);
  jsonwr_str_prop(args.json_output, "source", bg_format_names[args.bg_format]);
  if (args.bg_format == FILE_BG)
    jsonwr_str_prop(args.json_output, "file", args.bg_filename);
  jsonwr_property(args.json_output, "frequencies");
  jsonwr_start_array_value(args.json_output);
  for (i = 0; i < alph_size_core(alph); i++) {
    jsonwr_dbl_value(args.json_output, get_array_item(i, motifs.bg_freqs));
  }
  jsonwr_end_array_value(args.json_output);
  jsonwr_end_object_value(args.json_output);

  // reverse complement the originals adding to the original read in list
  num_motifs_before_rc = arraylst_size(read_motifs);

  // Add reverse complement motifs, interleaving with originals.
  // Check if alphabet uses RCs:
  use_rc = alph_has_complement(alph);
  if (use_rc) {
    add_reverse_complements(read_motifs);        
    motifs.revcomp = TRUE;
  }

  //Allocate array for the motifs
  motif_list_to_array(read_motifs, &(motifs.motifs), &(num_motifs_after_rc));
  motifs.num_b4_rc = num_motifs_before_rc;

  //free the list of motifs
  free_motifs(read_motifs);
  
  // check reverse complements.
  if (use_rc) {
    assert(num_motifs_after_rc / 2 == num_motifs_before_rc); 
  }

  // Allocate a pssm, odds matrix and db index for each
  // including reverse complements.
  motifs.pssms = malloc(sizeof(PSSM_T) * num_motifs_after_rc);
  motifs.odds = malloc(sizeof(MATRIX_T*) * num_motifs_after_rc); 
  motifs.db_idx = malloc(sizeof(int) * num_motifs_after_rc);
  memset(motifs.pssms, 0, sizeof(PSSM_T) * num_motifs_after_rc);
  memset(motifs.odds, 0, sizeof(MATRIX_T*) * num_motifs_after_rc);
  memset(motifs.db_idx, 0, sizeof(int) * num_motifs_after_rc);

  //Allocate our pv_lookup space if required.
  if (TOTAL_HITS == args.scoring) {
    motifs.pv_lookup = malloc(sizeof(ARRAY_T*) * num_motifs_before_rc); 
  }

  // We need to create odds matrices if we're doing AMA-type scoring
  // Or, we create PSSMs and pv_lookup tables for the total-hits style instead.
  for (i = 0, j = 0, before_count = 0; i < num_motifs_after_rc; i++) {
    // setup db index tracking
    if ((i - before_count) >= db_sizes[j]) {
      before_count += db_sizes[j];
      j++;
    }
    motifs.db_idx[i] = j;
    // create odds or log-odds matrix
    if (TOTAL_HITS == args.scoring) {
      // Build and scale PSSM for each motif
      PSSM_T *newpssm = build_motif_pssm(motif_at(motifs.motifs, i),
          motifs.bg_freqs, motifs.bg_freqs, NULL, 1.0, range, gcbins, FALSE);
      motifs.pssms[i] = *newpssm;
      // If using reverse complements check that RC matrix has same number or rows
      if (use_rc && (i % 2 != 0)) {
        assert(get_num_rows(motifs.pssms[i].matrix) == get_num_rows(motifs.pssms[i-1].matrix)); 
      }
    } else {
      //convert_to_odds_matrix(motif_at(motifs.motifs, i), motifs.bg_freqs);
      motifs.odds[i] = create_odds_matrix(motif_at(motifs.motifs, i), motifs.bg_freqs);
      // If using reverse complements check that RC matrix has same number or rows
      if (use_rc && (i % 2 != 0)) {
        assert(get_num_rows(motifs.odds[i]) == get_num_rows(motifs.odds[i-1])); 
      }
    }
  }

  // cleanup
  free(db_sizes);
}

void rsc_scan_sequences() {
  FILE* seq_file = NULL;
  FILE* control_file = NULL;
  SEQ_T* sequence = NULL;
  int i;
  int j;
  SEQ_T** seq_list;
  int num_seqs;
  int seq_len;
  //For the bdb_bg mode:
  ARRAY_T* seq_bg_freqs;
  MATRIX_T* pssm;
  MATRIX_T* rev_pssm;
  double atcontent;
  double roundatcontent;
  char atcontentstr[6];
  char keydata[16];
  double scaled_pvalue_threshold = args.pvalue_threshold;
  double avg_seq_length = 0;
  ALPH_T *alph;
  BOOLEAN_T use_rc = motifs.revcomp;

  // get the alphabet from the global motifs
  alph = NULL;
  if (motifs.num_b4_rc > 0) {
    alph = get_motif_alph(motif_at(motifs.motifs, 0));
  }
  assert(alph != NULL);

  //Open the (primary) sequence file.
  if (open_file(args.sequence_filename, "r", FALSE, "FASTA", "sequences", &seq_file) == 0) {
    fprintf(stderr, "Couldn't open the (primary) sequences file %s.\n", args.sequence_filename);
    rsc_terminate(1);
  }

  //Open the control sequence file.
  if (args.control_filename &&
      open_file(args.control_filename, "r", FALSE, "FASTA", "sequences", &control_file) == 0) {
    fprintf(stderr, "Couldn't open the control sequences file %s.\n", args.control_filename);
    rsc_terminate(1);
  }

  //Read in the primary sequences.
  num_seqs = 0;
  read_many_fastas(alph, seq_file, MAX_SEQ_LENGTH, &num_seqs, &seq_list);

  // Get width of widest motif.
  int max_width = 0;
  for (i=0;i<motifs.num_b4_rc;i++) {
    int imotif = use_rc ? 2*i : i;
    int w = get_motif_length(motif_at(motifs.motifs, imotif));
    if (w > max_width) max_width = w;
  }

  // Remove primary sequences that are shorter than the widest motif.
  for (i=j=0; i<num_seqs; i++) {
    int length = get_seq_length(seq_list[i]);
    if (length >= max_width) {
      seq_list[j++] = seq_list[i];
    } else {
      char *name = get_seq_name(seq_list[i]);
      if (args.verbose >= NORMAL_VERBOSE) {
	fprintf(stderr, "Removing primary sequence %s since it is shorter (%d) than widest motif (%d)\n", 
	  name, length, max_width);
      }
    }
  }
  if (j == 0) {
    fprintf(stderr, "Error: All control sequences are shorter than the longest motif (%d).\n", 
      max_width);
    rsc_terminate(1);
  }
  if (j != num_seqs) {
    if (args.verbose >= NORMAL_VERBOSE) 
      fprintf(stderr, "Warning: Removed %d primary sequence(s) that were shorter than the longest motif (%d)\n", 
        num_seqs-j, max_width);
  }
  num_seqs = j;			// too short sequences removed

  // write to html data section
  jsonwr_property(args.json_output, "sequence_db");
  jsonwr_start_object_value(args.json_output);
  jsonwr_str_prop(args.json_output, "source", args.sequence_filename);
  jsonwr_lng_prop(args.json_output, "count", num_seqs);
  jsonwr_end_object_value(args.json_output);

  //Read in the control sequences, appending them to the list.
  if (args.control_filename) {
    args.fix_partition = num_seqs;
    int old_num_seqs = num_seqs;
    read_many_fastas(alph, control_file, MAX_SEQ_LENGTH, &num_seqs, &seq_list);
    // Remove control sequences that are shorter than the widest motif.
    for (i=j=old_num_seqs; i<num_seqs; i++) {
      int length = get_seq_length(seq_list[i]);
      if (length >= max_width) {
	seq_list[j++] = seq_list[i];
      } else {
	char *name = get_seq_name(seq_list[i]);
	if (args.verbose >= NORMAL_VERBOSE) {
	  fprintf(stderr, "Removing control sequence %s since it is shorter (%d) than widest motif (%d)\n", 
	    name, length, max_width);
        }
      }
    }
    if (j == old_num_seqs) {
      fprintf(stderr, "Error: All control sequences are shorter than the longest motif (%d).\n", 
	max_width);
      rsc_terminate(1);
    }
    if (j != num_seqs) {
      if (args.verbose >= NORMAL_VERBOSE) 
        fprintf(stderr, "Warning: Removed %d control sequence(s) that were shorter than the longest motif (%d)\n", 
  	  num_seqs-j, max_width);
    }
    num_seqs = j;			// too short sequences removed

    // write to html data section
    jsonwr_property(args.json_output, "control_db");
    jsonwr_start_object_value(args.json_output);
    jsonwr_str_prop(args.json_output, "source", args.control_filename);
    jsonwr_lng_prop(args.json_output, "count", num_seqs-old_num_seqs);
    jsonwr_end_object_value(args.json_output);
  }

  //Set fix-partition if in Spearman mode
  if (args.pvalue_method == SPEARMAN_METHOD) {
    args.fix_partition = num_seqs;
  }

  //Check that fix_partition value is legal.
  if (args.fix_partition > num_seqs) {
    die("'--fix-partition' value cannot be larger than the total number of input sequences.\n");
  }

  seq_ids = new_string_list();
  seq_fscores = allocate_array(num_seqs);

  //FIXME: This is extremely memory inefficient--saving the score of every
  // sequence using every motif in one big array.
  // Score all the sequences using all the arrays and store in Global results[motif][seq].

  //Allocate the required space for results
  results = malloc(sizeof(double*) * motifs.num_b4_rc);
  for (i=0;i<motifs.num_b4_rc;i++) {
    results[i] = malloc(sizeof(double)*num_seqs);
  }

  // Score each sequence.
  for (j=0;j<num_seqs;j++) {

    //copy the pointer in
    sequence = seq_list[j];
    //Store some sequence information for later
    add_string(get_seq_name(sequence),seq_ids);
    seq_len = get_seq_length(sequence);
    set_array_item(j,atof(get_seq_description(sequence)),seq_fscores);

    BOOLEAN_T needs_postprocessing; // ama_sequence_scan always sets this to TRUE
    // Score current sequence with each motif.
    for (i = 0;i<motifs.num_b4_rc;i++) {
      results[i][j] = 0;
      //Adjust motif_index for possible reverse complements.
      int adj_motif_index = use_rc ? i*2 : i;
      if (args.scoring == TOTAL_HITS) {
        //Configure the length correction, if any
        if (args.length_correction) {
          //Robert's backwards bonferroni (of course somebody else would have invented this)
          scaled_pvalue_threshold = rsc_bonferroni_correction(args.pvalue_threshold, 1.0/seq_len);
        }
        results[i][j] = ramen_sequence_scan(sequence, 
            motif_at(motifs.motifs, adj_motif_index),
            use_rc ? motif_at(motifs.motifs, adj_motif_index+1) : NULL,
            &motifs.pssms[adj_motif_index],
            use_rc ? &motifs.pssms[adj_motif_index+1] : NULL,
            NULL, NULL, //no need to pass odds
            args.scoring, 0, use_rc,
            scaled_pvalue_threshold, motifs.bg_freqs);
      } else { 
        results[i][j] = ramen_sequence_scan(sequence, 
            motif_at(motifs.motifs, adj_motif_index),
            use_rc ? motif_at(motifs.motifs, adj_motif_index+1) : NULL,
            NULL, NULL, //no need to pass pssm	
            motifs.odds[adj_motif_index],
            use_rc ? motifs.odds[adj_motif_index+1] : NULL,
            args.scoring, 0, use_rc,
            0, motifs.bg_freqs);

        if (TRUE == args.linreg_normalise) {
          int k, sym, len;
          double maxscore;
          MATRIX_T *odds;

          odds = motifs.odds[adj_motif_index];
          len = get_num_rows(odds);
          for (k = 0, maxscore = 1; k < len; k++) {
            double maxp = 0;
            for (sym = 0; sym < alph_size_core(alph); sym++) {
              if (maxp < get_matrix_cell(k, sym, odds)) {
                maxp = get_matrix_cell(k, sym, odds);
              }
            }
            maxscore *= maxp;
          }
          results[i][j] /= maxscore;
          if (j == 0 && args.verbose >= HIGHER_VERBOSE)
            fprintf(stderr, "\nScaling %s down by 1 / %.4g", 
                get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), maxscore);
        } // linreg normalize
      } // not totalhits
    } // motif
  } // sequence 

  // Write the absolute min, max and maximum of median over motifs of sequence scores in the HTML.
  double min_sequence_score = +1e300;
  double max_sequence_score = -1e300;
  double max_median_sequence_score = -1e300;
  double *all_scores = malloc(sizeof(double) * num_seqs); 
  for (i=0;i<motifs.num_b4_rc;i++) {
    for (j=0;j<num_seqs;j++) {
      if (results[i][j] < min_sequence_score) min_sequence_score = results[i][j];
      if (results[i][j] > max_sequence_score) max_sequence_score = results[i][j];
      all_scores[j] = results[i][j]; 
    } // i
    // compute the median for this motif and save the largest
    qsort(all_scores, num_seqs, sizeof(double), rsc_compare_doubles);
    double median = all_scores[num_seqs/2];
    if (median > max_median_sequence_score) max_median_sequence_score = median;
    //fprintf(stderr, "motif %d %s median %.2e\n", i, get_motif_id(motif_at(motifs.motifs, 2*i)), median);
  } // j
  //fprintf(stderr, "maximum median %.2e\n", max_median_sequence_score);
  free(all_scores);

  jsonwr_dbl_prop(args.json_output, "min_sequence_score", min_sequence_score);
  jsonwr_dbl_prop(args.json_output, "max_sequence_score", max_sequence_score);
  jsonwr_dbl_prop(args.json_output, "max_median_sequence_score", max_median_sequence_score);
} // rsc_scan_sequences

// Do significance testing using scores in Global results[][]
void rsc_get_scores(int argc, char *argv[]) {
  int i;
  int j;
  double sorted_scores[get_num_strings(seq_ids)];
  int seq_num;
  rsc_result_t** rsr;
  rsc_result_t* result;
  rsc_rank_t** rankings;

  seq_num = get_num_strings(seq_ids); //number of sequences
  //allocate space for final one result per motif array.
  rsr = malloc(sizeof(rsc_result_t*)*motifs.num_b4_rc);

  for(i=0;i<motifs.num_b4_rc;i++) {
    double highest_score = 0; //We use this to make sure that a PWM has scored on at least 1 seq.

    //Create the list of ranks
    rankings = malloc(sizeof(rsc_rank_t*)*seq_num);
    for (j=0; j < seq_num; j++) {
      rankings[j] = malloc(sizeof(rsc_rank_t));
      rankings[j]->pwm_score = results[i][j];
      if (results[i][j] > highest_score)
        highest_score = results[i][j]; //keep track of the PWM high score.
      rankings[j]->f_score = get_array_item(j, seq_fscores); //we don't  use this yet
      rankings[j]->f_rank = j+1;	// Ranks according to input order.
    }

    //Try shuffling here.
    // This seems to be to handle the case of tied scores
    // which are ordered arbitrarily by the call to qsort.
    // Ideally ties should be handled correctly in, for example,
    // the ranksum test.
    SHUFFLE(rankings, seq_num);

    //Sort it (ties are handled arbitrarily.
    qsort (rankings, seq_num, sizeof(rsc_rank_t*), rsc_compare_ranks_pwm_score);

    //Assign pwm ranks
    for (j=0; j < seq_num; j++) {
      rankings[j]->pwm_rank = j+1;	// Ranks when sorted by PWM.
    }

    //Positive FL is the default.
    if (POS_FL == args.positive_list) {
      //resort via FASTA score rank; actually rank is just original order of sequences.
      qsort (rankings, seq_num, sizeof(rsc_rank_t*), rsc_compare_ranks_f_rank);
    }

    ///Debugging code.
    if (args.verbose >= HIGHER_VERBOSE) {
      fprintf(stdout, "\n");
      for (j=0; j < seq_num; j++) {
        fprintf(stdout, "M2: %s - Seq: %s Rankings[%i] -\tpwm: %.8f\tprank: %i\tf: %.8f\tfrank: %i\n", 
            get_motif_st_id(motif_at(motifs.motifs,i*2)), get_nth_string(j,seq_ids), j,
            rankings[j]->pwm_score, rankings[j]->pwm_rank,
            rankings[j]->f_score, rankings[j]->f_rank);
      }
    }

    if (highest_score > 0) {
      //Determine the lowest subset pvalue
      if (RANKSUM_METHOD == args.pvalue_method) {
        rsr[i] = rsc_do_ranksum_test(rankings, i);
      } else if (FISHER_METHOD == args.pvalue_method) {
        rsr[i] = rsc_do_fisher_test(rankings, i);
      } else if (MULTIHG_METHOD == args.pvalue_method || 
          LONG_MULTIHG_METHOD == args.pvalue_method) {
        rsr[i] = rsc_do_multihg_test(rankings, i);
      } else if (LINREG_METHOD == args.pvalue_method) {
        rsr[i] = rsc_do_linreg_test(rankings,i);
      } else if (SPEARMAN_METHOD == args.pvalue_method) {
        rsr[i] = rsc_do_spearman_test(rankings,i);
      }
    } else {
      //If no sequence has scored at all, then we give a null result.
      // Bugfix, rsr[i] was not assigned the null result.
      rsr[i] = init_result (NULL, motifs.db_idx[i], motif_at(motifs.motifs, i), seq_num, 1, 1, 1, -1);
    }

    //Free up some space - TODO: Move this into a more appropriate place
    for (j=0; j < seq_num - 1; j++) {
      free(rankings[j]);
    }
    free(rankings);
  }

  if (LINREG_METHOD == args.pvalue_method) {
    qsort(rsr, motifs.num_b4_rc, sizeof(rsc_result_t*), rsc_compare_mse);
  } else {
    qsort(rsr, motifs.num_b4_rc, sizeof(rsc_result_t*), rsc_compare_scores);
  }

  if (LINREG_METHOD != args.pvalue_method && SPEARMAN_METHOD != args.pvalue_method) {
    int thresh = 1;
    if (args.fix_partition <= 0) thresh = seq_num;
    rsc_print_results("Motif p-values are corrected by #Motifs * "
        "#ThresholdsTested - (%i * %i = %i)\n\n",
        motifs.num_b4_rc, thresh, motifs.num_b4_rc * thresh);  
  } 
  jsonwr_property(args.json_output, "motifs");
  jsonwr_start_array_value(args.json_output);
  for (i = 0; i < motifs.num_b4_rc; ++i) {
    double number_of_tests;
    double corrected_p_val = 0;

    result = rsr[i];

    if (args.fix_partition > 0) {
      number_of_tests = motifs.num_b4_rc * 1;
    } else {
      number_of_tests = motifs.num_b4_rc * seq_num;
    }

    if (RANKSUM_METHOD == args.pvalue_method) {
      //corrected_p_val = rsc_bonferroni_correction(result->pboth, number_of_tests);
      corrected_p_val = rsc_bonferroni_correction(result->pright, number_of_tests);
      if (corrected_p_val <= args.pvalue_report_threshold) {
        rsc_print_results("%i. Ranksum p-values of motif %s %s (%s) top %i seqs "
            "(left,right,twotailed): %.4g %.4g %.4g U-value: %.4g (Corrected "
            "p-values: %.4g %.4g %.4g)\n",
            i+1, 
            get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif), 
            result->split, result->pleft, 
            result->pright, result->pboth, result->u, 
            rsc_bonferroni_correction(result->pleft,number_of_tests), 
            rsc_bonferroni_correction(result->pright,number_of_tests), 
            corrected_p_val);
        record_result(result, corrected_p_val); 
      }
    } else if (FISHER_METHOD == args.pvalue_method) {
      corrected_p_val = rsc_bonferroni_correction(result->pright,number_of_tests);
      if (corrected_p_val <= args.pvalue_report_threshold) {
        rsc_print_results("%i. Fisher's exact test p-value of motif %s %s (%s) top %i seqs: %.4g (Corrected p-value: %.4g)\n",
            i+1, 
            get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif),
	    result->split, result->pright, corrected_p_val);
        record_result(result, corrected_p_val);
      }
    } else if (MULTIHG_METHOD == args.pvalue_method) {
      corrected_p_val = rsc_bonferroni_correction(result->pright,number_of_tests);
      if (corrected_p_val <= args.pvalue_report_threshold) {
        rsc_print_results("%i. MultiHG p-value of motif %s %s (%s) top %i seqs: %.4g (Corrected p-value: %.4g)\n",
            i+1, 
            get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif),
            result->split, result->pright, corrected_p_val);
        record_result(result, corrected_p_val);
      }
    } else if (LONG_MULTIHG_METHOD == args.pvalue_method) {
      corrected_p_val = rsc_bonferroni_correction(result->pright,number_of_tests);
      if (corrected_p_val <= args.pvalue_report_threshold) {
        rsc_print_results("%i. 4D-MultiHG p-value of motif %s %s (%s) top %i seqs: %.4g "
            "(Corrected p-value: %.4g)\n", i+1, 
            get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif),
            result->split, result->pright, corrected_p_val);
        record_result(result, corrected_p_val);
      }
    } else if (LINREG_METHOD == args.pvalue_method) {
      rsc_print_results("%i. LinReg MSE of motif %s %s (%s) top %i seqs: %.4e (m: %.4e b: %.4e)\n",
          i+1, 
          get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif),
          result->split, result->pboth, result->pleft, result->pright);
      record_result(result, 0.0);
    } else if (SPEARMAN_METHOD == args.pvalue_method) {
      rsc_print_results("%i. Spearman Rho of motif %s %s (%s) top %i seqs: %.4e\n",
          i+1, 
          get_motif_id(result->motif), get_motif_id2(result->motif), get_motif_consensus(result->motif),
          result->split, result->pboth);
      record_result(result, 0.0);
    }
  }
  jsonwr_end_array_value(args.json_output);
}

/*
 * Using the spearman rank correlation test,
 *
 */
rsc_result_t* rsc_do_spearman_test(rsc_rank_t** rankings, int motif_index) {
  //Assorted vars
  int seq_num;
  int j;

  //Adjust motif_index for possible reverse complements.
  int adj_motif_index = motifs.revcomp ? motif_index*2 : motif_index;

  //Vars for the test
  double* x;
  double* y;

  //Vars for scoring
  double rank_score;
  rsc_result_t* lowest_motif_result;

  //Allocate memory or set initial values
  seq_num = get_num_strings(seq_ids); //number of sequences
  lowest_motif_result = malloc(sizeof(rsc_result_t)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.
  x = malloc(sizeof(double)*seq_num);
  y = malloc(sizeof(double)*seq_num);

  //Now we need to copy the scores into two double arrays
  //Use LOG macro so that log(0) 'works'
  for (j=0; j < seq_num; j++) {
    if (args.log_fscores == TRUE) {
      x[j] = LOG(rankings[j]->f_score);
    } else {
      x[j] = rankings[j]->f_score;
    }

    if (args.log_pwmscores == TRUE) {
      y[j] = LOG(rankings[j]->pwm_score);
    } else {
      y[j] = rankings[j]->pwm_score;
    }
    if (args.verbose == 5) {
      fprintf(stderr, "Rank %i: LR F-Score %.4g (%.4g) LR PWM-Score: %.4g "
          "(%.4g)\n", j, x[j], rankings[j]->f_score,y[j], rankings[j]->pwm_score);
    }
  }

  //Get the result
  rank_score = spearman_rank_correlation(seq_num,x,y);

  if (args.verbose >= HIGH_VERBOSE) {
    rsc_print_results("Spearman MSE of motif %s top %i seqs: %.4g\n", 
        get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), seq_num, 
        rank_score);
  }

  init_result (lowest_motif_result, motifs.db_idx[adj_motif_index], motif_at(motifs.motifs, adj_motif_index),
      seq_num,
      rank_score, rank_score,  //Not really p-values, but they'll do...
      rank_score, -1);

  return lowest_motif_result;
}
/*
 * Using the linreg test,
 *
 * this method returns the lowest scoring subdivision of a set of sequences for a given motif.
 * It's not self-contained, as it requires to hook into the global variables results, motifs, seq_ids.
 */
rsc_result_t* rsc_do_linreg_test(rsc_rank_t** rankings, int motif_index) {
  //Assorted vars
  int seq_num;
  int j;

  //Adjust motif_index for possible reverse complements.
  int adj_motif_index = motifs.revcomp ? motif_index*2 : motif_index;

  //Vars for the regression
  double* x;
  double* y;
  double m = 0;
  double b = 0;

  //Vars for scoring
  double lowest_mse = 10000;
  rsc_result_t* lowest_motif_result;

  //Allocate memory or set initial values
  seq_num = get_num_strings(seq_ids); //number of sequences
  lowest_motif_result = malloc(sizeof(rsc_result_t)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.
  x = malloc(sizeof(double)*seq_num);
  y = malloc(sizeof(double)*seq_num);

  //Now we need to copy the scores into two double arrays
  //Use LOG macro so that log(0) 'works'
  for (j=0; j < seq_num; j++) {
    if (args.log_fscores == TRUE) {
      x[j] = LOG(rankings[j]->f_score);
    } else {
      x[j] = rankings[j]->f_score;
    }

    if (args.log_pwmscores == TRUE) {
      y[j] = LOG(rankings[j]->pwm_score);
    } else {
      y[j] = rankings[j]->pwm_score;
    }
    if (args.verbose == 5) {
      fprintf(stderr, "Rank %i: LR F-Score %.4g (%.4g) LR PWM-Score: %.4g "
          "(%.4g)\n", j, x[j], rankings[j]->f_score,y[j], 
          rankings[j]->pwm_score);
    }
  }

  // TODO: Remove this for production
  if(args.linreg_dump_dir != NULL) {
    FILE *fh;
    char* filename;
    filename = malloc(sizeof(char)*(strlen(args.linreg_dump_dir) + 50));
    sprintf(filename, "%s/%s.tsv", args.linreg_dump_dir, 
        get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)));

    fh = fopen(filename, "w");
    fputs("X\tY\n", fh);
    for (j=0; j < seq_num; j++) {
      fprintf(fh, "%.10e %.10e\n", x[j], y[j]);
    }
    fclose(fh);
    free(filename);
  }


  //We start with a minimum of three sequences so that the data
  //is over-described.
  int min = 3;
  int max = seq_num;
  if (args.fix_partition) {min = args.fix_partition-1; max = min+1;}
  if (min < 3) min = 3;
  if (max > seq_num) max = seq_num;
  for (j=min; j < max; j++) {

    double mse = 0;
    if (args.linreg_switchxy) {
      mse = regress(j+1, y, x, &m, &b);
    } else {
      mse = regress(j+1, x, y, &m, &b);
    }

    // fix NANs
    if (isnan(mse)) {
      mse = 9999;	// smaller than initial value to make sure a value is stored in motif name
      m = b = 0;
    }

    if (args.verbose >= HIGH_VERBOSE) {
      rsc_print_results("LinReg MSE of motif %s top %i seqs: %.4g (m: %.4g b: %.4g)\n",
        get_motif_st_id(motif_at(motifs.motifs,adj_motif_index)), j+1, mse, m, b);
    }

    //Add to our motif list if lowest MSE
    if (lowest_mse > mse) {
      lowest_mse = mse;
      init_result (lowest_motif_result, motifs.db_idx[adj_motif_index], motif_at(motifs.motifs,adj_motif_index), j+1,
          m, b,  //Not really p-values, but they'll do...
          mse, -1);
    }

  }

  return lowest_motif_result;
}

/*
 * Using the Ranksum test,
 *
 * This method returns the lowest scoring subdivision of a set of sequences for a given motif.
 * It's not self-contained, as it requires to hook into the global variables results, motifs, seq_ids.
 */
rsc_result_t* rsc_do_ranksum_test(rsc_rank_t** rankings, int motif_index) {
  int n,na; // for rs stats test.
  double ta_obs; // for rs stats test.
  int seq_num;
  int j;
  double lowest_pval;
  RSR_T* r;
  rsc_result_t* lowest_motif_result;

  //Adjust motif_index for possible reverse complements.
  int adj_motif_index = motifs.revcomp ? motif_index*2 : motif_index;

  seq_num = get_num_strings(seq_ids); //number of sequences
  n = seq_num;
  na = ta_obs = 0;
  lowest_pval = 1; //1 is highest possible pval
  lowest_motif_result = malloc(sizeof(rsc_result_t)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.

  // By default "rankings" is sorted in the original order of the sequences.
  // The top of the rankings list are considered "positives".
  for (j=0; j < seq_num - 1; j++) {
    /* We now are not going to compare about ties for right now - instead we'll just
     * count them up.
     * We need to keep track of:
     *      int n,          number of samples
            int na,         number of positives
            double ta_obs   sum of ranks of positives
     */
    na++;
    if (POS_PWM == args.positive_list) {
      ta_obs += rankings[j]->f_rank; //We consider FASTA score to tell us what's positive
    } else {
      ta_obs += rankings[j]->pwm_rank; //We consider the PWM to tell us what's positive
    }

    // Compute p-value if at allowed partition point
    // Fixed bug: was stopping at j==args.fix_partition.
    if (args.fix_partition <= 0 || j == args.fix_partition-1) {
      if (QUICK_RS == args.rs_method) {
	//This thing is working backwards.
	// Should be: r = ranksum_from_stats(n, na, ta_obs);
	// Note: 1.0 is needed to force float arithmetic to avoid overflow
	r = ranksum_from_stats(n, n - na, n*(n+1.0)/2 - ta_obs);
      } else {
	// Beware 2d array pointer arithmetic.
	// This needs to be fixed badly.
	r = ranksum_from_sets(results[adj_motif_index] + (j+1), seq_num-(j+1), 
	    results[adj_motif_index], j+1);
      }

      if (args.verbose >= HIGH_VERBOSE) {
	fprintf(stderr, "Ranksum p-values of motif %s top %i seqs "
	    "(left,right,twotailed): %g %g %g U-value: %.4g \n",
	    get_motif_st_id(motif_at(motifs.motifs,adj_motif_index)), j+1, 
	    RSR_get_p_left(r), RSR_get_p_right(r), RSR_get_p_twotailed(r), RSR_get_u(r));
      }
      //Add to our motif list
      if (lowest_pval >= RSR_get_p_right(r)) {
	lowest_pval = RSR_get_p_right(r);
	init_result (lowest_motif_result, motifs.db_idx[adj_motif_index], motif_at(motifs.motifs, adj_motif_index), j+1,
	    RSR_get_p_left(r), RSR_get_p_right(r), RSR_get_p_twotailed(r),
	    RSR_get_u(r));
      }
      // Exit loop if fixed partition
      if (args.fix_partition) break;  
    } // possible partition

  }

  return lowest_motif_result;
}

/*
 * Using Fisher's exact test that the number of positives is as great as observed (one-tailed).
 *
 * This method returns the lowest scoring subdivision of a set of sequences for a given motif.
 * It's not self-contained, as it requires to hook into the global variables results, motifs, seq_ids.
 */
rsc_result_t* rsc_do_fisher_test(rsc_rank_t** rankings, int motif_index) {
  int seq_num;
  int i,j;
  double lowest_pval;
  double p=1,pleft=1,pright=1,pboth=1;
  RSR_T* r;
  rsc_result_t* lowest_motif_result;
  int tp,fn,fp,tn;

  //Adjust motif_index for possible reverse complements.
  int adj_motif_index = motifs.revcomp ? motif_index*2 : motif_index;

  seq_num = get_num_strings(seq_ids); //number of sequences
  lowest_pval = 100; //a large number.
  lowest_motif_result = malloc(sizeof(rsc_result_t)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.

  //Used for testing so that we can be sure that we've initialised this variable
  lowest_motif_result->u = 0;

  //fisher_exact_init(seq_num); //initialise our factorial datastructure

  if (args.fix_partition > 0) {
    /*
     * We do just one partition.
     */
    int i = args.fix_partition - 1;
    tp = fn = fp = tn = 0;
    for (j=0;j< seq_num; j++) {
      if (j <= i) {
        //We're doing above the threshold
        if (POS_FL == args.positive_list) {
          //Well, we're sorted by FL score, so we check the PWM score, and if it's over
          //our threshold, then we add one to the tp.
          (rankings[j]->pwm_score >= args.fisher_pwm_threshold) ? tp++ : fn++;
        } else {
          //Well, we're sorted by PWM score, so we check the FL score, and if it's under
          //our threshold, then we add one to the tp.
          (rankings[j]->f_score <= args.fisher_pwm_threshold) ? tp++ : fn++;
        }
      } else {
        //We're doing below the threshold
        if (POS_FL == args.positive_list) {
          //Well, we're sorted by FL score, so we check the PWM score, and if it's over
          //our threshold, then we add one to the tn.
          (rankings[j]->pwm_score < args.fisher_pwm_threshold) ? tn++ : fp++;
        } else {
          //Well, we're sorted by PWM score, so we check the FL score, and if it's under
          //our threshold, then we add one to the tn.
          (rankings[j]->f_score > args.fisher_fasta_threshold) ? tn++ : fp++;
        }
      }
    }

    pright = exp(getLogFETPvalue(tp, tp+fn, fp, tn+fp, 0));
    lowest_pval = pright;
    init_result (lowest_motif_result, motifs.db_idx[adj_motif_index], motif_at(motifs.motifs,adj_motif_index), 
        i+1, lowest_pval, lowest_pval, lowest_pval, -1);
    if (args.verbose >= HIGH_VERBOSE) {
      fprintf(stderr, "M3: %s Threshold: %i tp: %i fn: %i tn: %i fp: %i P %g\n",
	  get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), 
	  i, tp,fn,tn,fp, pright);
    }
  } else {
    /*
     * We do partition maximization
     */
    for (i=0; i < seq_num - 1; i++) {
      //TODO: Make this not n^2
      tp = fn = fp = tn = 0;
      for (j=0;j< seq_num; j++) {
        if (j <= i) {
          //We're doing above the threshold
          if (POS_FL == args.positive_list) {
            //Well, we're sorted by FL score, so we check the PWM score, and if it's over
            //our threshold, then we add one to the tp.
            (rankings[j]->pwm_score >= args.fisher_pwm_threshold) ? tp++ : fn++;
          } else {
            //Well, we're sorted by PWM score, so we check the FL score, and if it's under
            //our threshold, then we add one to the tp.
            (rankings[j]->f_score <= args.fisher_pwm_threshold) ? tp++ : fn++;
          }
        } else {
          //We're doing below the threshold
          if (POS_FL == args.positive_list) {
            //Well, we're sorted by FL score, so we check the PWM score, and if it's over
            //our threshold, then we add one to the tn.
            (rankings[j]->pwm_score < args.fisher_pwm_threshold) ? tn++ : fp++;
          } else {
            //Well, we're sorted by PWM score, so we check the FL score, and if it's under
            //our threshold, then we add one to the tn.
            (rankings[j]->f_score > args.fisher_fasta_threshold) ? tn++ : fp++;
          }
        }
      }

      pright = exp(getLogFETPvalue(tp, tp+fn, fp, tn+fp, 0));
      if (args.verbose >= HIGH_VERBOSE) {
        fprintf(stderr, "Fisher's exact test p-value of motif %s top %i seqs (tp >= observed): %g\n",
            get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), i+1, pright);
      }

      //Add to our motif list
      if (lowest_pval >= pright) {
        lowest_pval = pright;
        init_result (lowest_motif_result, motifs.db_idx[adj_motif_index], motif_at(motifs.motifs,adj_motif_index),
            i+1, lowest_pval, lowest_pval, lowest_pval, -1);
      }
      if (args.verbose >= HIGHER_VERBOSE) {
	fprintf(stderr, "M: %s Threshold: %i tp: %i fn: %i tn: %i fp: %i P %g\n",
	    get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), 
	    i, tp,fn,tn,fp, pright);
      }
    }

  }

  assert(lowest_motif_result != NULL);
  assert(lowest_motif_result->u != 0);
  return lowest_motif_result;
}

/*
 * Using the Ranksum test,
 *
 * this method returns the lowest scoring subdivision of a set of sequences for a given motif.
 * It's not self-contained, as it requires to hook into the global variables results, motifs, seq_ids.
 */
rsc_result_t* rsc_do_multihg_test(rsc_rank_t** rankings, int motif_index) {
  int seq_num;
  int i,j,k;
  double lowest_pval;
  double p;
  RSR_T* r;
  rsc_result_t* lowest_motif_result;
  int b0,b1,b2,b3;
  int n,N;
  int i0,i1,i2,i3,B0,B1,B2,B3;

  //Adjust motif_index for possible reverse complements.
  int adj_motif_index = motifs.revcomp ? motif_index*2 : motif_index;

  seq_num = get_num_strings(seq_ids); //number of sequences
  lowest_pval = 100; //a large number.
  lowest_motif_result = malloc(sizeof(rsc_result_t)); //allocate space, as a ptr to this will go in the array later
  //that's why we don't free it in this loop.

  //Used for testing so that we can be sure that we're initialised this variable
  lowest_motif_result->u = 0;

  //Get per-class totals
  B0 = B1 = B2 = B3 = 0;
  for (i=0; i < seq_num; i++) {
    if (rankings[i]->pwm_score == 0) {
      B0++;
    } else if (rankings[i]->pwm_score == 1) {
      B1++;
    } else if (rankings[i]->pwm_score == 2) {
      B2++;
    } else {
      B3++;
    }
  }
  // If we are only 3d rather than 4d...
  if (args.pvalue_method == MULTIHG_METHOD) {
    B2+=B3;
  }

  /*
   * Do it for a fixed partition only.
   */
  if (args.fix_partition > 0) {
    int i = args.fix_partition - 1;
    /* Need to do count table:
     *
     *     a0 - above threshold, 0 hits
     *     a1 - above threshold, 1 motif hit
     *     a2 - above threshold, >=2 motif hits
     *
     */
    b0 = b1 = b2 = b3 = 0;
    for (j=0;j<=i; j++) {
      if (rankings[j]->pwm_score == 0) {
        b0++;
      } else if (rankings[j]->pwm_score == 1) {
        b1++;
      } else if (rankings[j]->pwm_score == 2) {
        b2++;
      } else {
        b3++;
      }
    }

    if (args.pvalue_method == MULTIHG_METHOD) {

      b2 += b3;
      /* Formula details:
       *
       * Please see page 31 of Robert's lab book.
       *
       *
       * The notation is the same as Eden et al. 2007.
       *
       * If we count two or more. (3>2, so we add it to B2 etc), since this
       * piece of code has a max dimension of 3.
       *
       */
      n = i+1; //n is the threshold
      N = seq_num; //total set size;
      p = LOGZERO;

      for (i0=b0; i0<=B0 && i0<=n; i0++) {
        for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
          for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
            //fprintf(stderr, "b0: %i b1: %i b2: %i B0: %i B1: %i B2: %i n: %i N: %i\n",
            //                                           b0,b1,b2,B0,B1,B2, n, N);
            p = LOG_SUM(p, //We're in log space, remember.
                (
                 (
                  factorial_array[n] - (factorial_array[i0] + 
                    factorial_array[i1] + factorial_array[i2])
                 ) + (
                   (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                     factorial_array[B1-i1] + factorial_array[B2-i2])
                   )
                ) - (
                  (factorial_array[N] - (factorial_array[B0] + 
                      factorial_array[B1] + factorial_array[B2]))
                  )
                );
            //fprintf(stderr, "Current p for i0: %i i1 %i i2 %i: %g\n",i0,i1,i2,exp(p));
          }
        }
      }

    } else { // We're going to do the complicated method with four loops.

      /* Formula details:
       *
       * Please see page 31 of Robert's lab book.
       *
       *
       * The notation is the same as Eden et al. 2007.
       *
       *
       */
      n = i+1; //n is the threshold
      N = seq_num; //total set size;
      p = LOGZERO;

      for (i0=b0; i0<=B0 && i0<=n; i0++) {
        for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
          for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
            for (i3=b3; i3<=B3 && i3<=n-(i0+i1+i2); i3++) {
              p = LOG_SUM(p, //We're in log space, remember.
                  (
                   (
                    factorial_array[n] - (factorial_array[i0] + 
                      factorial_array[i1] + factorial_array[i2] + 
                      factorial_array[i3])
                   ) + (
                     (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                       factorial_array[B1-i1] + factorial_array[B2-i2] + 
                       factorial_array[B3-i3])
                     )
                  ) - (
                    (factorial_array[N] - (factorial_array[B0] + 
                            factorial_array[B1] + factorial_array[B2]) + 
                     factorial_array[B3])
                    )
                  );
              //fprintf(stderr, "Current p for i0: %i i1 %i i2 %i: %g\n",
              //  i0,i1,i2,exp(p));
            }
          }
        }
      }

    }


    //Add to our motif list
    if (args.verbose >= HIGH_VERBOSE) {
      fprintf(stderr, "Motif: %s Threshold: %i P-value: %g\n", 
          get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), i, exp(p));
    }
    lowest_pval = p;
    init_result (lowest_motif_result, motifs.db_idx[adj_motif_index],
        motif_at(motifs.motifs, adj_motif_index), i+1,
        lowest_motif_result->pleft, // exp(p)
        p,
        lowest_motif_result->pboth, // exp(p)
        -1);
  } else {

    //TODO: Make this not n^2
    for (i=0; i < seq_num - 1; i++) {
      /* Need to do count table:
       *
       *     a0 - above threshold, 0 hits
       *     a1 - above threshold, 1 motif hit
       *     a2 - above threshold, >=2 motif hits
       *
       */
      b0 = b1 = b2 = b3 = 0;
      for (j=0;j<=i; j++) {
        if (rankings[j]->pwm_score == 0) {
          b0++;
        } else if (rankings[j]->pwm_score == 1) {
          b1++;
        } else if (rankings[j]->pwm_score == 2) {
          b2++;
        } else {
          b3++;
        }
      }

      if (args.pvalue_method == MULTIHG_METHOD) {

        b2 += b3;
        /* Formula details:
         *
         * Please see page 31 of Robert's lab book.
         *
         *
         * The notation is the same as Eden et al. 2007.
         *
         * If we count two or more. (3>2, so we add it to B2 etc), since this
         * piece of code has a max dimension of 3.
         *
         */
        n = i+1; //n is the threshold
        N = seq_num; //total set size;
        p = LOGZERO;

        for (i0=b0; i0<=B0 && i0<=n; i0++) {
          for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
            for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
              //fprintf(stderr, "b0: %i b1: %i b2: %i B0: %i B1: %i B2: %i n: %i N: %i\n",
              //                                           b0,b1,b2,B0,B1,B2, n, N);
              p = LOG_SUM(p, //We're in log space, remember.
                  (
                   (
                    factorial_array[n] - (factorial_array[i0] + 
                      factorial_array[i1] + factorial_array[i2])
                   ) + (
                     (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                       factorial_array[B1-i1] + factorial_array[B2-i2])
                     )
                  ) - (
                    (factorial_array[N] - (factorial_array[B0] + 
                        factorial_array[B1] + factorial_array[B2]))
                    )
                  );
              //fprintf(stderr, "Current p for i0: %i i1 %i i2 %i: %g\n",
              //  i0,i1,i2,exp(p));
            }
          }
        }

      } else { // We're going to do the complicated method with four loops.

        /* Formula details:
         *
         * Please see page 31 of Robert's lab book.
         *
         *
         * The notation is the same as Eden et al. 2007.
         *
         *
         */
        n = i+1; //n is the threshold
        N = seq_num; //total set size;
        p = LOGZERO;

        for (i0=b0; i0<=B0 && i0<=n; i0++) {
          for (i1=b1; i1<=B1 && i1<=n-i0; i1++) {
            for (i2=b2; i2<=B2 && i2<=n-(i0+i1); i2++) {
              for (i3=b3; i3<=B3 && i3<=n-(i0+i1+i2); i3++) {
                p = LOG_SUM(p, //We're in log space, remember.
                    (
                     (
                      factorial_array[n] - (factorial_array[i0] + 
                        factorial_array[i1] + factorial_array[i2] + 
                        factorial_array[i3])
                     ) + (
                       (factorial_array[N-n]) - (factorial_array[B0-i0] + 
                         factorial_array[B1-i1] + factorial_array[B2-i2] + 
                         factorial_array[B3-i3])
                       )
                    ) - (
                      (factorial_array[N] - (factorial_array[B0] + 
                              factorial_array[B1] + factorial_array[B2]) + 
                       factorial_array[B3])
                      )
                    );
                //fprintf(stderr, "Current p for i0: %i i1 %i i2 %i: %g\n",
                //    i0,i1,i2,exp(p));
              }
            }
          }
        }

      }


      //TODO :fix below here.
      //Add to our motif list
      if (args.verbose >= HIGH_VERBOSE) {
        fprintf(stderr, "Motif: %s Threshold: %i P-value: %g\n", 
            get_motif_st_id(motif_at(motifs.motifs, adj_motif_index)), i, exp(p));
      }
      if (lowest_pval >= p) {
        lowest_pval = p;
        init_result (lowest_motif_result, motifs.db_idx[adj_motif_index],
            motif_at(motifs.motifs, adj_motif_index),
            i+1,
            lowest_motif_result->pleft, // exp(p)
            p,
            lowest_motif_result->pboth, // exp(p)
            -1);
      }

    }
  }

  assert(lowest_motif_result != NULL);
  assert(lowest_motif_result->u != 0);
  lowest_motif_result->pright = exp(lowest_motif_result->pright); //Small increase in speed.
  return lowest_motif_result;
}

void rsc_terminate(int status) {
  /* Display time of execution */
  if (verbosity >= HIGH_VERBOSE) {
    t1 = time(NULL);
    c1 = clock();
    fprintf (stderr, "Elapsed wall clock time: %ld seconds\n", (long)(t1 - t0));
    fprintf (stderr, "Elapsed CPU time:        %f seconds\n", 
        (float)(c1 - c0) / CLOCKS_PER_SEC);
  }
  // don't risk closing the files if something went wrong: 
  //  should be more specific on when not to do this
  if (!status && args.sequence_filename) final_print_results();
  exit(status);
}

/*
 * COMPARISON METHODS FOLLOW
 */

int rsc_compare_scores (const void *a, const void *b)
{
  rsc_result_t r1 = **(rsc_result_t**)a;
  rsc_result_t r2 = **(rsc_result_t**)b;
  if (r2.pright - r1.pright < 0.0) {
    return 1;
  } else if (r2.pright - r1.pright > 0.0){
    return -1;
  } else {
    return 0;
  }
  //This compares it in reverse order
  //  return (int) (r2.pright - r1.pright);
}

int rsc_compare_mse (const void *a, const void *b)
{
  rsc_result_t r1 = **(rsc_result_t**)a;
  rsc_result_t r2 = **(rsc_result_t**)b;
  if (r2.pboth - r1.pboth < 0.0) {
    return 1;
  } else if (r2.pboth - r1.pboth > 0.0){
    return -1;
  } else {
    return 0;
  }
  //This compares it in reverse order
  //  return (int) (r2.pright - r1.pright);
}

int rsc_compare_doubles (const void *a, const void *b)
{
  //This does it in reverse order
  return (int) (*(double *)b - *(double *)a);
}

int rsc_compare_ranks_f_rank (const void *a, const void *b) {
  rsc_rank_t r1 = **(rsc_rank_t**)a;
  rsc_rank_t r2 = **(rsc_rank_t**)b;

  //fprintf(stderr, "Comparing: %i to %i\n", r1.f_rank, r2.f_rank);

  if (r2.f_rank - r1.f_rank < 0.0) {
    return 1;
  } else if (r2.f_rank - r1.f_rank > 0.0){
    return -1;
  } else {
    return 0;
  }
  //This compares it smallest first

}

int rsc_compare_ranks_pwm_score (const void *a, const void *b) {
  rsc_rank_t r1 = **(rsc_rank_t**)a;
  rsc_rank_t r2 = **(rsc_rank_t**)b;

  //fprintf(stderr, "\nComparing: %g to %g", r1.pwm_score, r2.pwm_score);

  if (r1.pwm_score - r2.pwm_score < 0.0) {
    return 1;
  } else if (r1.pwm_score - r2.pwm_score > 0.0){
    return -1;
  } else {
    return 0;
  }
  //This compares it largest first
}

/*
 * DEBUGGING METHODS FOLLOW.
 */

void rsc_dump_motif_freqs(MOTIF_T* m, MATRIX_T* freqs) {
  int i, j;
  ALPH_T *alph;

  alph = get_motif_alph(m);
  if (freqs == NULL) freqs = get_motif_freqs(m);
  // print alphabet headers
  rsc_print_results("%s", get_motif_st_id(m));
  for (j = 0; j < alph_size_core(alph); j++) {
    rsc_print_results("\t%c", alph_char(alph, j));
  }
  rsc_print_results("\n");
  // print matrix freqs
  for (i = 0; i < get_num_rows(freqs); i++) {
    rsc_print_results("%s", get_motif_st_id(m));
    for (j = 0; j < alph_size_core(alph); j++) {
      rsc_print_results("\t%.4f", get_matrix_cell(i, j,freqs));
    }
    rsc_print_results("\n");
  }
}

unsigned long long choose(unsigned n, unsigned k) {
  unsigned i;
  long double accum;

  //fprintf(stderr, "%i choose %i: ",n,k);

  if (k > n)
    return 0;

  if (k > n/2)
    k = n-k; // faster

  accum = 1;
  for (i = 1; i <= k; i++)
    accum = accum * (n-k+i) / i;

  //fprintf(stderr, "%i\n", (int)(accum+0.5));
  return accum + 0.5; // avoid rounding error
}

long double logchoose(unsigned n, unsigned k) {
  unsigned i;
  long double accum;

  //fprintf(stderr, "%i choose %i: ",n,k);

  if (k > n)
    return 0;

  if (k > n/2)
    k = n-k; // faster

  accum = 0;
  for (i = 1; i <= k; i++)
    accum += log((n-k+i) / i);

  //fprintf(stderr, "%i\n", (int)(accum+0.5));
  return accum; // avoid rounding error
}

// if p too small, the power forumula rounds to 0
double rsc_bonferroni_correction(double p, double numtests) {
  /* double correction = 1 - pow(1-p, numtests);
     if (correction)
     return correction;
     else
     return p*numtests; */
  return exp(LOGEV(log(numtests), log(p)));
}

// copy a string into new memory of exactly the right size; create a "" string
// if the original is a NULL pointer
char * my_strdup(const char *s1) {
  if (s1) {
    return strdup (s1);
  } else {
    char * newstr = malloc(sizeof(char));
    newstr[0] = '\0';
    return newstr;
  }
}

// create a new instance of a result, with a given motif and initial values
// if a non-null pointer is passed in, use its value, otherwise pass overwrite
// the location pointed to by new_result; in either case return new_result
rsc_result_t * init_result (rsc_result_t *new_result,
    int db_idx,
    MOTIF_T *motif,
    int split,
    double pleft,
    double pright,
    double pboth,
    double u) {
  if (!new_result) new_result = malloc(sizeof(rsc_result_t));
  new_result->u = u;
  new_result->split = split;
  new_result->pleft = pleft;
  new_result->pright = pright;
  new_result->pboth = pboth;
  new_result->db_idx = db_idx;
  new_result->motif = duplicate_motif(motif);
  return new_result;
}


