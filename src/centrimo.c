/********************************************************************
 * FILE: centrimo.c
 * AUTHOR: Timothy Bailey, Philip Machanick
 * CREATE DATE: 06/06/2011
 * PROJECT: MEME suite
 * COPYRIGHT: 2011, UQ
 ********************************************************************/

#define DEFINE_GLOBALS
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"
#include "alphabet.h"
#include "binomial.h"
#include "config.h"
#include "dir.h"
#include "fasta-io.h"
#include "html-monolith.h"
#include "io.h"
#include "motif-in.h"
#include "projrel.h"
#include "pssm.h"
#include "seq-reader-from-fasta.h"
#include "string-list.h"
#include "utils.h"
#include "seq.h"
#include "fisher_exact.h"

#include <sys/resource.h>
#include <unistd.h>

const double DEFAULT_PSEUDOCOUNT = 0.1;
const double DEFAULT_SCORE_THRESH = 5.0;
const double DEFAULT_EVALUE_THRESH = 10.0;
const double ALPHA = 1.0; // Non-motif specific scale factor.
const char* TEMPLATE_FILENAME = "centrimo_template.html";
const char* SITES_FILENAME = "site_counts.txt"; // Name of plain-text sites
const char* TEXT_FILENAME = "centrimo.txt"; // Name of plain-text centrimo output
const char* HTML_FILENAME = "centrimo.html"; // Name of HTML centrimo output

VERBOSE_T verbosity = NORMAL_VERBOSE;
// output a message provided the verbosity is set appropriately
#define DEBUG_MSG(debug_level, debug_msg) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg); \
  } \
}

#define DEBUG_FMT(debug_level, debug_msg_format, ...) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg_format, __VA_ARGS__); \
  } \
}

// Structure for tracking centrimo command line parameters.
typedef struct options {
  ALPH_T* alphabet; // Alphabet

  bool allow_clobber; // Allow overwriting of files in output directory.
  bool scan_both_strands; // Scan forward and reverse strands.
  bool scan_separately; // Scan with reverse complement motif separately
  bool flip; // score given and RC strand with given motif
  bool local; // Compute the enrichment of all bins
  bool noseq; // Do not store sequences in HTML output
  bool neg_sequences; // Use a set of negative sequences
  bool disc; // Use Fisher's Exact Test
  bool mcc; // Calculate Matthews Corrlelation Coefficent
  bool optimize_score;
  RBTREE_T *selected_motifs; // IDs of requested motifs.

  char* description; // description of job
  char* desc_file; // file containing description
  char* alph_file; // Name of file containing alphabet definition (for converting motifs)
  char* bg_source; // Name of file file containing background freq.
  char* output_dirname; // Name of the output directory
  char* seq_source; // Name of file containing sequences.
  char* negseq_source; // Name of file containing sequences.

  ARRAYLST_T* motif_sources; // Names of files containing motifs.
  int seq_len; // Only accept sequences with this length

  double score_thresh; // Minimum value of site score (bits).
  double pseudocount; // Pseudocount added to motif PSPM.
  double evalue_thresh; // Don't report results with worse evalue

  int min_win; // Minimum considered central window size
  int max_win; // Maximum considered central window size
} CENTRIMO_OPTIONS_T;

typedef struct motif_db {
  int id;
  char* source;
  ARRAYLST_T* motifs;
} MOTIF_DB_T;

typedef struct site {
  int start;
  char strand;
} SEQ_SITE_T;

typedef struct sites {
  double best;
  int allocated;
  int used;
  SEQ_SITE_T *sites;
} SEQ_SITES_T;

typedef struct best_site {
  double score;
  int sequence_number;
  int offset; // this is an array index
  int used; // number of items in array at above index
} BEST_SITE_T;

typedef struct best_sites {
  int *sites;
  int allocated_sites;
  int used_sites;
  BEST_SITE_T *scores;
  int allocated_scores;
  int used_scores;
} BEST_SITES_T;


typedef struct counts {
  long total_sites;
  int allocated;
  double *sites; // there are ((2 * seqlen) - 1) sites
} SITE_COUNTS_T;

typedef struct win_stats {
  // what is the minimum site score
  double site_score_threshold;
  // where is the window
  int center;
  int spread; // this value is >= 1 and odd
  // how many bins are in the window
  int n_bins;
  // how many sites appear in the window?
  double sites;
  double neg_sites;
  // how significant is the window? (possibly compared to the negatives)
  double log_pvalue;
  double log_adj_pvalue;
  // significance of the same window in the negative set if it exists
  double neg_log_pvalue;
  double neg_log_adj_pvalue;
  // Fisher's Exact test p-value comparing positive and negative enrichment
  double fisher_log_pvalue;
  double fisher_log_adj_pvalue;
  // how different is the negative set?
  double mcc; // Matthews correlation coefficient
} WIN_STATS_T;

typedef struct motif_stats {
  MOTIF_DB_T* db; // motif database
  MOTIF_T* motif;
  ARRAYLST_T *windows;
  RBTREE_T* seq_ids;
  // how many tests have been done?
  int n_tests;
  // how many bins could this motif land in
  int n_bins;
  // what was the scoring threshold used?
  double score_threshold;
  // how many sites are found for this motif?
  long sites;
  long neg_sites;
} MOTIF_STATS_T;


typedef struct buffer {
  SEQ_SITES_T* sites;
  BEST_SITES_T* pos_sites;
  BEST_SITES_T* neg_sites;
  SITE_COUNTS_T* pos_counts;
  SITE_COUNTS_T* neg_counts;
  SITE_COUNTS_T* pos_c_counts;
  SITE_COUNTS_T* neg_c_counts;
} CENTRIMO_BUFFER_T;

#define SCORE_BLOCK 20
#define BEST_SITE_BLOCK 100
#define BEST_SCORE_BLOCK 100

/*************************************************************************
 * Routines to destroy data
 *************************************************************************/
static void destroy_window(void *w) {
  WIN_STATS_T *win = (WIN_STATS_T*)w;
  memset(win, 0, sizeof(WIN_STATS_T));
  free(win);
}
static void destroy_stats(void *s) {
  MOTIF_STATS_T *stats;
  if (!s) return;
  stats = (MOTIF_STATS_T*)s;
  arraylst_destroy(destroy_window, stats->windows);
  if (stats->seq_ids) rbtree_destroy(stats->seq_ids);
  memset(stats, 0, sizeof(MOTIF_STATS_T));
  free(stats);
}

/*************************************************************************
 * Compare two window statistics in an arraylst
 * Takes pointers to pointers to WIN_STATS_T.
 *************************************************************************/
static int window_stats_compare_pvalue(const void *w1, const void *w2) {
  WIN_STATS_T *win1, *win2;
  double pv_diff;
  int spread_diff;
  win1 = *((WIN_STATS_T**) w1);
  win2 = *((WIN_STATS_T**) w2);
  // sort by p-value ascending
  pv_diff = win1->log_adj_pvalue - win2->log_adj_pvalue;
  if (pv_diff != 0) return (pv_diff < 0 ? -1 : 1);
  // sort by window size ascending
  spread_diff = win1->spread - win2->spread;
  if (spread_diff != 0) return spread_diff;
  // sort by position
  return win1->center - win2->center;
}

/*************************************************************************
 * Routines for item comparison
 *************************************************************************/
static int compare_best_site(const void *a, const void *b) {
  const BEST_SITE_T *a_site, *b_site;
  a_site = (BEST_SITE_T*)a;
  b_site = (BEST_SITE_T*)b;
  if (a_site->score < b_site->score) {
    return 1;
  } else if (b_site->score < a_site->score) {
    return -1;
  }
  if (a_site->sequence_number > b_site->sequence_number) {
    return 1;
  } else if (a_site->sequence_number < b_site->sequence_number) {
    return -1;
  }
  return 0;
}

/***********************************************************************
 Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options(CENTRIMO_OPTIONS_T *options) {
  alph_release(options->alphabet);
  rbtree_destroy(options->selected_motifs);
  arraylst_destroy(NULL, options->motif_sources);
}

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage = 
    "\n"
    "Usage: centrimo [options] <sequence file> <motif file>+\n"
    "\n"
    "   Options:\n"
    "     --o <output dir>         output directory; default: 'centrimo_out'\n"
    "     --oc <output dir>        allow overwriting; default: 'centrimo_out'\n"
    "     --neg <fasta file>       plot motif distributions in this set as well\n"
    "                               in the <sequence file> (positive sequences)\n"
    "     --disc                   use Fisher's exact test instead of the\n"
    "                               binomial test to compare the control and\n"
    "                               the treatment (requires --neg <fasta file>)\n"
    "     --xalph <alph file>     Name of the file containing a custom alphabet,\n"
    "                              which specifies that motifs should be converted;\n"
    "     --bfile <background>    background frequency model for PWMs;\n"
    "                               default: base frequencies in input sequences\n"
    "     --local                  compute the enrichment of all regions;\n"
    "                               default: enrichment of central regions only\n"
    "     --maxreg <maxreg>        maximum width of any region to consider;\n"
    "                               default: use the sequence length\n"
    "     --minreg <minreg>        minimum width of any region to consider;\n"
    "                               must be less than <maxreg>;\n"
    "                               default: 1 bp\n"
    "     --score <S>              score threshold for PWMs, in bits;\n"
    "                               sequences without a site with score >= <S>\n"
    "                               are ignored;\n"
    "                               default: %.1g\n"
    "     --optimize_score         search for optimized score above the threshold\n"
    "                               given by '--score' argument. Slow computation\n"
    "                               due to multiple tests.\n"
    "                               default: use fixed score threshold\n"
    "     --norc                   do not scan with the reverse complement motif\n"
    "     --sep                    scan separately with reverse complement motif;\n"
    "                               (requires --norc)\n"
    "     --flip                   'flip' sequences so that matches on the \n"
    "                               reverse strand are 'reflected' around center;\n"
    "                               default: do not flip sequences\n"
    "     --seqlen <length>        use sequences with this length; default: use\n"
    "                               sequences with the same length as the first\n"
    "     --motif <ID>             only scan with this motif; options may be\n"
    "                               repeated to specify more than one motif;\n"
    "                               default: scan with all motifs\n"
    "     --motif-pseudo <pseudo>  pseudo-count to use creating PWMs;\n"
    "                               default: %.3g\n"
    "     --ethresh <thresh>       evalue threshold for including in results;\n"
    "                               default: %g\n"
    "     --desc <description>     include the description in the output;\n"
    "                               default: no description\n"
    "     --dfile <desc file>      use the file content as the description;\n"
    "                               default: no descriptionn\n"
    "     --noseq                  do not store sequence IDs in HTML output;\n"
    "                               default: IDs are stored in the HTML output\n"
    "     --verbosity [1|2|3|4]    verbosity of output: 1 (quiet) - 4 (dump);\n"
    "                               default: %d\n"
    "     --version                print the version and exit\n"
    ;
  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, usage, DEFAULT_SCORE_THRESH, DEFAULT_PSEUDOCOUNT,
      DEFAULT_EVALUE_THRESH, NORMAL_VERBOSE);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {OPT_BFILE, OPT_O, OPT_OC, OPT_SEQLEN, OPT_SCORE, OPT_MOTIF,
  OPT_MOTIF_PSEUDO, OPT_ETHRESH, OPT_MINREG, OPT_MAXREG, OPT_NORC, OPT_SEP,
  OPT_FLIP, OPT_NOFLIP, OPT_DESC, OPT_DFILE, OPT_LOCAL, OPT_NOSEQ, OPT_NEG,
  OPT_DISC, OPT_VERBOSITY, OPT_OPTIMIZE_SCORE, OPT_VERSION, OPT_XALPH};

/***********************************************************************
 Process command line options
 ***********************************************************************/
static void process_command_line(int argc, char* argv[], CENTRIMO_OPTIONS_T *options) {
  struct option centrimo_options[] = {
    {"xalph", required_argument, NULL, OPT_XALPH},
    {"bfile", required_argument, NULL, OPT_BFILE},
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"seqlen", required_argument, NULL, OPT_SEQLEN},
    {"score", required_argument, NULL, OPT_SCORE},
    {"motif", required_argument, NULL, OPT_MOTIF},
    {"motif-pseudo", required_argument, NULL, OPT_MOTIF_PSEUDO},
    {"ethresh", required_argument, NULL, OPT_ETHRESH},
    {"minreg", required_argument, NULL, OPT_MINREG},
    {"maxreg", required_argument, NULL, OPT_MAXREG},
    {"norc", no_argument, NULL, OPT_NORC},
    {"sep", no_argument, NULL, OPT_SEP},
    {"flip", no_argument, NULL, OPT_FLIP},
    {"noflip", no_argument, NULL, OPT_NOFLIP},
    {"desc", required_argument, NULL, OPT_DESC},
    {"dfile", required_argument, NULL, OPT_DFILE},
    {"local", no_argument, NULL, OPT_LOCAL},
    {"noseq", no_argument, NULL, OPT_NOSEQ},
    {"neg", required_argument, NULL, OPT_NEG},
    {"disc", no_argument, NULL, OPT_DISC},
    {"optimize_score", no_argument, NULL, OPT_OPTIMIZE_SCORE},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };
  // Make sure options are set to defaults.
  memset(options, 0, sizeof(CENTRIMO_OPTIONS_T));
  options->alphabet = NULL;
  options->allow_clobber = true;
  options->scan_both_strands = true;
  options->scan_separately = false;
  options->flip = false;
  options->local = false;
  options->noseq = false;
  options->neg_sequences = false;
  options->disc = false;
  options->mcc = false;
  options->optimize_score = false;
  options->description = NULL;
  options->desc_file = NULL;
  options->alph_file = NULL;
  options->bg_source = NULL;
  options->output_dirname = "centrimo_out";
  options->seq_source = NULL;
  options->motif_sources = arraylst_create();
  options->seq_len = 0; // get the length from the first sequence
  options->score_thresh = DEFAULT_SCORE_THRESH;
  options->pseudocount = DEFAULT_PSEUDOCOUNT;
  options->evalue_thresh = DEFAULT_EVALUE_THRESH;
  options->min_win = 0;
  options->max_win = 0;
  options->selected_motifs = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
  verbosity = NORMAL_VERBOSE;
  // process arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", centrimo_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_XALPH:
        options->alph_file = optarg;
        break;
      case OPT_BFILE:
        options->bg_source = optarg;
        break;
      case OPT_O: // Set output directory with no clobber
        options->output_dirname = optarg;
        options->allow_clobber = false;
        break;
      case OPT_OC: // Set output directory with clobber
        options->output_dirname = optarg;
        options->allow_clobber = true;
        break;
      case OPT_SEQLEN:
        options->seq_len = atoi(optarg);
        break;
      case OPT_SCORE:
        options->score_thresh = atof(optarg);
        break;
      case OPT_MOTIF:
        rbtree_put(options->selected_motifs, optarg, NULL);
        break;
      case OPT_MOTIF_PSEUDO:
        options->pseudocount = atof(optarg);
        break;
      case OPT_ETHRESH:
        options->evalue_thresh = atof(optarg);
        break;
      case OPT_MINREG:
        options->min_win = atoi(optarg);
        break;
      case OPT_MAXREG:
        options->max_win = atoi(optarg);
        break;
      case OPT_NORC:
        options->scan_both_strands = false;
        break;
      case OPT_SEP:
        options->scan_separately = true;
        break;
      case OPT_FLIP:
        options->flip = true;
        break;
      case OPT_NOFLIP: // deprecated & intentionally undocumented
        options->flip = false;
        break;
      case OPT_DESC:
        options->description = optarg;
        break;
      case OPT_DFILE:
        options->desc_file = optarg;
        break;
      case OPT_LOCAL:
        options->local = true;
        break;
      case OPT_NOSEQ:
        options->noseq = true;
        break;
      case OPT_NEG:
        options->negseq_source = optarg;
        options->neg_sequences = true;
        break;
      case OPT_DISC:
        options->disc = true;
        options->mcc = true;
        break;
      case OPT_OPTIMIZE_SCORE:
        options->optimize_score = true;
        break;
      case OPT_VERBOSITY:
        verbosity = atoi(optarg);
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case '?':
        usage(NULL);
        break;
      default: // just in case we forget to handle a option
        die("Unhandled option %d", opt);
    }
  }
  // Must have sequence and motif file names
  if (argc < optind + 2) usage("Sequences and motifs are both required");
  if (options->scan_separately) options->scan_both_strands = false;
  if (options->disc && !options->neg_sequences)
    usage("You must supply a FASTA file with the negative dataset (--neg <fasta file>)"
        " when you use discriminative mode (--disc).");
  if (options->mcc && !options->neg_sequences)
    usage("You must supply a FASTA file with the negative dataset (--neg <fasta file>)"
        " to calculate Matthews correlation coefficient (--mcc).");
  if (options->max_win && options->min_win >= options->max_win)
    usage("You must specify --minreg smaller than --maxreg.");
  // Record the input file names
  options->seq_source = argv[optind++];
  for (; optind < argc; optind++)
    arraylst_add(argv[optind], options->motif_sources);

}

/*************************************************************************
 * score_motif_site
 *
 * Calculate the log odds score for a single motif-sized window.
 * Store scores that are better than what we've seen before.
 *************************************************************************/
static inline void score_motif_site(const char strand, const int start,
    int8_t *seq, PSSM_T *pssm, SEQ_SITES_T* seq_sites) {
  int motif_position;
  SEQ_SITE_T *site;
  double score = 0.0;
  const MATRIX_T* pssm_matrix = pssm->matrix;
  const int w = pssm->w;
  // For each position in the site
  for (motif_position = 0; motif_position < w; motif_position++, seq++) {
    // Check for gaps and ambiguity codes at this site
    if (*seq == -1) return;
    score += get_matrix_cell(motif_position, *seq, pssm_matrix);
  }
  // don't bother recording worse scores
  if (score < seq_sites->best) return;
  // better scores clear the list
  if (score > seq_sites->best) {
    seq_sites->used = 0;
    seq_sites->best = score;
  }
  // allocate memory on demand
  if (seq_sites->allocated <= seq_sites->used) {
    seq_sites->allocated += SCORE_BLOCK;
    mm_resize(seq_sites->sites, seq_sites->allocated, SEQ_SITE_T);
  }
  // store the site
  site = seq_sites->sites + (seq_sites->used++);
  site->start = start;
  site->strand = strand;
}

/*************************************************************************
 * score_sequence
 *
 * Calculate the log-odds score for each possible motif site in the
 * sequence and record the sites of the best. Apply a count to each
 * best site and increment the total site count.
 *************************************************************************/
static void score_sequence(CENTRIMO_OPTIONS_T *options, SEQ_T* sequence, 
    PSSM_T* pssm, PSSM_T* rev_pssm, SEQ_SITES_T* seq_sites, int seq_number,
    BEST_SITES_T* best) {
  int8_t *iseq;
  int i, pos;
  int *a_site;
  double count;
  SEQ_SITE_T *site;
  BEST_SITE_T *a_score;
  // get some constants
  const int L = get_seq_length(sequence);
  const int w = pssm->w;
  const int end = L - w + 1;
  // check we got passed stuff
  assert(options != NULL);
  assert(sequence != NULL);
  assert(pssm != NULL);
  assert(seq_sites != NULL);
  // Reset the sequence stats structure
  seq_sites->best = -BIG;
  seq_sites->used = 0;
  // Read and score each position in the sequence.
  iseq = get_isequence(sequence);
  if (rev_pssm) {
    // score both strands
    for (i = 0; i < end; i++, iseq++) {
      score_motif_site('+', i, iseq, pssm, seq_sites);
      score_motif_site('-', i, iseq, rev_pssm, seq_sites);
    }
  } else {
    // score forward strand
    for (i = 0; i < end; i++, iseq++) 
      score_motif_site('+', i, iseq, pssm, seq_sites);
  }
  // test that sites were found
  if (seq_sites->used == 0) return;
  // I assume that unscaling a score will get the same result for the
  // reverse complement as otherwise I can't do this optimisation.
  assert(rev_pssm == NULL || 
      (get_pssm_scale(rev_pssm) == get_pssm_scale(pssm) && 
       get_pssm_offset(rev_pssm) == get_pssm_offset(pssm)));
  // Handle scores that are out of range
  if ((int)seq_sites->best >= get_array_length(pssm->pv)) {
    seq_sites->best = (double)(get_array_length(pssm->pv) - 1);
  }
  // convert score to unscaled
  seq_sites->best = get_unscaled_pssm_score(seq_sites->best, pssm);
  // test that it passes the threshold
  if (seq_sites->best < options->score_thresh) return;
  // Record the position of best site
  if (best->used_scores >= best->allocated_scores) {
    best->allocated_scores *= 2;
    best->scores = mm_realloc(best->scores, sizeof(BEST_SITE_T) * best->allocated_scores);
  }
  a_score = best->scores+(best->used_scores);
  a_score->score = seq_sites->best;
  a_score->sequence_number = seq_number;
  a_score->offset = best->used_sites;
  a_score->used = seq_sites->used;
  best->used_scores++;
  // allocate memory for the sites, should be very rare
  if ((best->used_sites + seq_sites->used) > best->allocated_sites) {
    do {
      best->allocated_sites *= 2;
    } while ((best->used_sites + seq_sites->used) >= best->allocated_sites);
    best->sites = mm_realloc(best->sites, sizeof(int) * best->allocated_sites);
  }
  for (i = 0; i < seq_sites->used; i++) {
    site = seq_sites->sites + i;
    if (!options->flip || site->strand == '+') {
      //pos = 2 * (site->start + w/2 - 1/2); // a motif of width 1 can have sites at the first index
      pos = 2 * site->start + w - 1; // a motif of width 1 can have sites at the first index
    } else {
      //pos = 2 * (L - (site->start + w/2) - 1; // a motif of width 1 can have sites at the first index
      pos = 2 * (L - site->start) - w - 1;
    }
    best->sites[best->used_sites++] = pos;
  }
}

/*************************************************************************
 * score_sequences
 *
 * Score all the sequences and store the best location(s) in each sequence
 * that score better than the threshold in a list. Sort and return the list
 * so the best scores are at the start.
 **************************************************************************/
static void score_sequences(CENTRIMO_OPTIONS_T *options,
    BEST_SITES_T *best_sites, SEQ_SITES_T *sites, SEQ_T** sequences, int seqN, 
    PSSM_T* pos_pssm, PSSM_T* rev_pssm) {
  int i;
  best_sites->used_sites = 0;
  best_sites->used_scores = 0;
  for (i = 0; i < seqN; i++) {
    score_sequence(options, sequences[i], 
        pos_pssm, rev_pssm, sites, i, best_sites);
  }
  qsort(best_sites->scores, best_sites->used_scores, sizeof(BEST_SITE_T), compare_best_site);
}

/*************************************************************************
 * window_binomial
 *
 * Compute the binomial enrichment of a window
 *************************************************************************/
static inline double window_binomial(double window_sites, long total_sites,
    int bins, int max_bins) {
  // p-value vars
  double log_p_value, n_trials, n_successes, p_success;
  // calculate the log p-value
  if (window_sites == 0 || bins == max_bins) {
    log_p_value = 0; // pvalue of 1
  } else {
    n_trials = (double) total_sites;
    n_successes = window_sites;
    p_success = (double) bins / (double) max_bins;
    log_p_value = log_betai(n_successes, n_trials - n_successes + 1.0, p_success);
  }
  return log_p_value;
}

/*************************************************************************
 * window_FET
 *
 * Compute the fisher's exact test on a window
 *************************************************************************/
static inline double window_FET(double window_sites, long total_sites, 
    double neg_window_sites, long neg_total_sites) {
  // I'm not sure the rounds below are required, and if they are,
  // shouldn't they be in the function itself as it takes double type?
  // They were in Tom's implementation so for now I'm keeping them.
  return getLogFETPvalue(
      ROUND(window_sites),
      total_sites, 
      ROUND(neg_window_sites), 
      neg_total_sites, 
      1 // skips expensive calculation for p-values > 0.5, just returns 0.
    );
}

/*************************************************************************
 * window_MCC
 *
 * Compute the Matthews Correlation Coefficient of a window
 *************************************************************************/
static inline double window_MCC(double pos_window_sites, double neg_window_sites,
    int seqN, int neg_seqN) {
  double TP = (int) pos_window_sites;  // round down
  double FP = (int) neg_window_sites;  // round down
  double TN = neg_seqN-FP;
  double FN = seqN - TP;
  double base_sqr;
  if ((base_sqr = (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) == 0) {
    return 0;
  } else {
    return ((TP * TN) - (FP * FN)) / sqrt(base_sqr);
  }
}


/*************************************************************************
 * Accumulate a counts struct
 *************************************************************************/
static inline void accumulate_site_counts(SITE_COUNTS_T* counts, 
    SITE_COUNTS_T* target) {
  int i;
  assert(counts->allocated == target->allocated);
  target->total_sites = counts->total_sites;
  target->sites[0] = counts->sites[0];
  for (i = 1; i < counts->allocated; i++) {
    target->sites[i] = target->sites[i-1] + counts->sites[i];
  }
}
/*************************************************************************
 * Reset a counts struct
 *************************************************************************/
static inline void reset_site_counts(SITE_COUNTS_T* counts) {
  int i;
  for (i = 0; i < counts->allocated; i++) counts->sites[i] = 0;
  counts->total_sites = 0;
}
/*************************************************************************
 * Create and zero a counts struct
 *************************************************************************/
static inline SITE_COUNTS_T* create_site_counts(int seqlen) {
  int i;
  SITE_COUNTS_T* counts;
  counts = mm_malloc(sizeof(SITE_COUNTS_T));
  counts->allocated = ((2 * seqlen) - 1);
  counts->sites = mm_malloc(sizeof(double) * counts->allocated);
  reset_site_counts(counts);
  return counts;
}

/*************************************************************************
 * Destroy the site counts struct
 *************************************************************************/
static inline void destroy_site_counts(SITE_COUNTS_T* counts) {
  if (counts == NULL) return;
  free(counts->sites);
  memset(counts, 0, sizeof(SITE_COUNTS_T));
  free(counts);
}

/*************************************************************************
 * remove_redundant
 *
 * Given a list of windows sorted best to worst, remove any that overlap
 * with a better one.
 *************************************************************************/
static void remove_redundant(ARRAYLST_T* windows) {
  int i, j, end, extent, removed;
  int best_left, best_right, current_left, current_right;
  WIN_STATS_T *current, *best;
  // a single stat can't be redundant
  if (arraylst_size(windows) <= 1) return;
  removed = 0;
  best = NULL;
  for (i = 1; i < arraylst_size(windows); i++) { // For each stats (tested stats)
    current = (WIN_STATS_T*) arraylst_get(i, windows);
    extent = (current->spread - 1) / 2;
    current_left = current->center - extent;
    current_right = current->center + extent;
    // For each stats better than the tested one (better stats)
    end = i - removed;
    for (j = 0; j < end; j++) {
      best = (WIN_STATS_T*) arraylst_get(j, windows);
      extent = (best->spread - 1) / 2;
      best_left = best->center - extent;
      best_right = best->center + extent;
      if (best_right >= current_left && current_right >= best_left) break;
    }
    if (j < end) { // the window overlaps with a better one
      DEBUG_FMT(HIGH_VERBOSE, "Discarding redundant region "
          "(center: %d spread: %d) because it overlaps with a better "
          "region (center: %d spread: %d).\n", current->center,
          current->spread, best->center, best->spread);
      destroy_window(current);
      removed++;
    } else if (removed) {
      arraylst_set(i-removed, arraylst_get(i, windows), windows);
    }
  }
  if (removed) {
    arraylst_remove_range(arraylst_size(windows) - removed, removed, NULL, windows);
  }
}

/*************************************************************************
 * Output motif site counts
 *************************************************************************/
static void output_site_counts(FILE* fh, bool scan_separately,
    int sequence_length, MOTIF_DB_T* db, MOTIF_T* motif, SITE_COUNTS_T* counts,
    SITE_COUNTS_T* neg_counts) {
  // vars
  int i, w, end;
  char *alt;
  fprintf(fh, "DB %d MOTIF\t%s", db->id, 
      (scan_separately ? get_motif_st_id(motif) : get_motif_id(motif)));
  alt = get_motif_id2(motif);
  if (alt[0]) fprintf(fh, "\t%s", alt);
  fprintf(fh, "\n");
  w = get_motif_length(motif);
  end = counts->allocated - (w - 1);
  for (i = (w - 1); i < end; i += 2) {
    if (neg_counts) {
      fprintf(fh, "% 6.1f\t%g\t%g\n", ((double) (i - sequence_length + 1)) / 2.0, counts->sites[i], neg_counts->sites[i]);
    } else {
      fprintf(fh, "% 6.1f\t%g\n", ((double) (i - sequence_length + 1)) / 2.0, counts->sites[i]);
    }
  }
}

/*************************************************************************
 * Setup the JSON writer and output a lot of pre-calculation data
 *************************************************************************/
static void start_json(CENTRIMO_OPTIONS_T* options, int argc, char** argv, 
    ARRAY_T* bg_freqs, SEQ_T** sequences, int seqN, int seq_skipped, int neg_seqN, 
    int neg_seq_skipped, MOTIF_DB_T** dbs, int motifN, int seqlen,
    HTMLWR_T** html_out, JSONWR_T** json_out) {
  int i;
  MOTIF_DB_T* db;
  HTMLWR_T *html;
  JSONWR_T *json;
  // setup html monolith writer
  json = NULL;
  if ((html = htmlwr_create(get_meme_etc_dir(), TEMPLATE_FILENAME, false))) {
    htmlwr_set_dest_name(html, options->output_dirname, HTML_FILENAME);
    htmlwr_replace(html, "centrimo_data.js", "data");
    json = htmlwr_output(html);
    if (json == NULL) die("Template does not contain data section.\n");
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "Failed to open html template file.\n");
    *html_out = NULL;
    *json_out = NULL;
    return;
  }
  // now output some json
  // output some top level variables
  jsonwr_str_prop(json, "version", VERSION);
  jsonwr_str_prop(json, "revision", REVISION);
  jsonwr_str_prop(json, "release", ARCHIVE_DATE);
  //jsonwr_str_prop(json, "program", options->local ? "LocoMo" : "CentriMo");
  jsonwr_str_prop(json, "program", "CentriMo");
  // output cmd, have
  jsonwr_args_prop(json, "cmd", argc, argv);
  jsonwr_property(json, "options");
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "bfile", options->bg_source);
  jsonwr_dbl_prop(json, "motif-pseudo", options->pseudocount);
  jsonwr_dbl_prop(json, "score", options->score_thresh);
  jsonwr_bool_prop(json, "optimize_score", options->optimize_score);
  jsonwr_dbl_prop(json, "ethresh", options->evalue_thresh);
  jsonwr_lng_prop(json, "minbin", options->min_win);
  jsonwr_lng_prop(json, "maxbin", options->max_win);
  jsonwr_bool_prop(json, "local", options->local);
  jsonwr_bool_prop(json, "norc", !options->scan_both_strands);
  jsonwr_bool_prop(json, "sep", !options->scan_separately);
  jsonwr_bool_prop(json, "flip", options->flip);
  jsonwr_bool_prop(json, "noseq", options->noseq);
  jsonwr_bool_prop(json, "neg_sequences", options->neg_sequences);
  jsonwr_bool_prop(json, "disc", options->disc);
  jsonwr_bool_prop(json, "mcc", options->mcc);
  jsonwr_end_object_value(json);
  // output description
  jsonwr_desc_prop(json, "job_description", options->desc_file, options->description);
  // output size metrics
  jsonwr_lng_prop(json, "seqlen", seqlen);
  jsonwr_lng_prop(json, "tested", motifN);
  // output alphabet
  ALPH_T *alph = options->alphabet;
  jsonwr_property(json, "alphabet");
  alph_print_json(alph, json);
  // output background
  jsonwr_property(json, "background");
  jsonwr_start_array_value(json);
  for (i = 0; i < alph_size_core(alph); i++) {
    jsonwr_dbl_value(json, get_array_item(i, bg_freqs));
  }
  jsonwr_end_array_value(json);
  // output the fasta db
  jsonwr_property(json, "sequence_db");
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "source", options->seq_source);
  jsonwr_lng_prop(json, "count", seqN);
  jsonwr_lng_prop(json, "skipped", seq_skipped);
  jsonwr_end_object_value(json);
  if (options->neg_sequences) {
    jsonwr_property(json, "negative_sequence_db");
    jsonwr_start_object_value(json);
    jsonwr_str_prop(json, "source", options->negseq_source);
    jsonwr_lng_prop(json, "count", neg_seqN);
    jsonwr_lng_prop(json, "skipped", neg_seq_skipped);
    jsonwr_end_object_value(json);
  }
  // output the motif dbs
  jsonwr_property(json, "motif_dbs");
  jsonwr_start_array_value(json);
  for (i = 0; i < arraylst_size(options->motif_sources); i++) {
    db = dbs[i];
    jsonwr_start_object_value(json);
    jsonwr_str_prop(json, "source", db->source);
    jsonwr_lng_prop(json, "count", arraylst_size(db->motifs));
    jsonwr_end_object_value(json);
  }
  jsonwr_end_array_value(json);
  if (!options->noseq) {
    // output the sequences ID
    jsonwr_property(json, "sequences");
    jsonwr_start_array_value(json);
      for (i = 0; i < seqN; i++) {
        jsonwr_str_value(json, get_seq_name(sequences[i]));
      }
    jsonwr_end_array_value(json);
  }
  // start the motif array
  jsonwr_property(json, "motifs");
  jsonwr_start_array_value(json);
  // return the html and json writers
  *html_out = html;
  *json_out = json;
}

/*************************************************************************
 * Output JSON data for a motif.
 *************************************************************************/
static void output_motif_json(JSONWR_T* json, int sequence_length,
    MOTIF_STATS_T *stats, SITE_COUNTS_T* counts, SITE_COUNTS_T* neg_counts,
    bool store_sequences, bool negative_sequences,
    bool discriminative, bool mcc, bool scan_separately) {
  WIN_STATS_T *window;
  MOTIF_T *motif;
  MATRIX_T *freqs;
  int i, j, mlen, asize, end, index;
  RBNODE_T *seq;
  motif = stats->motif;
  freqs = get_motif_freqs(motif);
  asize = alph_size_core(get_motif_alph(motif));
  jsonwr_start_object_value(json);
  jsonwr_lng_prop(json, "db", stats->db->id);
  jsonwr_str_prop(json, "id", scan_separately ? get_motif_st_id(motif) : get_motif_id(motif));
  if (*(get_motif_id2(motif))) jsonwr_str_prop(json, "alt", get_motif_id2(motif));
  if (get_motif_consensus(motif)) jsonwr_str_prop(json, "consensus", get_motif_consensus(motif));
  mlen = get_motif_length(motif);
  jsonwr_lng_prop(json, "len", mlen);
  jsonwr_log10num_prop(json, "motif_evalue", get_motif_log_evalue(motif), 1);
  if (get_motif_nsites(motif) > 0) jsonwr_dbl_prop(json, "motif_nsites", get_motif_nsites(motif));
  jsonwr_lng_prop(json, "n_tested", stats->n_tests);
  //jsonwr_dbl_prop(json, "max_prob", stats->max_prob);
  jsonwr_dbl_prop(json, "score_threshold", stats->score_threshold);
  if (get_motif_url(motif) && *get_motif_url(motif)) jsonwr_str_prop(json, "url", get_motif_url(motif));
  jsonwr_property(json, "pwm");
  jsonwr_start_array_value(json);
  for (i = 0; i < mlen; i++) {
    jsonwr_start_array_value(json);
    for (j = 0; j < asize; j++) {
      jsonwr_dbl_value(json, get_matrix_cell(i, j, freqs));
    }
    jsonwr_end_array_value(json);
  }
  jsonwr_end_array_value(json);
  jsonwr_lng_prop(json, "total_sites", counts->total_sites);
  jsonwr_property(json, "sites");
  jsonwr_start_array_value(json);
  end = counts->allocated - (mlen - 1);
  for (i = (mlen - 1); i < end; i += 2) {
    jsonwr_dbl_value(json, counts->sites[i]);
  }
  jsonwr_end_array_value(json);
  if (negative_sequences) {
    jsonwr_lng_prop(json, "neg_total_sites", neg_counts->total_sites);
    jsonwr_property(json, "neg_sites");
    jsonwr_start_array_value(json);
    end = neg_counts->allocated - (mlen - 1);
    for (i = (mlen - 1); i < end; i += 2) {
      jsonwr_dbl_value(json, neg_counts->sites[i]);
    }
    jsonwr_end_array_value(json);
  }
  jsonwr_property(json, "seqs");
  jsonwr_start_array_value(json);
  if (store_sequences) {
    for (seq = rbtree_first(stats->seq_ids); seq; seq = rbtree_next(seq)) {
      jsonwr_lng_value(json, (long)(*((int*)rbtree_key(seq))));
    }
  }
  jsonwr_end_array_value(json);
  jsonwr_property(json, "peaks"); // There are several possible peaks for LocoMo output
  jsonwr_start_array_value(json);
  for (index = 0; index < arraylst_size(stats->windows); index++) {
    jsonwr_start_object_value(json);
    window = (WIN_STATS_T *) arraylst_get(index, stats->windows);
    jsonwr_dbl_prop(json, "center", ((double) (window->center - sequence_length + 1)) / 2.0);
    jsonwr_lng_prop(json, "spread", (window->spread + 1) / 2);
    jsonwr_dbl_prop(json, "sites", window->sites);
    jsonwr_dbl_prop(json, "log_adj_pvalue", window->log_adj_pvalue);
    if (negative_sequences) {
      jsonwr_dbl_prop(json, "neg_sites", window->neg_sites);
      jsonwr_dbl_prop(json, "neg_log_adj_pvalue", window->neg_log_adj_pvalue);
      // Compute Fisher's Exact Test for peak
      window->fisher_log_pvalue = window_FET(window->sites, counts->total_sites, window->neg_sites,
        neg_counts->total_sites); 
      window->fisher_log_adj_pvalue = LOGEV(log(stats->n_tests), window->fisher_log_pvalue);
      jsonwr_dbl_prop(json, "fisher_log_adj_pvalue", window->fisher_log_adj_pvalue);
      if (mcc) jsonwr_dbl_prop(json, "mcc", window->mcc);
    }
    jsonwr_end_object_value(json);
  }
  jsonwr_end_array_value(json);
  jsonwr_end_object_value(json);
}

/*************************************************************************
 * Finish up JSON output and HTML output
 *************************************************************************/
static void end_json(HTMLWR_T* html, JSONWR_T* json) {
  // finish writing motifs
  if (json) jsonwr_end_array_value(json);
  // finish writing html file
  if (html) {
    if (htmlwr_output(html) != NULL) {
      die("Found another JSON replacement!\n");
    }
    htmlwr_destroy(html);
  }
}

/*************************************************************************
 * Output a log value in scientific notation
 *************************************************************************/
static void print_log_value(FILE *file, double loge_val, int prec) {
  double log10_val, e, m;
  log10_val = loge_val / log(10);
  e = floor(log10_val);
  m = pow(10.0, (log10_val - e));
  if ((m + (0.5 * pow(10, -prec))) >= 10) {
    m = 1;
    e += 1;
  }
  fprintf(file, "%.*fe%04.0f", prec, m, e);
}

/*************************************************************************
 * CentriMo sites (create file)
 *************************************************************************/
static FILE* start_centrimo_sites(CENTRIMO_OPTIONS_T *options) {
  char *sites_path;
  FILE *sites_file;
  sites_path = make_path_to_file(options->output_dirname, SITES_FILENAME);
  sites_file = fopen(sites_path, "w");
  free(sites_path);
  return sites_file;
}

/*************************************************************************
 * CentriMo text (create file)
 *************************************************************************/
static FILE* start_centrimo_text(CENTRIMO_OPTIONS_T *options) {
  char *file_path;
  FILE *text_file;
  // open centrimo text file
  file_path = make_path_to_file(options->output_dirname, TEXT_FILENAME);
  text_file = fopen(file_path, "w");
  free(file_path);
  fputs("# WARNING: this file is not sorted!\n"
      "# db\tid                         alt\t"
      "consensus\t"
      " E-value\tadj_p-value\tlog_adj_p-value\t"
      "bin_location\tbin_width\ttotal_width\tsites_in_bin\ttotal_sites\t"
      "p_success\t p-value\tmult_tests", text_file);
  if (options->neg_sequences) {
    fputs("\tneg_sites_in_bin\tneg_sites\tneg_adj_pvalue\tlog_neg_adj_pvalue\tfisher_adj_pvalue\tlog_fisher_adj_pvalue", text_file);
  }
  fputs("\n", text_file);
  return text_file;
}

/*************************************************************************
 * CentriMo text (output motif stats)
 *************************************************************************/
static void output_centrimo_text(FILE* text_file, bool scan_separately,
   bool negative_sequences, MOTIF_STATS_T* stats, int motifN, int sequence_length) {
  //vars
  MOTIF_DB_T *db;
  MOTIF_T *motif;
  WIN_STATS_T *window;
  double log_motifN, actual_center, window_prob;
  int pad;
  char* (*get_id)(MOTIF_T*);

  get_id = (scan_separately ? get_motif_st_id : get_motif_id);
  log_motifN = log(motifN);
  motif = stats->motif;
  window = arraylst_get(0, stats->windows);
  db = stats->db;
  // convert from the internal coordinate system
  // adjust so the center of the sequence is reported as 0
  actual_center = ((double) (window->center - sequence_length + 1)) / 2.0;
  // calculate the probability of landing in the window
  window_prob = (double) window->n_bins / stats->n_bins;

  // write the motif DB
  fprintf(text_file, "%4d\t", db->id + 1);
  // write the motif name and alternate name, pad it so things nicely align
  pad = 30 - strlen(get_id(motif)) - 1;
  if (pad < 0) pad = 0;
  char *alt = get_motif_id2(motif);
  fprintf(text_file, "%s %*s\t%s\t", get_id(motif), pad, alt[0] ? alt : get_motif_id(motif),
    get_motif_consensus(motif));
  // write the best window E-value
  print_log_value(text_file, window->log_adj_pvalue + log_motifN, 1);
  fputs("\t", text_file);
  // write the best window adjusted p-value
  fputs("   ", text_file);
  print_log_value(text_file, window->log_adj_pvalue, 1);
  fputs("\t", text_file);
  // write the best window log adjusted p-value
  fprintf(text_file, "%15.2f\t", window->log_adj_pvalue);
  // write the best window stats as follows:
  // center, included bins, total bins, included sites, total sites, probability
  // need to use round function on window sites because fprintf does bankers rounding
  // which is different to the javascript.
  fprintf(text_file, "%12.1f\t%9d\t%11d\t%12.0f\t%11ld\t%9.5f\t",
      actual_center, window->n_bins, stats->n_bins, 
      round(window->sites), stats->sites, window_prob);
  // write the best window log raw p-value
  print_log_value(text_file, window->log_pvalue, 1);
  fputs("\t", text_file);
  // write the number of windows tested
  fprintf(text_file, "%10d", stats->n_tests);
  if (negative_sequences) {
    fprintf(text_file, "\t%12.0f\t%11ld\t", window->neg_sites, stats->neg_sites);
    // write adjusted p-value and its log
    print_log_value(text_file, window->neg_log_adj_pvalue, 1);
    fputs("\t", text_file);
    fprintf(text_file, "%15.2f\t\t", window->neg_log_adj_pvalue);
    // write adjusted Fisher's p-value and its log
    print_log_value(text_file, window->fisher_log_adj_pvalue, 1);
    fputs("\t", text_file);
    fprintf(text_file, "%15.2f\t", window->fisher_log_adj_pvalue);
  }
  fprintf(text_file, "\n");
}

/*************************************************************************
 * Read all the sequences into an array of SEQ_T
 *************************************************************************/
static void read_sequences(ALPH_T *alph, char *seq_file_name,
    SEQ_T ***sequences, int *seq_num, int *seq_skipped, int *seq_len) {
  const int max_sequence = 32768; // unlikely to be this big
  int i, move;
  FILE * seq_fh = fopen(seq_file_name, "r");
  if (!seq_fh) die("failed to open sequence file `%s'", seq_file_name);
  *seq_num = 0;
  read_many_fastas(alph, seq_fh, max_sequence, seq_num, sequences);
  if (fclose(seq_fh) != 0) die("failed to close sequence file\n");
  if (*seq_len == 0) *seq_len = get_seq_length((*sequences)[0]);
  // remove sequences that don't match the seq_len value
  move = 0;
  for (i = 0; i < *seq_num; i++) {
    if (*seq_len == get_seq_length((*sequences)[i])) {
      // convert sequence into indexed form (don't keep information about ambiguous characters)
      index_sequence((*sequences)[i], alph, SEQ_NOAMBIG);
      // move to fill gaps in array
      if (move > 0) (*sequences)[i - move] = (*sequences)[i];
    } else {
      fprintf(stderr, "Skipping sequence %s as its length (%d) does not "
          "match the expected length (%d).\n", get_seq_name((*sequences)[i]), 
          get_seq_length((*sequences)[i]), *seq_len);
      free_seq((*sequences)[i]);
      move++;
    }
  }
  *seq_num -= move;
  *seq_skipped = move;
  for (i--; i >= *seq_num; i--) (*sequences)[i] = NULL;
  *sequences = mm_realloc(*sequences, sizeof(SEQ_T*) * (*seq_num));
  if (*seq_num == 0) die("Failed to find a single sequence with the"
      " expected length %d\n", *seq_len);  
}

/*************************************************************************
 * Read a motif database
 *************************************************************************/
static MOTIF_DB_T* read_motifs(CENTRIMO_OPTIONS_T *options, int id,
    char* motif_source, ARRAY_T** bg) 
{
  // vars
  int read_motifs, i;
  MOTIF_DB_T* motifdb;
  MREAD_T *mread;
  MOTIF_T *motif;
  ARRAYLST_T *motifs;
  RBTREE_T *seen;
  ALPH_T *db_alph;

  // Load the motifs and the background.
  mread = mread_create(motif_source, OPEN_MFILE);
  mread_set_pseudocount(mread, options->pseudocount);
  if (id == 0) { 			// Get background first time only
    if (options->bg_source) {           // background from file or motifs or uniform
      mread_set_bg_source(mread, options->bg_source);
      *bg = mread_get_background(mread);
      options->bg_source = mread_get_other_bg_src(mread);
    } else {				// background from sequences
      *bg = calc_bg_from_file(options->alphabet, options->seq_source, false);
      mread_set_background(mread, *bg);
      options->bg_source = strdup("--sequences--");
    }
  } else {				// Use previously read background.
    mread_set_background(mread, *bg);
  }

  if (!alph_equal(options->alphabet, mread_get_alphabet(mread))) {
    mread_set_conversion(mread, options->alphabet, *bg);
  }

  // load motifs
  read_motifs = 0;
  if (rbtree_size(options->selected_motifs) > 0) {
    motifs = arraylst_create();
    while (mread_has_motif(mread)) {
      motif = mread_next_motif(mread);
      read_motifs++;
      if (rbtree_find(options->selected_motifs, get_motif_id(motif))) {
        arraylst_add(motif, motifs);
      } else {
        DEBUG_FMT(NORMAL_VERBOSE, "Discarding motif %s in %s (not selected).\n",
            get_motif_id(motif), motif_source);
        destroy_motif(motif);
      }
    }
  } else {
    motifs = mread_load(mread, NULL);
    read_motifs = arraylst_size(motifs);    
  }
  // remove duplicate ID motifs
  seen = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, NULL);
  for (i = 0; i < arraylst_size(motifs); i++) {
    motif = (MOTIF_T*)arraylst_get(i, motifs);
    if (!rbtree_make(seen, get_motif_id(motif), NULL)) {
      DEBUG_FMT(NORMAL_VERBOSE, "Discarding motif %s in %s (non-unique ID).\n",
          get_motif_id(motif), motif_source);
      arraylst_remove(i--, motifs);
      destroy_motif(motif);
    }
  }
  rbtree_destroy(seen);

  // Add reverse complements if scanning strands separately
  if (options->scan_separately) add_reverse_complements(motifs);

  arraylst_fit(motifs);
  if (read_motifs > 0) {
    // check the alphabet
    if (!alph_equal(options->alphabet, get_motif_alph((MOTIF_T*)arraylst_peek(motifs)))) {
      die("Expected %s alphabet motifs\n", alph_name(options->alphabet));
    }
    // get the background
    if (*bg == NULL) *bg = mread_get_background(mread);
  } else {
    fprintf(stderr, "Warning: Motif file %s contains no motifs.\n", motif_source);
  }
  // clean up motif reader
  mread_destroy(mread);
  // create motif db
  motifdb = mm_malloc(sizeof(MOTIF_DB_T));
  memset(motifdb, 0, sizeof(MOTIF_DB_T));
  motifdb->id = id;
  motifdb->source = strdup(motif_source);
  motifdb->motifs = motifs;
  return motifdb;
}

/*************************************************************************
 * Free a motif database
 *************************************************************************/
static void free_db(MOTIF_DB_T* db) {
  free(db->source);
  arraylst_destroy(destroy_motif, db->motifs);
  memset(db, 0, sizeof(MOTIF_DB_T));
  free(db);
}

/*************************************************************************
 * Make a PSSM for centrimo
 *************************************************************************/
static PSSM_T* make_pssm(ARRAY_T* bg_freqs, MOTIF_T *motif) {
  PSSM_T *pssm;

// Build PSSM for motif and tables for p-value calculation.
// p-values are not used in centrimo, only scores
  pssm = build_motif_pssm(motif, bg_freqs, bg_freqs, //p-value background
      NULL, // prior distribution
      ALPHA, //non-motif specific scale factor
      PSSM_RANGE, 0, // no GC bins
      false // make log-likelihood pssm
      );
  return pssm;
}

/*************************************************************************
 * test_counts_sites
 *
 * For a given counts struct, make sure that the sum of each site is equals to total sites value
 *************************************************************************/
static void test_counts_sites(SITE_COUNTS_T *counts) {
  int i;
  double sum_check = 0, sum_diff;
  for (i = 0; i < counts->allocated; i++)
    sum_check += counts->sites[i];
  sum_diff = counts->total_sites - sum_check;
  if (sum_diff < 0)
    sum_diff = -sum_diff;
  if (sum_diff > 0.1) {
    fprintf(stderr, "Warning: site counts don't sum to accurate value! %g != %ld\n", sum_check, counts->total_sites);
  }
}

/*************************************************************************
 * compute_counts
 *
 *************************************************************************/
static void compute_counts(BEST_SITES_T* best, SITE_COUNTS_T *counts, int* start_index, double score_thresh) {
  int i, j;
  double count;
  BEST_SITE_T *score;
  for (i = (start_index ? *start_index : 0); i < best->used_scores; i++) { // For each score
    score = best->scores+i;
    if (score->score >= score_thresh) {
      count = 1.0 / score->used; // Compute count value to add at every position
      for (j = 0; j < score->used; j++) { // For each position of the current score
        counts->sites[best->sites[score->offset + j]] += count; // Set values into counts
      }
      counts->total_sites++; // Increase total sites
    } else {
      break;
    }
  }
  if (start_index) *start_index = i;
  test_counts_sites(counts);
}

/*************************************************************************
 * all_sequences_in_window
 *
 * Given a window, an arraylst of scores and a score threshold; this function
 * will find all sequences where the best motif site is in the window.
 *
 * Assumes the list of scores is sorted best to worst.
 *************************************************************************/
static RBTREE_T* all_sequences_in_window(int window_center, int window_spread,
    double min_score_threshold, BEST_SITES_T *best_sites) {
  RBTREE_T *seq_ids;
  int window_extent, window_left, window_right, i , j, pos;
  BEST_SITE_T *score;
  seq_ids = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);
  // Get the minimum and maximum locations of the window
  window_extent = (window_spread - 1) / 2;
  window_left = window_center - window_extent;
  window_right = window_center + window_extent;
  // now see if any positions in any scores fit in the window
  for (i = 0; i < best_sites->used_scores; i++) {
    score = best_sites->scores+i;
    if (score->score < min_score_threshold) {
      // as this list is sorted there won't be any better scores
      break;
    }
    for (j = 0; j < score->used; j++) { 
      pos = best_sites->sites[score->offset + j];
      if (pos >= window_left && pos <= window_right) { 
        // the position is inside the given window
        rbtree_make(seq_ids, &(score->sequence_number), NULL);
      }
    }
  }
  return seq_ids;
}

/*************************************************************************
 * smallest_score
 *
 * Helper function to get the smallest score from the positive and negative
 * scores. Assumes scores are sorted largest to smallest.
 *************************************************************************/
static inline double smallest_score(BEST_SITES_T* pos, BEST_SITES_T* neg) {
  int len;
  double value, value2;
  value = BIG;
  if (pos && (len = pos->used_scores) > 0) {
    value = pos->scores[len - 1].score;
  }
  if (neg && (len = neg->used_scores) > 0) {
    value2 = neg->scores[len - 1].score;
    value = MIN(value, value2);
  }
  return value;
}

/*************************************************************************
 * largest_score
 *
 * Helper function to get the largest score at an index from the 
 * positive and negative scores. Assumes scores are sorted largest to smallest.
 *************************************************************************/
static inline double largest_score(BEST_SITES_T* pos, int pos_i, BEST_SITES_T* neg, int neg_i) {
  double value, value2;
  value = -BIG;
  if (pos && pos_i < pos->used_scores) {
    value = pos->scores[pos_i].score;
  }
  if (neg && neg_i < neg->used_scores) {
    value2 = neg->scores[neg_i].score;
    value = MAX(value, value2);
  }
  return value;
}


/*************************************************************************
 * count_unique_scores
 *
 * Counts the number of unique scores. This information is then used to do
 * a multiple test correction.
 *************************************************************************/
static int count_unique_scores(BEST_SITES_T* pos, BEST_SITES_T* neg) {
  double last;
  int p_i, n_i, count;
  p_i = n_i = 0;
  count = 0;
  last = largest_score(pos, p_i, neg, n_i);
  while (last != -BIG) {
    count++;
    // skip over scores with the same value
    if (pos) {
      for (; p_i < pos->used_scores; p_i++) {
        if (pos->scores[p_i].score != last) break;
      }
    }
    if (neg) {
      for (; n_i < neg->used_scores; n_i++) {
        if (neg->scores[n_i].score != last) break;
      }
    }
    // get the next largest score
    last = largest_score(pos, p_i, neg, n_i);
  }
  return count;
}

/*************************************************************************
 * test_window
 *
 * Perform statistical tests on a window.
 *************************************************************************/
static void test_window(CENTRIMO_OPTIONS_T* options, double log_pvalue_thresh,
    int n_bins, double log_n_tests,
    SITE_COUNTS_T* pve_c_counts, SITE_COUNTS_T* neg_c_counts, 
    int center, int spread, 
    double* best_ignore, ARRAYLST_T* windows, int seqN, int neg_seqN) {
  WIN_STATS_T *window;
  int expand, left, right, bins;
  double pve_sites, neg_sites, log_pvalue, log_adj_pvalue;
  neg_sites = 0;
  expand = (spread - 1) / 2; // spread is guaranteed to be +ve and odd
  left = center - expand;
  right = center + expand;
  // Warning: assumes that a window with empty outermost bins will not be passed
  bins = expand + 1;

  // calculate the number of sites using the cumulative counts
  pve_sites = pve_c_counts->sites[right] - (left ? pve_c_counts->sites[left - 1] : 0);
  if (options->neg_sequences) {
    neg_sites = neg_c_counts->sites[right] - (left ? neg_c_counts->sites[left - 1] : 0);
  }
  // calculate window significance
  if (options->disc) {
    // use Fisher's Exact Test
    log_pvalue = window_FET(pve_sites, pve_c_counts->total_sites, neg_sites,
        neg_c_counts->total_sites); 
  } else {
    // use Binomial Test
    if (pve_sites <= *best_ignore) {
      // we've previously seen a smaller or equal window that had at least this
      // many sites and we skipped it because it did not meet the threshold
      // it follows that we can skip this one too because it won't do better
      return;
    }
    log_pvalue = window_binomial(pve_sites, pve_c_counts->total_sites, 
        bins, n_bins);
  }
  log_adj_pvalue = LOGEV(log_n_tests, log_pvalue);
  if (log_adj_pvalue > log_pvalue_thresh) {
    // Note: best_ignore only used as a shortcut for binomial testing
    *best_ignore = pve_sites;
    return;
  }
  // this looks like an interesting window! Record some statistics...
  window = mm_malloc(sizeof(WIN_STATS_T));
  memset(window, 0, sizeof(WIN_STATS_T));
  window->center = center;
  window->spread = spread;
  window->n_bins = bins;
  window->sites = pve_sites;
  window->log_pvalue = log_pvalue;
  window->log_adj_pvalue = log_adj_pvalue;
  if (options->neg_sequences) {
    window->neg_sites = neg_sites;
    window->neg_log_pvalue = window_binomial(neg_sites, 
      neg_c_counts->total_sites, bins, n_bins);
    window->neg_log_adj_pvalue = LOGEV(log_n_tests, window->neg_log_pvalue);
    if (options->mcc) {
      window->mcc = window_MCC(pve_sites, neg_sites, seqN, neg_seqN);
    }
  }
  arraylst_add(window, windows);
} // test_window

/*************************************************************************
 * calculate_best_windows
 *
 *************************************************************************/
static MOTIF_STATS_T* calculate_best_windows(CENTRIMO_OPTIONS_T* options, 
    CENTRIMO_BUFFER_T* buffers, double log_pvalue_thresh, 
    MOTIF_DB_T *db, MOTIF_T* motif, int seq_len, 
    BEST_SITES_T* pve_scores, BEST_SITES_T* neg_scores,
    int seqN, int neg_seqN) {
  int n_bins, n_minus, min_win, max_win, n_windows, n_scores, n_tests, n_actual_tests;
  int spread, max_spread, pos_first, pos_last, pos;
  int i, pve_index, neg_index;
  long best_total_sites, best_neg_total_sites;
  double score_thresh_min, best_score_thresh, log_n_tests;
  ARRAYLST_T *best_windows, *current_windows, *temp;
  WIN_STATS_T *win_stats;
  MOTIF_STATS_T *motif_stats;
  double best_log_pv, current_log_pv, best_ignore_count;
  // initilise vars
  best_total_sites = 0;
  best_neg_total_sites = 0;
  best_score_thresh = BIG;
  if (pve_scores->used_scores == 0) {
    // no scores in the positive set, hence no best windows
    return NULL;
  }
  //
  // calculate the multiple test correction
  //
  n_bins = seq_len - get_motif_length(motif) + 1;
  // Insure that we have at least one "negative bin", max_win < n_bins
  max_win = n_bins - 1;
  if (options->max_win && options->max_win >= n_bins) {
    fprintf(stderr, "--maxreg (%d) too large for motif; setting to %d\n", options->max_win, max_win);
  } else if (options->max_win) {
    max_win = options->max_win;
  }
  // Insure that min_win is less than max_win so that even/odd parity situations both work.
  min_win = 1;
  if (options->min_win && options->min_win >= max_win) {
    min_win = max_win - 1;
    fprintf(stderr, "--minreg (%d) too large; setting to %d\n", options->min_win, min_win);
  } else if (options->min_win) {
    min_win = options->min_win;
  }
  // If the sequence_length:motif_width parity is EVEN, then
  // the allowable window sizes are ODD: 1, 3, 5, ...
  // If the sequence_length:motif_width parity is ODD, then
  // the allowable window sizes are EVEN: 2, 4, 6, ...
  // So we adjust min_win UP and max_win DOWN so that they have correct ODD/EVENness.
  int even = (n_bins % 2);
  if (even == ((min_win+1) % 2)) min_win++;
  if (even == ((max_win+1) % 2)) max_win--;
  n_scores = (options->optimize_score ? count_unique_scores(pve_scores, neg_scores) : 1);
  if (options->local) {
    // A_w = n_bins - w + 1
    // n_windows = \sum_w={min_win}^{max_win} A_w
    //  = \sum_w={1}^{max_win} A_w + \sum_w={min_win-1}^{max_win} A_w
    // algebra gives the following:
    n_windows = (n_bins+1)*(max_win-min_win+1) - (max_win*(max_win+1) - min_win*(min_win-1))/2;
  } else {
    //n_windows = ((max_win - (n_bins % 2)) / 2) + (n_bins % 2);
    // adjust max_win, min_win down if seq_length:motif_width parity is even
    // add in central bin if parity is even and min_win == 1
    n_windows = ((max_win - even) / 2) - ((min_win - 1 - even) / 2) + ((even && min_win==1) ? 1 : 0);
  }
  n_tests = n_windows * n_scores;
  log_n_tests = log(n_tests);
  // fprintf(stderr, "w %d l %d min_win %d max_win %d n_windows %d n_scores %d n_tests %d\n", get_motif_length(motif), seq_len, min_win, max_win, n_windows, n_scores, n_tests);
  // reset vars
  pve_index = neg_index = 0;
  reset_site_counts(buffers->pos_counts);
  if (options->neg_sequences) reset_site_counts(buffers->neg_counts);
  best_windows = arraylst_create();
  current_windows = arraylst_create();
  best_log_pv = BIG;
  n_actual_tests = 0;
  max_spread = 2 * max_win - 1;
  // try each score threshold in turn
  while (1) {
    best_ignore_count = -1;
    // get the scoring threshold
    if (options->optimize_score) {
      score_thresh_min = largest_score(pve_scores, pve_index, neg_scores, neg_index);
      if (score_thresh_min == -BIG) break;
    } else {
      score_thresh_min = smallest_score(pve_scores, neg_scores);
    }
    // compute counts and advance the index to the first score worse than the threshold
    compute_counts(pve_scores, buffers->pos_counts, &pve_index, score_thresh_min);
    //convert to cumulative site counts
    accumulate_site_counts(buffers->pos_counts, buffers->pos_c_counts);
    // do same for negative set
    if (options->neg_sequences) {
      compute_counts(neg_scores, buffers->neg_counts, &neg_index, score_thresh_min);
      accumulate_site_counts(buffers->neg_counts, buffers->neg_c_counts);
    }

    if (options->local) {
      // calculate the first and last position for a window size of 1
      // pos_first = get_motif_length(motif) - 1;
      // pos_last = 2 * seq_len - get_motif_length(motif) - 1;
      // calculate the first and last position for a window size of min_win 
      pos_first = get_motif_length(motif) - 1 + (min_win-1);
      pos_last = 2 * seq_len - get_motif_length(motif) - 1 - (min_win-1);
      // loop over window sizes increasing
      //for (spread = 1; spread <= max_spread; spread += 2, pos_first++, pos_last--) {
      //fprintf(stderr, "pos_first %d pos_last %d max_spread %d n_bins %d min_win %d max_win %d\n", pos_first, pos_last, max_spread, n_bins, min_win, max_win);
      for (spread = 2*min_win-1; spread <= max_spread; spread += 2, pos_first++, pos_last--) {
        //loop over center positions
        for (pos = pos_first; pos <= pos_last; pos += 2) {
          test_window(options, log_pvalue_thresh, n_bins, log_n_tests,
              buffers->pos_c_counts, buffers->neg_c_counts, pos, spread,
              &best_ignore_count, current_windows, seqN, neg_seqN);
          n_actual_tests++;
        }
      }
    } else {
      pos = seq_len - 1;
      // loop over centered window sizes increasing
      for (spread=2*min_win-1; spread <= max_spread; spread += 4) {
        test_window(options, log_pvalue_thresh, n_bins, log_n_tests,
            buffers->pos_c_counts, buffers->neg_c_counts, pos, spread,
            &best_ignore_count, current_windows, seqN, neg_seqN);
        n_actual_tests++;
      }
    }
    if (!options->optimize_score) {
      // swap
      temp = best_windows;best_windows = current_windows;current_windows = temp;
      // keep track of best windows stats
      best_score_thresh = score_thresh_min;
      best_total_sites = buffers->pos_counts->total_sites;
      if (options->neg_sequences) {
        best_neg_total_sites = buffers->neg_counts->total_sites;
      }
      break;
    }

    // find the best score for the current threshold
    current_log_pv = BIG;
    for (i = 0; i < arraylst_size(current_windows); i++) {
      win_stats = (WIN_STATS_T*) arraylst_get(i, current_windows);
      if (win_stats->log_adj_pvalue < current_log_pv) {
        current_log_pv = win_stats->log_adj_pvalue;
      }
    }
    // see if we found something better!
    if (current_log_pv < best_log_pv) {
      // swap
      temp = best_windows;best_windows = current_windows;current_windows = temp;
      // update pvalue
      best_log_pv = current_log_pv;
      // keep track of best windows stats
      best_score_thresh = score_thresh_min;
      best_total_sites = buffers->pos_counts->total_sites;
      if (options->neg_sequences) {
        best_neg_total_sites = buffers->neg_counts->total_sites;
      }
    }
    arraylst_clear(destroy_window, current_windows);
  }
  // fprintf(stderr, "n_tests %d n_actual %d\n", n_tests, n_actual_tests);
  assert(n_tests == n_actual_tests);
  arraylst_destroy(destroy_window, current_windows);
  arraylst_qsort(window_stats_compare_pvalue, best_windows);
  // remove any weaker overlapping windows
  remove_redundant(best_windows);
  arraylst_fit(best_windows);
  // check if we found at least one good window
  if (arraylst_size(best_windows) == 0) {
    // nothing useful found
    arraylst_destroy(destroy_window, best_windows);
    return NULL;
  }
  // found some good windows, make motif stats
  motif_stats = mm_malloc(sizeof(MOTIF_STATS_T));
  memset(motif_stats, 0, sizeof(MOTIF_STATS_T));
  motif_stats->db = db;
  motif_stats->motif = motif;
  motif_stats->windows = best_windows;
  if (!options->noseq) { 
    win_stats = arraylst_get(0, best_windows);
    motif_stats->seq_ids = all_sequences_in_window(win_stats->center,
        win_stats->spread, best_score_thresh, pve_scores);
  }
  motif_stats->n_tests = n_tests;
  motif_stats->n_bins = n_bins;
  motif_stats->score_threshold = options->optimize_score ? best_score_thresh : options->score_thresh;
  motif_stats->sites = best_total_sites;
  motif_stats->neg_sites = best_neg_total_sites;
  return motif_stats;
} // calculate_best_windows


/*************************************************************************
 * Allocates memory for a best sites buffer
 *************************************************************************/
BEST_SITES_T* create_best_sites(int seqlen, int seqn) {
  BEST_SITES_T *bsites;
  bsites = mm_malloc(sizeof(BEST_SITES_T));
  memset(bsites, 0, sizeof(BEST_SITES_T));
  bsites->allocated_scores = 2;
  bsites->scores = mm_malloc(sizeof(BEST_SITE_T) * bsites->allocated_scores);
  bsites->allocated_sites = 2;
  bsites->sites = mm_malloc(sizeof(int) * bsites->allocated_sites);
  return bsites;
}

/*************************************************************************
 * Frees memory for a best sites buffer
 *************************************************************************/
void destroy_best_sites(BEST_SITES_T *bsites) {
  if (bsites == NULL) return;
  if (bsites->scores) free(bsites->scores);
  if (bsites->sites) free(bsites->sites);
  free(bsites);
}

/*************************************************************************
 * Allocates memory for some buffers that centrimo needs
 *************************************************************************/
void create_buffers(CENTRIMO_OPTIONS_T* options, CENTRIMO_BUFFER_T* buffers, int seqlen, int pos_seqN, int neg_seqN) {
  memset(buffers, 0, sizeof(CENTRIMO_BUFFER_T));
  // create sites
  buffers->sites = mm_malloc(sizeof(SEQ_SITES_T));
  memset(buffers->sites, 0, sizeof(SEQ_SITES_T));
  // create best sites
  buffers->pos_sites = create_best_sites(seqlen, pos_seqN);
  if (options->neg_sequences) buffers->neg_sites = create_best_sites(seqlen, neg_seqN);
  // create counts
  buffers->pos_counts = create_site_counts(seqlen);
  if (options->neg_sequences) buffers->neg_counts = create_site_counts(seqlen);
  buffers->pos_c_counts = create_site_counts(seqlen);
  if (options->neg_sequences) buffers->neg_c_counts = create_site_counts(seqlen);
}

/*************************************************************************
 * Frees memory for centrimo's buffers
 *************************************************************************/
void destroy_buffers(CENTRIMO_BUFFER_T* buffers) {
  destroy_site_counts(buffers->pos_counts);
  destroy_site_counts(buffers->neg_counts);
  destroy_site_counts(buffers->pos_c_counts);
  destroy_site_counts(buffers->neg_c_counts);
  destroy_best_sites(buffers->pos_sites);
  destroy_best_sites(buffers->neg_sites);
  free(buffers->sites->sites);
  free(buffers->sites);
  memset(buffers, 0, sizeof(CENTRIMO_BUFFER_T)); 
}

/*************************************************************************
 * Writes documentation at end of text output file.
 *************************************************************************/
void end_centrimo_text(FILE *text_file, bool neg_sequences) {
  fputs("##\n# Detailed descriptions of columns in this file:\n#\n"
      "# db:\tThe name of the database (file name) that contains the motif.\n"
      "# id:\tA name for the motif that is unique in the motif database file.\n"
      "# alt:\tAn alternate name of the motif that may be provided\n"
      "#\tin the motif database file.\n"
      "# consensus:\tA consensus sequence computed from the motif.\n"
      "# E-value:\tThe expected number motifs that would have least one.\n"
      "#\tregion as enriched for best matches to the motif as the reported region.\n"
      "#\tThe E-value is the p-value multiplied by the number of motifs in the\n"
      "#\tinput database(s).\n"
      "# adj_p-value:\tThe probability that any tested region would be as enriched for\n"
      "#\tbest matches to this motif as the reported region is.\n"
      "#\tBy default the p-value is calculated by using the one-tailed binomial\n"
      "#\ttest on the number of sequences with a match to the motif \n"
      "#\tthat have their best match in the reported region, corrected for\n"
      "#\tthe number of regions and score thresholds tested.\n"
      "#\tThe test assumes that the probability that the best match in a sequence\n"
      "#\tfalls in the region is the region width divided by the\n"
      "#\tnumber of places a motif\n"
      "#\tcan align in the sequence (sequence length minus motif width plus 1).\n"
      "#\tWhen CentriMo is run in discriminative mode with a negative\n"
      "#\tset of sequences, the p-value of a region is calculated\n"
      "#\tusing the Fisher exact test on the \n"
      "#\tenrichment of best matches in the positive sequences relative\n"
      "#\tto the negative sequences, corrected\n"
      "#\tfor the number of regions and score thresholds tested.\n"
      "#\tThe test assumes that the probability that the best match (if any)\n"
      "#\tfalls into a given region\n"
      "#\tis the same for all positive and negative sequences.\n"
      "# log_adj_p-value:\tLog of adjusted p-value.\n"
      "# bin_location:\tLocation of the center of the most enriched region.\n"
      "# bin_width:\tThe width (in sequence positions) of the most enriched region.\n"
      "#\tA best match to the motif is counted as being in the region if the\n"
      "#\tcenter of the motif falls in the region.\n"
      "# total_width:\tThe window maximal size which can be reached for this motif:\n"
      "#\t\trounded(sequence length - motif length +1)/2\n"
      "# sites_in_bin:\tThe number of (positive) sequences whose best match to the motif\n"
      "#\tfalls in the reported region.\n"
      "#\tNote: This number may be less than the number of\n"
      "#\t(positive) sequences that have a best match in the region.\n"
      "#\tThe reason for this is that a sequence may have many matches that score\n"
      "#\tequally best.\n"
      "#\tIf n matches have the best score in a sequence, 1/n is added to the\n"
      "#\tappropriate bin for each match.\n"
      "# total_sites:\tThe number of sequences containing a match to the motif\n"
      "#\tabove the score threshold.\n"
      "# p_success:\tThe probability of falling in the enriched window:\n"
      "#\t\tbin width / total width\n"
      "# p-value:\tThe uncorrected p-value before it gets adjusted to the\n"
      "#\tnumber of multiple tests to give the adjusted p-value.\n"
      "# mult_tests:\tThis is the number of multiple tests (n) done for this motif.\n"
      "#\tIt was used to correct the original p-value of a region for\n"
      "#\tmultiple tests using the formula:\n"
      "#\t\tp' = 1 - (1-p)^n where p is the uncorrected p-value.\n"
      "#\tThe number of multiple tests is the number of regions\n"
      "#\tconsidered times the number of score thresholds considered.\n"
      "#\tIt depends on the motif length, sequence length, and the type of\n"
      "#\toptimizations being done (central enrichment, local enrichment,\n"
      "#\tscore optimization).\n",
      text_file);
  if (neg_sequences) {
    fputs("# neg_sites_in_bin:\tThe number of negative sequences where the best\n"
        "#\tmatch to the motif falls in the reported region.\n"
        "#\tThis value is rounded but the underlying value may contain\n"
        "#\tfractional counts.\n"
        "#\tNote: This number may be less than the number of negative have a\n"
        "#\tbest match in the region.\n"
        "#\tThe reason for this is that a sequence may have many matches that\n"
        "#\tscore equally best.\n"
        "#\tIf n matches have the best score in a sequence, 1/n is added to the\n"
        "#\tappropriate bin for each match.\n"
        "# neg_sites:\tThe number of negative sequences containing a match to the\n"
        "#\tmotif above the minimum score threshold.\n"
        "#\tWhen score optimization is enabled the score threshold may be raised\n"
        "#\thigher than the minimum.\n"
        "# neg_adj_pvalue:\tThe probability that any tested region in the negative\n"
        "#\tsequences would be as enriched for best matches to this motif\n"
        "#\taccording to the Binomial test.\n"
        "# log_neg_adj_pvalue:\tLog of negative adjusted p-value.\n"
        "# fisher_adj_pvalue:\tFisher adjusted p-value before it gets adjusted to the\n"
        "#\tnumber of motifs in the input database(s).\n"
        "#\tRefers to the E-value definition using the discriminative mode.\n"
        "# log_fisher_adj_pvalue:\tLog of Fisher adjusted p-value.\n",text_file);
  }
}  // end_centrimo_text

/*************************************************************************
 * Entry point for centrimo
 *************************************************************************/
int main(int argc, char *argv[]) {
  CENTRIMO_OPTIONS_T options;
  CENTRIMO_BUFFER_T buffers;
  int seqlen, seqN, seq_skipped, neg_seqN, neg_seq_skipped, motifN;
  int i, db_i, motif_i;
  double log_pvalue_thresh, log_half;
  SEQ_T **sequences, **neg_sequences;
  ARRAY_T* bg_freqs;
  MOTIF_DB_T **dbs, *db;
  MOTIF_STATS_T *motif_stats;
  MOTIF_T *motif, *rev_motif;
  PSSM_T *pos_pssm, *rev_pssm;
  char *sites_path, *desc;
  FILE *sites_file, *text_out;
  HTMLWR_T *html;
  JSONWR_T *json;
  
  // command line processing
  process_command_line(argc, argv, &options);

  options.alphabet = NULL;
  if (options.alph_file != NULL) {
    options.alphabet = alph_load(options.alph_file, true);
    if (options.alphabet == NULL) exit(EXIT_FAILURE);
  }

  // Read all the alphabets and make sure they are the same.
  read_motif_alphabets(options.motif_sources, options.alph_file != NULL, &(options.alphabet));

  if (!alph_has_complement(options.alphabet)) {
    if (options.scan_separately)
      usage("The option --sep can not be used with an unstranded alphabet.");
    options.scan_both_strands = false;
  }

  // calculate the background (if not provided)
  bg_freqs = NULL;

  // load the positive sequences
  DEBUG_MSG(NORMAL_VERBOSE, "Loading sequences.\n");
  seqlen = options.seq_len;
  read_sequences(options.alphabet, options.seq_source, &sequences, &seqN, 
      &seq_skipped, &seqlen);

  // load the negative sequences
  if (options.neg_sequences) {
    read_sequences(options.alphabet, options.negseq_source, &neg_sequences,
        &neg_seqN, &neg_seq_skipped, &seqlen);
  } else {
    neg_seqN = 0;
    neg_seq_skipped = 0;
  }

  // load the motifs
  DEBUG_MSG(NORMAL_VERBOSE, "Loading motifs.\n");
  dbs = mm_malloc(sizeof(MOTIF_DB_T*) * arraylst_size(options.motif_sources));
  for (motifN = 0, i = 0; i < arraylst_size(options.motif_sources); i++) {
    char* db_source;
    db_source = (char*) arraylst_get(i, options.motif_sources);
    dbs[i] = read_motifs(&options, i, db_source, &bg_freqs);
    motifN += arraylst_size(dbs[i]->motifs);
  }
  // convert the evalue threshold into a pvalue threshold
  log_pvalue_thresh = log(options.evalue_thresh) - log(motifN);
  // if p-value threshold would be 1.0, reduce it slightly to
  // prevent jillions of absolutely non-significant peaks being printed
  //if (log_pvalue_thresh >= 0) log_pvalue_thresh = log(0.999999999);
  // The above made things really confusing if E=1 and there was one
  // motif with a p-value=1.  Nothing got printed.  So I changed it
  // to allow the pvalue_thresh to be 1.0.
  if (log_pvalue_thresh > 0) log_pvalue_thresh = 0;

  // Setup some things for double strand scanning
  if (options.scan_both_strands == true) {
    // Correct background by averaging on freq. for both strands.
    average_freq_with_complement(options.alphabet, bg_freqs);
    normalize_subarray(0, alph_size_core(options.alphabet), 0.0, bg_freqs);
  }

  // Create output directory
  if (create_output_directory(options.output_dirname, options.allow_clobber,
        (verbosity >= NORMAL_VERBOSE))) {
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }

  // open output files
  sites_file = start_centrimo_sites(&options);
  text_out = start_centrimo_text(&options);
  start_json(&options, argc, argv, bg_freqs, sequences, seqN, seq_skipped, neg_seqN, 
      neg_seq_skipped, dbs, motifN, seqlen, &html, &json);

  // initialize local variables
  create_buffers(&options, &buffers, seqlen, seqN, neg_seqN);
  motif = rev_motif = NULL;
  pos_pssm = rev_pssm = NULL;

  // calculate and output the best windows for each motif
  for (db_i = 0, i = 1; db_i < arraylst_size(options.motif_sources); db_i++) {
    db = dbs[db_i];
    for (motif_i = 0; motif_i < arraylst_size(db->motifs); motif_i++, i++) {
      motif = (MOTIF_T *) arraylst_get(motif_i, db->motifs);
      DEBUG_FMT(NORMAL_VERBOSE, "Using motif %s (%d/%d) of width %d.\n", 
          options.scan_separately ? get_motif_st_id(motif) : get_motif_id(motif), i, motifN, get_motif_length(motif));
      // create the pssm
      pos_pssm = make_pssm(bg_freqs, motif);
      // If required, do the same for the reverse complement motif.
      if (options.scan_both_strands) {
        rev_motif = dup_rc_motif(motif);
        rev_pssm = make_pssm(bg_freqs, rev_motif);
      }
      // score the sequences
      score_sequences(&options, buffers.pos_sites, buffers.sites, sequences, seqN, pos_pssm, rev_pssm);
      if (options.neg_sequences) {
        score_sequences(&options, buffers.neg_sites, buffers.sites, neg_sequences, neg_seqN, pos_pssm, rev_pssm);
      }
      // calculate the best windows 
      motif_stats = calculate_best_windows(&options, &buffers,
          log_pvalue_thresh, db, motif, seqlen, buffers.pos_sites, buffers.neg_sites, seqN, neg_seqN);
      // If there is a result to ouput
      if (motif_stats) {
        // recalculate the site counts for this specific score
        reset_site_counts(buffers.pos_counts);
        compute_counts(buffers.pos_sites, buffers.pos_counts, NULL, motif_stats->score_threshold);
        if (options.neg_sequences) {
          reset_site_counts(buffers.neg_counts);
          compute_counts(buffers.neg_sites, buffers.neg_counts, NULL, motif_stats->score_threshold);
        }
        // Output JSON results
        if (json) {
          output_motif_json(json, seqlen, motif_stats, buffers.pos_counts, 
              buffers.neg_counts, !options.noseq, options.neg_sequences, 
              options.disc, options.mcc, options.scan_separately); // Write values into HTML file
        }
        output_site_counts(sites_file, options.scan_separately, seqlen, db, motif, buffers.pos_counts, buffers.neg_counts);
        // Note: Fisher's stuff is computed in output_motif_json so it better be called!
        output_centrimo_text(text_out, options.scan_separately, options.neg_sequences,
          motif_stats, motifN, seqlen);
      }
      // Free memory associated with this motif.
      destroy_stats(motif_stats);
      free_pssm(pos_pssm);
      free_pssm(rev_pssm);
      destroy_motif(rev_motif);
    } // looping over motifs
  } // looping over motif databases

  end_json(html, json);
  end_centrimo_text(text_out, options.neg_sequences);
  // finish writing files
  fclose(sites_file);
  fclose(text_out);
  // Clean up.
  destroy_buffers(&buffers);
  for (i = 0; i < seqN; ++i) {
    free_seq(sequences[i]);
  }
  free(sequences);
  if (options.neg_sequences) {
    for (i = 0; i < neg_seqN; ++i) {
      free_seq(neg_sequences[i]);
    }
    free(neg_sequences);
  }
  for (i = 0; i < arraylst_size(options.motif_sources); i++) {
    free_db(dbs[i]);
  }
  free(dbs);
  free_array(bg_freqs);
  cleanup_options(&options);
  DEBUG_MSG(NORMAL_VERBOSE, "Program ends correctly\n");
  return EXIT_SUCCESS;
} // main
