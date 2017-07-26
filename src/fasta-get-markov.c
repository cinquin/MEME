#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/select.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "alphabet.h"
#include "string-builder.h"

#define BAD -1
#define NEWLINE -2
#define WHITESPACE -3

// GLOBAL VARIABLE
volatile bool update_status = false;

// STRUCTURES
typedef struct fmkv_options {
  char* alph; // alphabet file, if NULL then use atype
  char atype; // alphabet type 'n' = DNA, 'p' = Protein, 'f' = Full list of seen symbols, 'u' = Unknown (guess DNA or Protein based on counts)
  char* fasta;
  char* bg;
  char* seqc; // sequence count output file
  int order;
  double pseudo;
  bool norc;
  bool counts; // output sparse counts instead of frequencies
  bool status;
  bool machine;
  bool summary;
} FMKV_OPTIONS_T;

typedef struct fmkv_count FMKV_COUNT_T;

struct fmkv_count {
  int64_t count;
  int size;
  FMKV_COUNT_T **counts;
};

typedef struct fmkv_stats {
  off_t fasta_size;
  off_t fasta_proc;
  int64_t total_seqs;
  int64_t total_slen;
  int64_t slen_min;
  int64_t slen_max;
  bool bad;
  off_t bad_pos;
  STR_T *symbols; // all the symbols we've seen
  off_t *symbols_offset;
  FMKV_COUNT_T counts; // the counts of all the symbols we've seen
} FMKV_STATS_T;

typedef struct fmkv_state FMKV_STATE_T;

struct fmkv_state {
  // tracks if the previous letter was a newline
  bool nl;
  // tracks the current sequence length
  int64_t slen;
  // tracks the current method of processing the file
  size_t (*routine)(FMKV_OPTIONS_T*, FMKV_STATS_T*, FMKV_STATE_T*, char*, size_t);
  // does quick mapping of a symbol to an index. All symbols out of range
  // are mapped to index 1. Any unassigned symbols have the index 0 until they
  // are seen and assigned an index.
  int8_t hasher[UCHAR_MAX + 1];
  // number of symbols seen since last reset, 
  int64_t seen;
  // pointers into the hierarchical letter counts
  // history[0] == &(stats->counts)
  // history[1] is the pointer to history[0]->counts[x] where x is the index
  //    assigned to the last seen symbol
  // To update the history we start at the end and work back to the start.
  FMKV_COUNT_T **history;
};

static void cleanup_options(FMKV_OPTIONS_T *options) {
}

static void cleanup_counts(FMKV_COUNT_T *entry) {
  FMKV_COUNT_T *count;
  int i;
  for (i = 0; i < entry->size; i++) {
    if (entry->counts[i] != NULL) {
      cleanup_counts(entry->counts[i]);
      free(entry->counts[i]);
      entry->counts[i] = NULL;
    }
  }
  free(entry->counts);
  memset(entry, 0, sizeof(FMKV_COUNT_T));
}

static void cleanup_stats(FMKV_OPTIONS_T *options, FMKV_STATS_T *stats) {
  cleanup_counts(&(stats->counts));
}

static void usage(FMKV_OPTIONS_T *options, char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n"
    "    fasta-get-markov [options] [sequence file] [background file]\n"
    "Options:\n"
    "    [-m <order>]        order of Markov model to use; default 0\n"
    "    [-alph <alphabet file>] use the specified custom alphabet;\n"
    "                        default: guess alphabet\n"
    "    [-dna]              use DNA alphabet; default: guess alphabet\n"
    "    [-protein]          use protein alphabet; default: guess alphabet\n"
    "    [-rna]              use rna alphabet; default: guess alphabet\n"
    "    [-full]             use full list of seen symbols as the alphabet;\n"
    "                        default: guess alphabet\n"
    "    [-counts]           instead of a traditional markov model output counts\n"
    "                        and skip entries with no counts; implies \"-norc\";\n"
    "                        default: output frequencies\n"
    "    [-norc]             do not combine forward and reverse complement freqs;\n"
    "                        this is highly recommended for RNA sequences.\n"
    "    [-pseudo <count>]   pseudocount added to avoid probabilities of zero;\n"
    "                        default: use a pseudocount of 0.1 .\n"
    "    [-nostatus]         do not print status messages.\n"
    "    [-nosummary]        do not print the summary report even when a\n"
    "                        background file is specified.\n"
    "    [-help]             display this help message.\n"
    "Description:\n"
    "    Estimate a Markov model from a FASTA file of sequences.\n"
    "    Skips tuples containing ambiguous symbols.\n"
    "    Combines both strands of complementable alphabets unless -norc is set.\n\n"
    "    Reads standard input if a sequence file is not specified.\n"
    "    Writes standard output if a background file is not specified.\n\n"
    "    When the background file is specified the following report is made:\n"
    "    <sequence count> <min length> <max length> <average length> <summed length>\n\n";

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
  cleanup_options(options);
  if (format) exit(EXIT_FAILURE);
  exit(EXIT_SUCCESS);
}

#define MRKV 0
#define PROT 1
#define NORC 2
#define PSDO 3
#define NSTA 4
#define NSUM 5
#define HELP 6
#define SEQC 7
#define MECH 8
#define ALPH 9
#define DNAB 10
#define FULL 11
#define COUN 12
#define RNAB 13
static void process_arguments(FMKV_OPTIONS_T *options, int argc, char **argv) {
  struct option fmkv_options[] = {
    {"m", required_argument, NULL, MRKV}, // level of markov model
    {"seqc", required_argument, NULL, SEQC}, // sequence count, undocumented
    {"alph", required_argument, NULL, ALPH}, // Custom alphabet
    {"dna", no_argument, NULL, DNAB}, // DNA alphabet
    {"n", no_argument, NULL, DNAB}, // DNA alphabet
    {"protein", no_argument, NULL, PROT}, // Protein alphabet
    {"p", no_argument, NULL, PROT}, // Protein alphabet
    {"rna", no_argument, NULL, RNAB}, // RNA alphabet
    {"full", no_argument, NULL, FULL}, // Alphabet from full list of seen symbols
    {"counts", no_argument, NULL, COUN}, // output counts instead of frequencies
    {"norc", no_argument, NULL, NORC}, // average complementary chains
    {"pseudo", required_argument, NULL, PSDO}, // apply a pseudocount before calculation of frequencies
    {"nostatus", no_argument, NULL, NSTA}, // do not output status messages
    {"machine", no_argument, NULL, MECH}, // machine readable status, undocumented
    {"nosummary", no_argument, NULL, NSUM}, // do not display a summary
    {"help", no_argument, NULL, HELP}, // show the usages message
    {NULL, 0, NULL, 0} //boundary indicator
  };
  bool bad_argument = false;

  options->bg = "-";
  options->fasta = "-";
  options->alph = NULL;
  options->seqc = NULL;
  options->order = 0;
  options->pseudo = 0.1;
  options->atype = 'u'; // unknown
  options->counts = false;
  options->norc = false;
  options->status = true;
  options->machine = false;
  options->summary = true;

  while (1) {
    int opt = getopt_long_only(argc, argv, "", fmkv_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case MRKV: //-m <order>
        options->order = atoi(optarg);
        break;
      case SEQC: //-seqc <nseqs file>
        options->seqc = optarg;
        break;
      case ALPH: //-alph <alphabet file>
        options->alph = optarg;
        break;
      case DNAB: //-n
        options->atype = 'n';
        break;
      case RNAB: //-rna
        options->atype = 'r';
        break;
      case PROT: //-p
        options->atype = 'p';
        break;
      case FULL: //-f
        options->atype = 'f';
        break;
      case COUN: //-counts
        options->counts = true;
        break;
      case NORC: //-norc
        options->norc = true;
        break;
      case PSDO: //-pseudo <count>
        options->pseudo = strtod(optarg, NULL);
        break;
      case NSTA: //-nostatus
        options->status = false;
        break;
      case MECH: //-machine (machine readable status)
        options->machine = true;
        break;
      case NSUM: //-nosummary
        options->summary = false;
        break;
      case HELP: //-help
        usage(options, NULL); // does not return
      case '?':           //unrecognised or ambiguous argument
        bad_argument = true;
    }
  }
  if (bad_argument) usage(options, "One or more unknown or ambiguous options were supplied.");

  if (options->order < 0) usage(options, "Model order can not be less than zero.");
  if (options->pseudo < 0) usage(options, "Pseudo-count must be positive.");

  if (optind < argc) options->fasta = argv[optind++];
  if (optind < argc) options->bg = argv[optind++];
}

static void display_status_update(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats) {
  long double total, processed, fraction;
  if (opts->machine) {
    fprintf(stderr, "%jd\n", (intmax_t)stats->fasta_proc);
  } else if (stats->fasta_size) {
    total = stats->fasta_size;
    processed = stats->fasta_proc;
    fraction = (processed / total) * 100;
    fprintf(stderr, "\rprocessed: %.1Lf%%", fraction);
  } else {
    fprintf(stderr, "\rsequences: %" PRId64, stats->total_seqs);
  }
  update_status = false;
}

// predeclare routines so they can reference each-other.
static size_t routine_seq_data(FMKV_OPTIONS_T*, FMKV_STATS_T*, FMKV_STATE_T*, char*, size_t);
static size_t routine_seq_name(FMKV_OPTIONS_T*, FMKV_STATS_T*, FMKV_STATE_T*, char*, size_t);
static size_t routine_find_start(FMKV_OPTIONS_T*, FMKV_STATS_T*, FMKV_STATE_T*, char*, size_t);

// routine to read sequence data
static size_t routine_seq_data(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats,
    FMKV_STATE_T *state, char *buffer, size_t len) {
  size_t i, j, k, chain;
  int pos, order_plus, ascii_sym;
  FMKV_COUNT_T *now, *next;
  order_plus = opts->order + 1;
  if (state->nl) {
    state->nl = false;
    if (buffer[0] == '>') {
      state->routine = routine_seq_name;
      if (state->slen < stats->slen_min) stats->slen_min = state->slen;
      if (state->slen > stats->slen_max) stats->slen_max = state->slen;
      stats->total_slen += state->slen;
      return 1;
    }
  }
  for (i = 0; i < len; i++) {
    switch ((pos = state->hasher[(uint8_t)buffer[i]])) {
      case BAD: // illegal but just set flag and treat as a specal character
        if (!stats->bad) {
          // record the position of the first bad character sighted
          stats->bad_pos = stats->fasta_proc + i;
        }
        stats->bad = true;
        // chain interrupted -> erase history
        state->seen = 0;
        state->slen++;
        break;
      case NEWLINE: // newline character
        if ((i + 1) < len) {
          // don't bother using the flag as we can look ahead
          if (buffer[i+1] == '>') {
            state->nl = false;
            state->routine = routine_seq_name;
            if (state->slen < stats->slen_min) stats->slen_min = state->slen;
            if (state->slen > stats->slen_max) stats->slen_max = state->slen;
            stats->total_slen += state->slen;
            return i + 2;
          }
        } else { // last character in buffer so we need to use the flag
          state->nl = true;
        }
        // fall through
      case WHITESPACE: // whitespace charater
        break;
      case 0: // unassigned symbol
        ascii_sym = (buffer[i] & 0x7F);
        str_append_code(stats->symbols, ascii_sym);
        pos = str_len(stats->symbols);
        state->hasher[ascii_sym] = pos;
        stats->symbols_offset = mm_realloc(stats->symbols_offset, sizeof(off_t) * pos);
        stats->symbols_offset[pos - 1] = stats->fasta_proc + i;
        // fall through
      default:
        // process the symbol
        for (j = (state->seen < opts->order ? state->seen : opts->order) + 1; j > 0; j--) {
          now = state->history[j - 1];
          // extend the counts list to support the symbol
          if (now->size < pos) {
            int new_size;
            new_size = str_len(stats->symbols);
            now->counts = mm_realloc(now->counts, sizeof(FMKV_COUNT_T*) * new_size);
            for (k = now->size; k < new_size; k++) now->counts[k] = NULL;
            now->size = new_size;
          }
          // get the next symbol counter
          next = now->counts[pos - 1];
          // allocate a counter for the symbol if we didn't find one
          if (next == NULL) {
            next = mm_malloc(sizeof(FMKV_COUNT_T));
            next->count = 0;
            next->size = 0;
            next->counts = NULL;
            now->counts[pos - 1] = next;
          }
          next->count++;
          state->history[j] = next;
        }
        state->seen++;
        state->slen++;
        break;
    }
  }
  return len;
}

// routine to skip over the name/description line of the sequence
static size_t routine_seq_name(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats,
    FMKV_STATE_T *state, char *buffer, size_t len) {
  size_t i, j;
  for (i = 0; i < len; i++) {
    if (buffer[i] == '\n' || buffer[i] == '\r') {
      state->nl = true;
      state->routine = routine_seq_data;
      state->slen = 0;
      // reset history
      state->seen = 0;
      // increment sequence count
      stats->total_seqs++;
      return i + 1;
    }
  }
  return len;
}

// routine to find the first sequence
static size_t routine_find_start(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats,
    FMKV_STATE_T *state, char *buffer, size_t len) {
  size_t i;
  for (i = 0; i < len;) {
    if (state->nl && buffer[i] == '>') {
      state->routine = routine_seq_name;
      state->nl = false;
      return i + 1;
    } else {
      state->nl = false;
      for (; i < len; i++) {
        if (buffer[i] == '\n' || buffer[i] == '\r') {
          state->nl = true;
          i++;
          break;
        }
      }
    }
  }
  return len;
}

static void process_fasta2(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats, FILE *fasta_fh, char *buffer, size_t buffer_len) {
  int i;
  size_t bytes_read, offset, consumed;
  FMKV_STATE_T state;
  size_t (*prev_routine)(FMKV_OPTIONS_T*, FMKV_STATS_T*, FMKV_STATE_T*, char*, size_t);
  // init state
  state.nl = true;
  state.slen = 0;
  state.routine = routine_find_start;
  // setup the hasher
  for (i = 0; i <= UCHAR_MAX; i++) {
    // set everything to the bad index
    state.hasher[i] = BAD;
  }
  for (i = 33; i <= 126; i++) {
    // zero the visible ascii to allow them to be used
    state.hasher[i] = 0;
  }
  // add newlines
  state.hasher['\r'] = NEWLINE;
  state.hasher['\n'] = NEWLINE;
  // add whitespace
  state.hasher[' '] = WHITESPACE;
  state.hasher['\t'] = WHITESPACE;
  state.hasher['\v'] = WHITESPACE; //vertical tab, rather unlikely

  // setup the history
  state.seen = 0;
  state.history = mm_malloc(sizeof(FMKV_COUNT_T*) * (opts->order + 2));
  state.history[0] = &(stats->counts);
  for (i = 1; i <= (opts->order + 1); i++) state.history[i] = NULL;
  //load data into buffer
  while ((bytes_read = fread(buffer, sizeof(char), buffer_len, fasta_fh)) > 0) {
    offset = 0;
    while (offset < bytes_read) {
      prev_routine = state.routine;
      consumed = state.routine(opts, stats, &state, buffer+offset, bytes_read - offset);
      if (consumed == 0 && state.routine == prev_routine) {
        fprintf(stderr, "Infinite loop detected!\n");
        abort(); // we'd want to debug this so an abort makes more sense than exit.
      }
      offset += consumed;
      stats->fasta_proc += consumed;
      if (update_status) display_status_update(opts, stats);
    }
    if (feof(fasta_fh) || ferror(fasta_fh)) break;
  }
  if (state.routine == routine_seq_name) {
    stats->total_seqs++;
    stats->slen_min = 0;
  } else if (state.routine == routine_seq_data) {
    if (state.slen < stats->slen_min) stats->slen_min = state.slen;
    if (state.slen > stats->slen_max) stats->slen_max = state.slen;
    stats->total_slen += state.slen;
  } else if (state.routine == routine_find_start) {
    stats->slen_min = 0;
  }
  free(state.history);
}

static void process_fasta(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats) {
  FILE *fasta_fh;
  size_t buffer_len, i;
  char *buffer;
  struct stat fasta_stat;
  off_t fasta_size;
  // get file handle
  if (strcmp(opts->fasta, "-") == 0) {
    // As fasta-get-markov allows reading from stdin we can't just display
    // an error message when no input file is given and so people who have never
    // used it before may not know what is going on.
    // So when we're reading from a terminal we try to read for half a second
    // and if nothing happens we display the usage message.
    if (isatty(STDIN_FILENO)) {
      fd_set input_set;
      struct timeval timeout;
      FD_ZERO(&input_set);
      FD_SET(STDIN_FILENO, &input_set);
      timeout.tv_sec = 0;
      timeout.tv_usec = 500000; // half a second
      if (select(1, &input_set, NULL, NULL, &timeout) != 1) {
        usage(opts, "No input was received after half a second."); // does not return
      }
    }
    // read from stdin
    fasta_fh = stdin;
    fasta_size = 0;
  } else {
    // open file
    fasta_fh = fopen(opts->fasta, "r");
    // check it worked
    if (fasta_fh == NULL) {
      // report error and exit
      fprintf(stderr, "Failed to open fasta file \"%s\" for reading: ", opts->fasta);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
    // get the file size
    if (fstat(fileno(fasta_fh), &fasta_stat) == -1) {
      perror("Failed to get size of the fasta file");
      exit(EXIT_FAILURE);
    }
    fasta_size = fasta_stat.st_size;
  }
  // allocate 1MB buffer
  buffer_len =  (1 << 20);
  buffer = mm_malloc(sizeof(char) * buffer_len);
  // init stats
  stats->fasta_size = fasta_size;
  stats->fasta_proc = 0;
  stats->bad = false;
  stats->total_seqs = 0;
  stats->total_slen = 0;
  stats->slen_min = 0x7FFFFFFFFFFFFFFF;
  stats->slen_max = 0;
  stats->symbols = str_create(10);
  stats->symbols_offset = NULL;
  stats->counts.count = 0;
  stats->counts.size = 0;
  stats->counts.counts = NULL;
  // read and process sequences
  process_fasta2(opts, stats, fasta_fh, buffer, buffer_len);
  if (opts->status) {
    display_status_update(opts, stats);
    fputs("\r\033[K", stderr);
  }
  // clean up
  free(buffer);
  if (strcmp(opts->fasta, "-") != 0) {
    fclose(fasta_fh);
  }
}

static ALPH_T* determine_alphabet(FMKV_OPTIONS_T *opts, FMKV_STATS_T *stats) {
  if (opts->alph) {
    ALPH_T *alph;
    alph = alph_load(opts->alph, true);
    if (alph == NULL) exit(EXIT_FAILURE);
    return alph;
  } else {
    switch (opts->atype) {
      case 'n':
        return alph_dna();
      case 'r':
        return alph_rna();
      case 'p':
        return alph_protein();
      case 'f':
        return alph_generic(str_internal(stats->symbols));
      case 'u':
      default: { // unknown - pick best standard alphabet
        ALPH_T *alph;
        int i;
        int64_t *counts;
        counts = mm_malloc(sizeof(int64_t) * str_len(stats->symbols));
        for (i = 0; i < str_len(stats->symbols); i++) {
          counts[i] = stats->counts.counts[i]->count;
        }
        alph = alph_guess(str_internal(stats->symbols), counts);
        free(counts);
        return alph;
      }
    }
  }
}

static void rc_chain(ALPH_T *alph, int *chain, int *chain_rc, int len) {
  int i, j;
  for (i = 0, j = len - 1; i < j; i++, j--) {
    chain_rc[i] = alph_complement(alph, chain[j]);
    chain_rc[j] = alph_complement(alph, chain[i]);
  }
  if (i == j) {
    chain_rc[i] = alph_complement(alph, chain[i]);
  }
}

static int64_t chain_count(FMKV_STATS_T *stats, int **mapping, int len, int *chain, int *pos_buf) {
  int64_t count;
  int i, j, k, *pos;
  FMKV_COUNT_T *lookup;
  // initilize count
  count = 0;
  // use the passed in buffer
  pos = pos_buf;
  // check that at least one alias was seen for each symbol in this chain
  for (j = 0; j < len; j++) {
    if (mapping[chain[j]][0] == 0) break;
  }
  if (j == len) {
    // iterate over every alias that was seen;
    // we know there is at least one because we checked
    for (j = 0; j < len; j++) pos[j] = 0;
    j = len - 1;
    while (j >= 0) {
      // lookup the chain alias count
      lookup = &(stats->counts);
      for (k = 0; k < len; k++) {
        i = mapping[chain[k]][pos[k]] - 1;
        if (lookup->size <= i) break;
        lookup = lookup->counts[i];
        if (lookup == NULL) break;
      }
      if (k == len) count += lookup->count;
      // next alias
      for (j = len - 1; j >= 0; j--) {
        pos[j]++;
        if (mapping[chain[j]][pos[j]] != 0) {
          break;
        } else {
          pos[j] = 0;
        }
      }
    }
  }
  return count;
}

static void output_stats(FMKV_OPTIONS_T *opts, ALPH_T *alph, FMKV_STATS_T *stats) {
  long double avg;
  FILE *bg_fh, *seqc_fh;
  int i, j, len, **mapping, *chain, *chain_rc, *pos;
  int64_t total, count, len_chains;
  long double pseudo, pseudo_frac, prob;
  
  if (strcmp(opts->bg, "-") == 0) {
    // write to stdout
    bg_fh = stdout;
  } else {
    // open file
    bg_fh = fopen(opts->bg, "w");
    // check it worked
    if (bg_fh == NULL) {
      // report error and exit
      fprintf(stderr, "Failed to open background file \"%s\" for writing: ", opts->bg);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  }
  if (strcmp(opts->fasta, "-") != 0) {
    fprintf(bg_fh, "# %d-order Markov frequencies from file %s\n", opts->order, opts->fasta);
  } else {
    fprintf(bg_fh, "# %d-order Markov frequencies\n", opts->order);
  }
  if (stats->bad) {
    fprintf(bg_fh, "# Warning: unexpected characters were detected, "
        "first instance at byte offset %jd\n", (intmax_t)stats->bad_pos);
  }
  if (stats->total_seqs > 0) {
    avg = (long double)stats->total_slen / (long double)stats->total_seqs;
  } else {
    avg = 0;
  }
  fprintf(bg_fh, "# seqs: %" PRId64 "    min: %" PRId64 "    max: %"
      PRId64 "    avg: %.1Lf    sum: %" PRId64"    alph: %s\n",
      stats->total_seqs, stats->slen_min, stats->slen_max, avg,
      stats->total_slen, alph_name(alph));
  // create a mapping from the alphabet core symbols to the positions in the stats
  mapping = mm_malloc(sizeof(int*) * alph_size_core(alph));
  for (i = 0; i < alph_size_core(alph); i++) {
    mapping[i] = mm_malloc(sizeof(int));
    mapping[i][0] = 0;
  }
  for (i = 0; i < str_len(stats->symbols); i++) {
    char sym;
    int symi;
    sym = str_char(stats->symbols, i);
    symi = alph_index(alph, sym);
    if (symi == -1) {
      fprintf(bg_fh, "# Warning: unexpected '%c' first instance at byte offset %jd\n",
          sym, (intmax_t) stats->symbols_offset[i]);
      continue;
    } else if (symi >= alph_size_core(alph)) {
      continue;
    }
    for (len = 0; mapping[symi][len] != 0;) len++;
    mapping[symi] == mm_realloc(mapping[symi], sizeof(int) * (len + 2));
    mapping[symi][len] = i+1;
    mapping[symi][len+1] = 0;
  }
  // iterate over all chains
  chain = mm_malloc(sizeof(int) * (opts->order + 1));
  if (alph_has_complement(alph) && !opts->norc && !opts->counts) {
    chain_rc = mm_malloc(sizeof(int) * (opts->order + 1));
  } else {
    chain_rc = NULL;
  }
  pos = mm_malloc(sizeof(int) * (opts->order + 1));
  for (len = 1, len_chains = alph_size_core(alph); len <= (opts->order + 1); len++, len_chains *= alph_size_core(alph)) {
    fprintf(bg_fh, "# order %d\n", len - 1);
    total = 0;
    pseudo = 0;
    pseudo_frac = 0;
    // calculate the total core symbol matches for this length
    if (!opts->counts) {
      // reset chain
      for (i = 0; i < len; i++) chain[i] = 0;
      // iterate over every combination of the core alphabet letters for this length
      i = len - 1;
      while (i >= 0) {
        total += chain_count(stats, mapping, len, chain, pos);
        // try next combination of the core symbols
        for (i = len - 1; i >= 0; i--) {
          chain[i]++;
          if (chain[i] < alph_size_core(alph)) break;
          chain[i] = 0;
        }
      }
      // when there are no entries for this length, force a pseudocount to avoid divide by zero
      pseudo = (total == 0 ?  1.0 : opts->pseudo);
      pseudo_frac = pseudo / len_chains;
      // when combining reverse complements double the total
      if (chain_rc) total *= 2;
    }
    // reset chain
    for (i = 0; i < len; i++) chain[i] = 0;
    // iterate over every combination of the core alphabet letters for this length
    i = len - 1;
    while (i >= 0) {
      count = chain_count(stats, mapping, len, chain, pos);
      if (count == 0 && opts->counts) goto next_chain;
      // write out the chain entry
      for (j = 0; j < len; j++) fputc(alph_char(alph, chain[j]), bg_fh);
      fputs(" ", bg_fh);
      if (opts->counts) {
        fprintf(bg_fh, "%" PRId64 "\n", count);
      } else {
        if (chain_rc) {
          rc_chain(alph, chain, chain_rc, len);
          count += chain_count(stats, mapping, len, chain_rc, pos);
        }
        prob = (pseudo_frac + count) / (pseudo + total);
        fprintf(bg_fh, "%9.3Le\n", prob);
      }
next_chain:
      // try next combination of the core symbols
      for (i = len - 1; i >= 0; i--) {
        chain[i]++;
        if (chain[i] < alph_size_core(alph)) break;
        chain[i] = 0;
      }
    }
  }
  // write out the sequence count
  if (opts->seqc != NULL) {
    // open file
    seqc_fh = fopen(opts->seqc, "w");
    // check it worked
    if (seqc_fh == NULL) {
      // report error and exit
      fprintf(stderr, "Failed to open sequence count file \"%s\" for writing: ", opts->seqc);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
    // write out the sequence count
    fprintf(seqc_fh, "%" PRId64 "\n", stats->total_seqs);
    // close the file
    fclose(seqc_fh);
  }

  if (strcmp(opts->bg, "-") != 0 && opts->summary) {
    fclose(bg_fh);
    // when we don't have to print the background to stdout
    // print a report like getsize except without the alphabet
    printf("%" PRId64 " %" PRId64 " %" PRId64 " %.1Lf %" PRId64 "\n", 
      stats->total_seqs, stats->slen_min, stats->slen_max, avg,
      stats->total_slen);
  }
}

void status_update_alarm(int signum) {
  update_status = true;
}

static void start_status_update_alarm() {
  struct sigaction alarm_action;
  struct itimerval msg;

  alarm_action.sa_handler = status_update_alarm;
  sigemptyset(&alarm_action.sa_mask);
  alarm_action.sa_flags = SA_RESTART;

  sigaction(SIGALRM, &alarm_action, NULL);

  msg.it_interval.tv_usec = 500000; // half a second
  msg.it_interval.tv_sec = 0;
  msg.it_value.tv_usec = 500000; // half a second
  msg.it_value.tv_sec = 0;

  setitimer(ITIMER_REAL, &msg, NULL);

  update_status = true;
}

static void stop_status_update_alarm() {
  struct itimerval msg;
  msg.it_interval.tv_usec = 0;
  msg.it_interval.tv_sec = 0;
  msg.it_value.tv_usec = 0;
  msg.it_value.tv_sec = 0;
  setitimer(ITIMER_REAL, &msg, NULL);
}

int main(int argc, char **argv) {
  FMKV_OPTIONS_T opts;
  FMKV_STATS_T stats;
  ALPH_T *alph;
  // run
  process_arguments(&opts, argc, argv);
  if (opts.status) start_status_update_alarm();
  process_fasta(&opts, &stats);
  if (opts.status) stop_status_update_alarm();
  alph = determine_alphabet(&opts, &stats);
  output_stats(&opts, alph, &stats);
  // cleanup
  alph_release(alph);
  cleanup_stats(&opts, &stats);
  cleanup_options(&opts);
  return EXIT_SUCCESS;
}
