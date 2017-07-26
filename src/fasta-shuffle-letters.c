#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "alphabet.h"
#include "config.h"
#include "mtwist.h"
#include "string-builder.h"
#include "ushuffle.h"

/***********************************************************************
 CONSTANTS
 ***********************************************************************/
// different types of symbol
#define BAD 1
#define WHITESPACE 2
#define NEWLINE 4
#define VISIBLE 8
#define CHEVRON 16
// different alphabet types
#define DNA 1
#define PROTEIN 2
#define RNA 3

/***********************************************************************
 STRUCTURES
 ***********************************************************************/
// store the command line options
typedef struct options {
  char *fasta_in;
  char *fasta_out;
  int kmer;
  int copies;
  int line;
  char *tag;
  int alphabet_type;
  char *alphabet_file;
  uint32_t seed;
} OPTIONS_T;

// store the program state
typedef struct state STATE_T;
struct state {
  uint8_t hasher[UCHAR_MAX + 1];
  size_t (*routine)(OPTIONS_T*, STATE_T*, char*, size_t);
  bool nl;
  STR_T *id;
  STR_T *seq;
  ALPH_T *alph;
  FILE *out;
};

/***********************************************************************
 GLOBALS
 ***********************************************************************/
// state of pseudo-random number generator (mersenne twister)
mt_state prng;

/***********************************************************************
 Get positive random numbers from the GLOBAL mersenne twister prng
 between 0 and 2^31-1.
 ***********************************************************************/
static long rand_mt() {
  long value;
  // get unsigned value and convert to positive signed by zeroing the top bit
  value = mts_lrand(&prng) & 0x7FFFFFFF;
  return value;
}

/***********************************************************************
 Display usage message and exit.
 ***********************************************************************/
static void usage(char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n"
    "    fasta-shuffle-letters [options] <sequence file> [<output file>]\n\n"
    "Options:\n"
    "    [-kmer <num>]       size of k-mer to shuffle; default 1.\n"
    "    [-copies <num>]     number of copies to create for each sequence in\n"
    "                         the source; default: 1.\n"
    "    [-dna]              convert symbols assuming DNA.\n"
    "    [-protein]          convert symbols assuming protein.\n"
    "    [-rna]              convert symbols assuming RNA.\n"
    "    [-alph <file>]      convert symbols to their primary representation;\n"
    "                         in the given alphabet; unknown symbols will be\n"
    "                         converted to the alphabet's wildcard.\n"
    "    [-line <num>]       the line length used to output sequences;\n"
    "                         default 100.\n"
    "    [-tag <text>]       the text to append to the name of shuffled\n"
    "                         sequences; default: \"_shuf\".\n"
    "    [-seed <num>]       the seed to the random number generator.\n"
    "    [-help]             display this help message.\n"
    "    \nShuffles FASTA sequences maintaining k-mers.\n"
    "    \nNote: if the sequences contain aliases (for example soft masked sequence)\n"
    "    then specifying the alphabet will ensure a better shuffle.\n"
    "    \nReads the sequence file or stdin if \"-\" is specified.\n"
    "    Writes standard output if a output file is not specified.\n";
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
enum Opts {OPT_ALPH, OPT_DNA, OPT_PROTEIN, OPT_RNA, OPT_KMER, OPT_COPIES,
  OPT_SEED, OPT_LINE, OPT_TAG, OPT_HELP, OPT_VERSION};

/***********************************************************************
 Process command line options
 ***********************************************************************/
static void process_command_line(int argc, char* argv[], OPTIONS_T *config) {
  char *endptr;
  struct option options[] = {
    {"alph", required_argument, NULL, OPT_ALPH},
    {"dna", no_argument, NULL, OPT_DNA},
    {"protein", no_argument, NULL, OPT_PROTEIN},
    {"rna", no_argument, NULL, OPT_RNA},
    {"kmer", required_argument, NULL, OPT_KMER},
    {"copies", required_argument, NULL, OPT_COPIES},
    {"seed", required_argument, NULL, OPT_SEED},
    {"line", required_argument, NULL, OPT_LINE},
    {"tag", required_argument, NULL, OPT_TAG},
    {"help", no_argument, NULL, OPT_HELP},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };
  // Make sure options are set to defaults.
  memset(config, 0, sizeof(OPTIONS_T));
  config->fasta_in = NULL;
  config->fasta_out = NULL;
  config->kmer = 1;
  config->copies = 1;
  config->line = 100;
  config->tag = "_shuf";
  config->alphabet_type = 0;
  config->alphabet_file = NULL;
  config->seed = mts_seed(&prng); // note: referring to GLOBAL here
  // process arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_ALPH:
        config->alphabet_file = optarg;
        break;
      case OPT_DNA:
        config->alphabet_type = DNA;
        break;
      case OPT_PROTEIN:
        config->alphabet_type = PROTEIN;
        break;
      case OPT_RNA:
        config->alphabet_type = RNA;
        break;
      case OPT_KMER:
        config->kmer = (int)strtol(optarg, &endptr, 10);
        if (*endptr != '\0') {
          usage("K-mer \"%s\" was not a number", optarg);
        } else if (config->kmer <= 0) {
          usage("K-mer \"%s\" must be positive", optarg);
        }
        break;
      case OPT_COPIES:
        config->copies = (int)strtol(optarg, &endptr, 10);
        if (*endptr != '\0') {
          usage("Copies \"%s\" was not a number", optarg);
        } else if (config->copies <= 0) {
          usage("Copies \"%s\" must be positive", optarg);
        }
        break;
      case OPT_SEED:
        if (sscanf(optarg, "%" SCNu32, &(config->seed)) == 1) {
          mts_seed32new(&prng, config->seed); // note: referring to GLOBAL here
        } else {
          usage("Seed \"%s\" could not be interpreted as an unsigned 32bit integer", optarg);
        }
        break;
      case OPT_LINE:
        config->line = (int)strtol(optarg, &endptr, 10);
        if (*endptr != '\0') {
          usage("Line \"%s\" was not a number", optarg);
        } else if (config->line < 1) {
          usage("Line \"%s\" must be > 1", optarg);
        }
        break;
      case OPT_TAG:
        config->tag = optarg;
        break;
      case OPT_HELP:
        usage(NULL); // does not return
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case '?':
        usage(NULL);
        break;
      default: // just in case we forget to handle a option
        fprintf(stderr, "Unhandled option %d", opt);
        exit(EXIT_FAILURE);
    }
  }
  // Must have sequence file names
  if (optind >= argc) usage("Sequences not specified");
  config->fasta_in = argv[optind++];
  if (optind < argc) {
    config->fasta_out = argv[optind++];
  }
  if (optind < argc) usage("Unused arguments provided");
}

/***********************************************************************
 Output a sequence
 ***********************************************************************/
static void output_seq(FILE *out, ALPH_T *alph, int line, int kmer, int copies, char *tag, STR_T *id, STR_T *seq) {
  char *shuf, *seq2;
  int i_copy, offset, idx;
  if (alph != NULL) {
    // convert to prime representation inplace
    seq2 = str_internal(seq);
    for (; *seq2 != '\0'; seq2++) {
      idx = alph_index(alph, *seq2);
      if (idx == -1) {
        // convert unknown symbols to wildcard
        *seq2 = alph_wildcard(alph);
      } else {
        // convert known symbols to their prime representation
        *seq2 = alph_char(alph, idx);
      }
    }
  }
  // Create space for shuffled copy of sequence.
  shuf = malloc(sizeof(char) * (str_len(seq) + 1));
  if (shuf == NULL) {
    fprintf(stderr, "Failed to allocate memory for shuffled sequence");
    exit(EXIT_FAILURE);
  }
  shuf[str_len(seq)] = '\0';
  // setup the globals in ushuffle
  ushuffle1(str_internal(seq), str_len(seq), kmer);     
  for (i_copy = 1; i_copy <= copies; i_copy++) {
    // Create a shuffled copy of the sequence from the globals in ushuffle
    ushuffle2(shuf);
    if (copies > 1) {
      fprintf(out, ">%s%s_%d\n", str_internal(id), tag, i_copy);
    } else {
      fprintf(out, ">%s%s\n", str_internal(id), tag);
    }
    for (offset = 0; offset < str_len(seq); offset += line) {
      fprintf(out, "%.*s\n", line, shuf+offset);
    }
  } // copies
  free(shuf);
}

// predeclare routines so they can reference each-other.
static size_t routine_seq_data(OPTIONS_T*, STATE_T*, char*, size_t);
static size_t routine_seq_whitespace(OPTIONS_T*, STATE_T*, char*, size_t);
static size_t routine_seq_name(OPTIONS_T*, STATE_T*, char*, size_t);
static size_t routine_seq_desc(OPTIONS_T*, STATE_T*, char*, size_t);
static size_t routine_find_start(OPTIONS_T*, STATE_T*, char*, size_t);

/***********************************************************************
 routine to read sequence data
 ***********************************************************************/
static size_t routine_seq_data(OPTIONS_T *opts, STATE_T *state, char *buffer, size_t len) {
  size_t i, first_visible;
  for (i = 0; i < len; i++) {
    if ((state->hasher[(unsigned char)buffer[i]] & VISIBLE) == 0) {
      if (i > 0) str_append(state->seq, buffer, i);
      state->routine = routine_seq_whitespace;
      return i;
    }
  }
  str_append(state->seq, buffer, len);
  return len;
}

/***********************************************************************
 routine to skip the non-sequence parts of the sequence (like whitespace)
 ***********************************************************************/
static size_t routine_seq_whitespace(OPTIONS_T *opts, STATE_T *state, char *buffer, size_t len) {
  size_t i, first_visible;
  uint8_t type;
  for (i = 0; i < len; i++) {
    type = state->hasher[(unsigned char)buffer[i]];
    if (state->nl && buffer[i] == '>') {
      output_seq(state->out, state->alph, opts->line, opts->kmer, opts->copies, opts->tag, state->id, state->seq);
      str_clear(state->id);
      str_clear(state->seq);
      state->routine = routine_seq_name;
      return i+1;
    }
    state->nl = false;
    if ((type & NEWLINE) != 0) {
      state->nl = true;
    } else if ((type & VISIBLE) != 0) {
      state->routine = routine_seq_data;
      return i;
    }
  }
  return len;
}

/***********************************************************************
 routine to store the sequence name
 ***********************************************************************/
static size_t routine_seq_name(OPTIONS_T *opts, STATE_T *state, char *buffer, size_t len) {
  size_t i;
  for (i = 0; i < len; i++) {
    if ((state->hasher[(unsigned char)buffer[i]] & WHITESPACE) != 0) {
      if (i > 0) str_append(state->id, buffer, i);
      state->routine = routine_seq_desc;
      return i;
    }
  }
  str_append(state->id, buffer, len);
  return len;
}

/***********************************************************************
 routine to skip over the description
 ***********************************************************************/
static size_t routine_seq_desc(OPTIONS_T *opts, STATE_T *state, char *buffer, size_t len) {
  size_t i;
  for (i = 0; i < len; i++) {
    if ((state->hasher[(unsigned char)buffer[i]] & NEWLINE) != 0) {
      state->nl = true;
      state->routine = routine_seq_data;
      return i + 1;
    }
  }
  return len;
}

/***********************************************************************
 routine to find the first sequence
 ***********************************************************************/
static size_t routine_find_start(OPTIONS_T *opts, STATE_T *state, char *buffer, size_t len) {
  size_t i;
  for (i = 0; i < len;) {
    if (state->nl && buffer[i] == '>') {
      state->routine = routine_seq_name;
      state->nl = false;
      return i + 1;
    } else {
      state->nl = false;
      for (; i < len; i++) {
        if ((state->hasher[(unsigned char)buffer[i]] & NEWLINE) != 0) {
          state->nl = true;
          i++;
          break;
        }
      }
    }
  }
  return len;
}

/***********************************************************************
 process all the sequences
 ***********************************************************************/
static void process(OPTIONS_T *opts, ALPH_T *alph, FILE *in_fh, FILE *out_fh) {
  char *buffer;
  int i;
  size_t bytes_read, offset, consumed, buffer_len;
  STATE_T state;
  size_t (*prev_routine)(OPTIONS_T*, STATE_T*, char*, size_t);
  // init state
  state.routine = routine_find_start;
  // setup the hasher
  for (i = 0; i <= UCHAR_MAX; i++) {
    // set everything to the bad index
    state.hasher[i] = BAD;
  }
  for (i = 33; i <= 126; i++) {
    // zero the visible ascii to allow them to be used
    state.hasher[i] = VISIBLE;
  }
  // disallow chevron
  state.hasher['>'] = CHEVRON;
  // add newlines
  state.hasher['\r'] = NEWLINE | WHITESPACE;
  state.hasher['\n'] = NEWLINE | WHITESPACE;
  // add whitespace
  state.hasher[' '] = WHITESPACE;
  state.hasher['\t'] = WHITESPACE;
  state.hasher['\v'] = WHITESPACE; //vertical tab, rather unlikely
  // track new lines so we know when to look for '>'
  state.nl = true;
  // create storage for sequence parts.
  state.id = str_create(100);
  state.seq = str_create(1000);
  // track out alphabet
  state.alph = alph;
  // track out fh
  state.out = out_fh;
  // allocate 1MB buffer
  buffer_len =  (1 << 20);
  buffer = malloc(sizeof(char) * buffer_len);
  if (buffer == NULL) {
    fprintf(stderr, "Unable to allocate buffer\n");
    exit(EXIT_FAILURE);
  }
  //load data into buffer
  while ((bytes_read = fread(buffer, sizeof(char), buffer_len, in_fh)) > 0) {
    offset = 0;
    while (offset < bytes_read) {
      prev_routine = state.routine;
      consumed = state.routine(opts, &state, buffer+offset, bytes_read - offset);
      if (consumed == 0 && state.routine == prev_routine) {
        fprintf(stderr, "Infinite loop detected!\n");
        abort(); // we'd want to debug this so an abort makes more sense than exit.
      }
      offset += consumed;
    }
    if (feof(in_fh) || ferror(in_fh)) break;
  }
  if (state.routine != routine_find_start) {
    output_seq(out_fh, alph, opts->line, opts->kmer, opts->copies, opts->tag, state.id, state.seq);
  }
  // clean up
  str_destroy(state.id, false);
  str_destroy(state.seq, false);
  free(buffer);
}

/***********************************************************************
 program start point
 ***********************************************************************/
int main(int argc, char **argv) {
  OPTIONS_T config;
  FILE *in_fh, *out_fh;
  ALPH_T *alph;
  // process command line
  process_command_line(argc, argv, &config);
  // set ushuffle's random number function to use the mtwist library with our local copy of state
  set_randfunc(rand_mt);
  // load the alphabet
  if (config.alphabet_file != NULL) {
    alph = alph_load(config.alphabet_file, true);
    if (alph == NULL) exit(EXIT_FAILURE);
  } else {
    switch (config.alphabet_type) {
      case DNA:
        alph = alph_dna();
        break;
      case PROTEIN:
        alph = alph_protein();
        break;
      case RNA:
        alph = alph_rna();
        break;
      default:
        alph = NULL;
    }
  }
  // determine the source file
  if (strcmp(config.fasta_in, "-") == 0) {
    in_fh = stdin;
  } else {
    // open file
    in_fh = fopen(config.fasta_in, "r");
    // check it worked
    if (in_fh == NULL) {
      // report error and exit
      fprintf(stderr, "Failed to open fasta file \"%s\" for reading: ", config.fasta_in);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  }
  // determine the destination file
  if (config.fasta_out == NULL || strcmp(config.fasta_out, "-") == 0) {
    out_fh = stdout;
  } else {
    // open file
    out_fh = fopen(config.fasta_out, "w");
    // check it worked
    if (out_fh == NULL) {
      // report error and exit
      fprintf(stderr, "Failed to open file \"%s\" for writing: ", config.fasta_out);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  }
  // process file
  process(&config, alph, in_fh, out_fh);
  // close destination file
  if (config.fasta_out != NULL && strcmp(config.fasta_out, "-") != 0) {
    fclose(out_fh);
  }
  // close source file
  if (strcmp(config.fasta_in, "-") != 0) {
    fclose(in_fh);
  }
  return EXIT_SUCCESS;
}
