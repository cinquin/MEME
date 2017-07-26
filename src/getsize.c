/***********************************************************************
*                                                                      *
*       MEME
*       Copyright 1994--2014, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
/* getsize.c */
/*      
        getsize <datafile> [options]

        Reads in a sequence file (Pearson/FASTA format).
        Prints the number of sequences, min L, max L, mean L, total residues
        and letters in alphabet used.

        Type "getsize" to see options.
*/

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
        
#include "alphabet.h"
#include "config.h"
#include "read_sequence.h"
#include "hash_table.h"
#include "utils.h"

// JIJ: the hash size was 100000 but I changed it to 100003 to make it prime
// as that should reduce hash collisions
#define DATA_HASH_SIZE 100003

VERBOSE_T verbosity = NORMAL_VERBOSE;

typedef struct options {
  char *datafile;
  char *alphfile;
  int bfile;                        // do not print bfile format
  bool l_only;
  bool print_duplicates;
  bool print_frequencies;
  bool print_table;
  bool xlate_dna;
  bool print_codons;
} OPTIONS_T;

static void usage(char *format, ...) {
  va_list argp;

  char *usage = 
"Usage:\n"
"        getsize <datafile> [options]\n"
"\n"
"        <datafile>              file containing sequences in FASTA format\n"
"        [-dna]                  print DNA frequencies in bfile format\n"
"        [-rna]                  print RNA frequencies in bfile format\n"
"        [-prot]                 print protein frequencies in bfile format\n"
"        [-alph <file>]          use the specified alphabet\n"
"        [-f]                    print letter frequencies as C array\n"
"        [-ft]                   print letter frequencies as latex table\n"
"        [-l]                    just print the length of each sequence\n"
"        [-nd]                   do not print warnings about duplicate sequences\n"
"        [-x]                    translate DNA in six frames and print freqs as C arrays\n"
"        [-codons]               print frame freqs and frame0 codon usage as C arrays\n"
"\n"
"        Measure statistics of a FASTA file.\n";
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
enum Opts {OPT_LENGTH, OPT_NO_DUPLICATES, OPT_DNA, OPT_PROTEIN, OPT_RNA,
  OPT_FREQUENCIES, OPT_FREQUENCIES_TABLE, OPT_TRANSLATE_DNA, OPT_CODONS,
  OPT_ALPH, OPT_SHUFFLE, OPT_COPIES, OPT_VERSION};

/***********************************************************************
 Process command line options
 ***********************************************************************/
static void process_command_line(int argc, char* argv[], OPTIONS_T *options) {
  struct option centrimo_options[] = {
    {"length", no_argument, NULL, OPT_LENGTH},
    {"nd", no_argument, NULL, OPT_NO_DUPLICATES},
    {"dna", no_argument, NULL, OPT_DNA},
    {"protein", no_argument, NULL, OPT_PROTEIN},
    {"rna", no_argument, NULL, OPT_RNA},
    {"f", no_argument, NULL, OPT_FREQUENCIES},
    {"ft", no_argument, NULL, OPT_FREQUENCIES_TABLE},
    {"x", no_argument, NULL, OPT_TRANSLATE_DNA},
    {"codons", no_argument, NULL, OPT_CODONS},
    {"alph", required_argument, NULL, OPT_ALPH},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };
  // Make sure options are set to defaults.
  memset(options, 0, sizeof(OPTIONS_T));
  options->datafile = NULL;
  options->alphfile = NULL;
  options->bfile = 0; // not generating background
  options->l_only = false;
  options->print_duplicates = true;
  options->print_frequencies = false;
  options->print_table = false;
  options->xlate_dna = false;
  options->print_codons = false;
  // process arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", centrimo_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_LENGTH:
        options->l_only = true;
        break;
      case OPT_NO_DUPLICATES:
        options->print_duplicates = false;
        break;
      case OPT_DNA:
        options->bfile = 1;
        break;
      case OPT_PROTEIN:
        options->bfile = 2;
        break;
      case OPT_RNA:
        options->bfile = 3;
        break;
      case OPT_FREQUENCIES:
        options->print_frequencies = true;
        break;
      case OPT_FREQUENCIES_TABLE:
        options->print_table = true;
        break;
      case OPT_TRANSLATE_DNA:
        options->xlate_dna = true;
        break;
      case OPT_CODONS:
        options->xlate_dna = true;
        options->print_frequencies = true;
        options->print_codons = true;
        break;
      case OPT_ALPH:
        options->alphfile = optarg;
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
  if (optind >= argc) usage("Sequences not specified");
  options->datafile = argv[optind++];
  if (optind < argc) usage("Unused arguments provided");
}

/**********************************************************************/
/*
        main 
*/
/**********************************************************************/

int main(
  int argc,
  char *argv[]
)
{
  OPTIONS_T options;
  long i, j, k, count;
  FILE *data_file;
  int n_samples; 
  long max_slength, min_slength, length;
  char *sample_name, *id, *seq;
  HASH_TABLE ht_seq_names;              // hash of dataset seq names
  double *codons;                       // counts of used codons
  double **alphabet;                    // counts of used letters (per frame)
  double *total_res;                    // total letters per frame
  char *letters = NULL;                 // string of used letters
  ALPH_T *alph = NULL;
  XLATE_T *xlate = NULL;
  ALPH_T *src_alph, *dest_alph;
  bool is_generic;

  process_command_line(argc, argv, &options);

  // Setup hashing function for encoding strings as integers.
  // If translating from DNA to protein, input alphabet must be DNAB;
  // otherwise it may be all 26 letters plus asterisk.
  alph = NULL;
  src_alph = NULL;
  dest_alph = NULL;
  is_generic = false;
  if (options.xlate_dna) {
    xlate = xlate_dna2protein();
    src_alph = xlate_src_alph(xlate);
    dest_alph = xlate_dest_alph(xlate);
  } else if (options.bfile == 1) {// get DNA frequencies
    alph = alph_dna();
    src_alph = alph;
    dest_alph = alph;
  } else if (options.bfile == 2) {// get protein frequencies
    alph = alph_protein();
    src_alph = alph;
    dest_alph = alph;
  } else if (options.bfile == 3) {// get RNA frequencies
    alph = alph_rna();
    src_alph = alph;
    dest_alph = alph;
  } else if (options.alphfile != NULL) {
    alph = alph_load(options.alphfile, true);
    if (alph == NULL) exit(EXIT_FAILURE);
    src_alph = alph;
    dest_alph = alph;
  } else {
    alph = alph_generic("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789*-.");
    src_alph = alph;
    dest_alph = alph;
    is_generic = true;
  }
  Resize(letters, MAXASCII, char);

  // create a hash table of sequence names
  ht_seq_names = hash_create(DATA_HASH_SIZE);

  // open data file
  if (options.datafile == NULL) {
    fprintf(stderr, "You must specify a data file or 'stdin'\n");
    exit(1);
  } else if ((strcmp(options.datafile, "stdin") == 0) || 
             (strcmp(options.datafile, "-") == 0)) {
    data_file = stdin;
  } else {
    data_file = fopen(options.datafile, "r");
    if (data_file == NULL) {
      fprintf(stderr, "Cannot open file '%s'\n", options.datafile);
      exit(1);
    }
  }

  // initialize counts of letters and codons used
  codons = NULL;
  if (xlate != NULL) {
    // work out how many frames
    count = xlate_src_nsyms(xlate) * (alph_has_complement(src_alph) ? 2 : 1);
    // initilize counts for each translated residue at each frame
    alphabet = mm_malloc(sizeof(double*) * count);
    for (i = 0; i < count; i++) {
      alphabet[i] = mm_malloc(sizeof(double) * (alph_size_full(dest_alph) + 1));
      for (j = 0; j <= alph_size_full(dest_alph); j++) alphabet[i][j] = 0;
    }
    // initilize counts for the total count of translated residues for each frame
    total_res = mm_malloc(sizeof(double) * count);
    for (i = 0; i < count; i++) total_res[i] = 0;
    // initilize counts for the frame 0 residue sets
    count = pow(alph_size_full(src_alph) + 1, xlate_src_nsyms(xlate));
    codons = mm_malloc(sizeof(double) * count);
    for (i = 0; i < count; i++) codons[i] = 0;
  } else {
    // when not translating there is only one frame
    // initilize counts for each residue
    alphabet = mm_malloc(sizeof(double*));
    alphabet[0] = mm_malloc(sizeof(double) * (alph_size_full(alph) + 1));
    for (j = 0; j <= alph_size_full(alph); j++) alphabet[0][j] = 0;
    // initilize counts for the total count of residue
    total_res = mm_malloc(sizeof(double));
    total_res[0] = 0;
  }

  // initialize maximum length of sequences
  max_slength = 0;
  min_slength = 1e12;

  n_samples = 0; // no samples yet

  while (read_sequence(src_alph, data_file, &sample_name, &id, &seq, &length)) {
    // Skip weights
    if (strcmp(sample_name, "WEIGHTS")==0) continue;

    if (options.print_duplicates) { // ignore duplicate (same sample name) sequences
      if (hash_lookup_str(sample_name, ht_seq_names) != NULL) {
        fprintf(stderr, "Duplicate sequence: %s.\n", 
          sample_name);
        // free up unneeded space
        myfree(sample_name);
        myfree(id);
        myfree(seq);
        continue;
      }
      hash_insert_str(sample_name, ht_seq_names); // put name in table
    } else {
      myfree(sample_name);
    }

    // Count letters used in sequence.
    if (!options.l_only) { 
      if (xlate != NULL) { // translate 
        int end = length - (xlate_src_nsyms(xlate) - 1);
        int strands = (alph_has_complement(xlate_src_alph(xlate)) ? 2 : 1);
        for (i = 0; i < strands; i++) { // positive then negative strands
          if (i) invcomp_seq(xlate_src_alph(xlate), seq, length); // negative strand
          for (j = 0; j < xlate_src_nsyms(xlate); j++) { // frame
            int f = j + (xlate_src_nsyms(xlate) * i);
            for (k = j; k < end; k += xlate_src_nsyms(xlate)) {
              // convert the sequence into index for translation
              int codon_pos = xlate_pos(xlate, false, seq+k);
              alphabet[f][xlate_index2(xlate, codon_pos) + 1]++;
              total_res[f]++;
              if (options.print_codons && f == 0) { // frame 0
                codons[codon_pos]++;
              }
            }
          } // frame
        } // strand
      } else { // don't translate
        for (i = 0; i < length; i++) alphabet[0][alph_index(alph, seq[i]) + 1]++;
        total_res[0] += length;
      } // xlate_dna
    } // !l_only

    // free up unneeded space
    myfree(id);
    myfree(seq);

    // record maximum length of actual sequences
    max_slength = MAX(length, max_slength);
    min_slength = MIN(length, min_slength);

    if (options.l_only) printf("%ld\n", length);

    n_samples++; // number of non-duplicate sequences
    if (options.print_frequencies && n_samples%1000 == 0) fprintf(stderr, "%d\r", n_samples);
  } // read_sequence

  /* print results */
  if (!options.l_only) {

    /* make string of letters used */
    if (xlate != NULL) {
      // calculate frame count
      count = xlate_src_nsyms(xlate) * (alph_has_complement(src_alph) ? 2 : 1);
      for (i = 1, j = 0; i <= alph_size_full(dest_alph); i++) {
        for (k = 0; k < count; k++) {
          if (alphabet[k][i]) {
            letters[j++] = alph_char(dest_alph, i - 1);
            break;
          }
        }
      }
      letters[j] = '\0';
    } else {
      for (i = 1, j = 0; i <= alph_size_full(alph); i++) {
        if (alphabet[0][i]) {
          letters[j++] = alph_char(alph, i - 1);
        }
      }
      letters[j] = '\0';
    }

    if (!options.bfile) {
      printf("%d %ld %ld %10.1f %.0f %s\n", n_samples, min_slength, max_slength,
        total_res[0]/n_samples, total_res[0], letters);
    }

    /*
      Print frequencies in bfile format
    */
    if (options.bfile) {
      double tot; // total for core letters
      printf("# 0-order Markov frequencies from file %s\n", options.datafile);
      for (i = 0, tot = 0; i < alph_size_core(dest_alph); i++) {
        tot += alphabet[0][i + 1];
      }
      for (i = 0; i < alph_size_core(dest_alph); i++) {                       /* letter */
        double freq;
        freq = tot ? alphabet[0][i + 1] / tot : 0;
        printf("%c %.3e\n", alph_char(dest_alph, i), freq);
      }
    }

    /*
      Print letter frequencies as C arrays.
    */
    if (options.print_frequencies) {
      int nframes;
      if (xlate != NULL) {
        nframes = xlate_src_nsyms(xlate) * (alph_has_complement(src_alph) ? 2 : 1);
      } else {
        nframes = 1;
      }
      for (i = 0; i < nframes; i++) {                       /* frame */
        if (nframes) {
          printf("  double frame%ld[] = {\n", i);
        } else {
          printf("  double freq[] = {\n");
        }
        if (is_generic) {
          for (j = 0; letters[j]; j++) { // letters used
            char c = letters[j]; // letter
            int k = alph_index(dest_alph, c); // index
            printf("  %11.8f /* %c */,\n", alphabet[i][k + 1] / total_res[i], c);
          }
        } else {
          for (j = 0; j < alph_size_core(dest_alph); j++) {
            printf("  %11.8f /* %c */,\n", alphabet[i][j + 1] / total_res[i], alph_char(dest_alph, j));
          }
        }
        printf("  };\n");
      } /* frame */
      if (options.print_codons) {
        int *indexes;
        char *codon;
        indexes = mm_malloc(sizeof(int) * xlate_src_nsyms(xlate));
        codon = mm_malloc(sizeof(int) * (xlate_src_nsyms(xlate) + 1));
        for (i = 0; i < xlate_src_nsyms(xlate); i++) {
          indexes[i] = 0;
          codon[i] = alph_char(src_alph, 0);
        }
        codon[xlate_src_nsyms(xlate)] = '\0';
        printf("  double fcodon[] = {\n");              /* start C array */
        while (i >= 0) {
          printf("  %11.8f /* %s */,\n", 
              codons[xlate_pos(xlate, false, codon)]/total_res[0], codon);
          // increment
          for (i = xlate_src_nsyms(xlate) - 1; i >= 0; i--) {
            indexes[i]++;
            if (indexes[i] < alph_size_core(src_alph)) {
              codon[i] = alph_char(src_alph, indexes[i]);
              break;
            }
            indexes[i] = 0;
            codon[i] = alph_char(src_alph, 0);
          }
        }
        printf("  };\n");
      } /* print_codons */
    } /* print_frequencies */

    /*
      Print letter frequencies as latex table.
    */
    if (options.print_table) {
      ARRAY_T *freqs;
      int nframes; // number of frames
      if (xlate != NULL) {
        nframes = xlate_src_nsyms(xlate) * (alph_has_complement(src_alph) ? 2 : 1);
        freqs = get_nrdb_frequencies(dest_alph, NULL);
        calc_ambigs(dest_alph, false, freqs);
      } else {
        dest_alph = alph;
        nframes = 1;
        freqs = NULL;
      }
      for (j = 0; letters[j]; j++) { // letters used
        char c;
        int k;
        c = letters[j]; // letter
        k = alph_index(dest_alph, c); // index
        if (options.xlate_dna) {
          printf("%c & %7.4f", c, get_array_item(k, freqs));
        } else {
          printf("%c", c);
        }
        for (i = 0; i < nframes; i++) {                     /* frame */
          printf(" & %7.4f", alphabet[i][k + 1] / total_res[i]);
        } /* frame */
        printf(" \\\\\n");
      } /* letters used */
      free_array(freqs);
    } /* print_table */
  }

  return(0); 
} /* getsize */

