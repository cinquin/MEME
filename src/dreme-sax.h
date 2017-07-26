
#ifndef DREME_SAX_H
#define DREME_SAX_H

#include <libxml/parser.h>

#include <stdarg.h>

/*****************************************************************************
 * Method for obtaining negative dataset enumeration
 ****************************************************************************/
enum dreme_neg {
  DREME_NEG_SHUFFLED, // apply a dinucleotide shuffle to the positive set
  DREME_NEG_FILE  // use a supplied file for the negative dataset
};
typedef enum dreme_neg DREME_NEG_EN;

/*****************************************************************************
 * Method for obtaining negative dataset enumeration
 ****************************************************************************/
enum dreme_strands {
  DREME_STRANDS_BOTH, // apply a dinucleotide shuffle to the positive set
  DREME_STRANDS_GIVEN,// use a supplied file for the negative dataset
  DREME_STRANDS_NONE
};
typedef enum dreme_strands DREME_STRANDS_EN;

/*****************************************************************************
 * Background sources enumeration
 ****************************************************************************/
enum dreme_bg {
  DREME_BG_FROM_DATASET, // use negative dataset frequencies
  DREME_BG_FROM_FILE  // use pre-calculated frequencies from a file
};
typedef enum dreme_bg DREME_BG_EN;

/*****************************************************************************
 * Stopping reasons enumeration
 ****************************************************************************/
enum dreme_stop {
  DREME_STOP_EVALUE, // couldn't find a motif with a smaller evalue
  DREME_STOP_COUNT, // already found the expected number of motifs
  DREME_STOP_TIME // dreme has run longer than the maximum expected time
};
typedef enum dreme_stop DREME_STOP_EN;

struct dreme_io_xml_callbacks {

/*****************************************************************************
 * Error message handler
 ****************************************************************************/
  void (*error)(void * ctx, const char *format, va_list args);

/*****************************************************************************
 * dreme
 * contains: model, motifs, run_time
 *
 *  release         the release date.
 *  version         the program version.
 ****************************************************************************/
  void (*start_dreme)(void *ctx, int major_version, int minor_version, 
      int patch_version, char *release_date);

/*****************************************************************************
 * /dreme
 * contains: model, motifs, run_time
 ****************************************************************************/
  void (*end_dreme)(void *ctx);

/*****************************************************************************
 * dreme > model
 * contains: command_line, positives, negatives, background, stop, ngen, 
 *          add_pv_thresh, seed, host, when, description
 ****************************************************************************/
  void (*start_model)(void *ctx);

/*****************************************************************************
 * dreme > /model
 * contains: command_line, positives, negatives, background, stop, ngen, 
 *          add_pv_thresh, seed, host, when, description?
 ****************************************************************************/
  void (*end_model)(void *ctx);

/*****************************************************************************
 * dreme > model > command_line
 * the command-line used to run the program.
 ****************************************************************************/
  void (*handle_command_line)(void *ctx, char *command_line);

/*****************************************************************************
 * dreme > model > positives
 *
 *  name            the name of the positive input set
 *  count           the number of sequences in the positive input set
 *  file            the file containing the positive input set
 *  last_mod_date   the time the positive input set was last modified
 ****************************************************************************/
  void (*handle_positives) (void *ctx, char *dataset_name, long sequence_count, 
      char *positives_file_path, char *positives_file_last_modified_date);

/*****************************************************************************
 * dreme > model > negatives
 *
 *  name            the name of the negative dataset
 *  count           the number of sequences in the negative dataset
 *  from            the source of the negative dataset (eg shuffled positives)
 *  file            the file containing the negative dataset (optional)
 *  last_mod_date   the last modified date of the file (optional)
 ****************************************************************************/
  void (*handle_negatives) (void *ctx, char *dataset_name, long sequence_count, 
      DREME_NEG_EN negatives_source, char *negatives_file_path, 
      char *negatives_file_last_modified_date);

/*****************************************************************************
 * dreme > model > alphabet
 *
 *  name            is the name of the alphabet
 *  extends_flags   is the alphabet being extended or 0
 ****************************************************************************/
  void (*start_alphabet) (void *ctx, char *name, int extends_flags);

/*****************************************************************************
 * dreme > model > /alphabet
 *
 *  name            is the name of the alphabet
 ****************************************************************************/
  void (*end_alphabet) (void *ctx);

/*****************************************************************************
 * dreme > model > alphabet > letter
 *
 *  id              the identifier used to represent the symbol as XML attributes
 *  symbol          the symbol
 *  aliases         additional symbols that resolve to the same meaning
 *  complement      the complement symbol in a double stranded alphabet
 *  equals          the comprising symbols
 *  name            the name of the negative dataset
 *  colour          the colour in RGB triplet (first 3 bytes)
 ****************************************************************************/
  void (*handle_alphabet_letter)(void *ctx, char *id, char symbol,
      char *aliases, char complement, char *equals, char *name, int color);

/*****************************************************************************
 * dreme > model > strands
 *
 *  strands         how is dreme handling the strands
 ****************************************************************************/
  void (*handle_strands) (void *ctx, DREME_STRANDS_EN strands);

/*****************************************************************************
 * dreme > model > background
 *
 *  nfreqs          is the count of frequencies (same as alphabet core size)
 *  freqs           is the array of frequencies (ordered as in alphabet core)
 *  from            from the negative dataset or a background file
 *  file            the background file (optional)
 *  last_mod_date   the last modified date of the background file (optional)
 ****************************************************************************/
  void (*handle_background) (void *ctx, int nfreqs, double* freqs, 
      DREME_BG_EN background_source, char *background_file_path, 
      char *background_file_last_modified_date);

/*****************************************************************************
 * dreme > model > stop
 *
 *  evalue          the stopping evalue (in log base 10) (optional).
 *  count           the stopping count (optional).
 *  time            the stopping time (optional).
 ****************************************************************************/
  void (*handle_stop) (void *ctx, double *log10_max_motif_evalue, 
      int *max_motif_count, int *max_time_elapsed);

/*****************************************************************************
 * dreme > model > ngen
 * the number of generations to check (or something like that).
 ****************************************************************************/
  void (*handle_ngen) (void *ctx, int num_gen);

/*****************************************************************************
 * dreme > model > add_pv_thresh
 * the p-value threshold (in log base 10) for adding a word to the list of 
 * possible motif seed values (this description is probably wrong).
 ****************************************************************************/
  void (*handle_add_pv_thresh)(void *ctx, double log10_add_pv_thresh);

/*****************************************************************************
 * dreme > model > seed
 * the seed value used by the pseudo-random number generator.
 ****************************************************************************/
  void (*handle_seed)(void *ctx, long seed);

/*****************************************************************************
 * dreme > model > host
 * the name of the computer which ran DREME.
 ****************************************************************************/
  void (*handle_host)(void *ctx, char *host_name);

/*****************************************************************************
 * dreme > model > when
 * the time DREME was run. 
 ****************************************************************************/
  void (*handle_when)(void *ctx, char *start_date);

/*****************************************************************************
 * dreme > model > description
 * an optional description of the experiment.
 ****************************************************************************/
  void (*handle_description)(void *ctx, char *experiment_description);

/*****************************************************************************
 * dreme > motifs
 * contains: motif*
 ****************************************************************************/
  void (*start_motifs)(void *ctx);

/*****************************************************************************
 * dreme > /motifs
 * contains: motif*
 ****************************************************************************/
  void (*end_motifs)(void *ctx);

/*****************************************************************************
 * dreme > motifs > motif
 * contains: pos+, match*
 *
 *  id                the identifier used by DREME
 *  alt               the alternate identifier used by DREME
 *  seq               the DNA iupac sequence representing the motif.
 *  length            the length of the motif
 *  nsites            the number of sites used to create the motif
 *  p                 the number of sequences in the positive set with the motif
 *  n                 the number of sequences in the negative set with the motif
 *  pvalue            the pvalue of the motif after erasing (in log base 10)
 *  evalue            the evalue of the motif after erasing (in log base 10)
 *  unerased_evalue   the evalue of the motif without erasing (in log base 10)
 ****************************************************************************/
  void (*start_motif)(void *ctx, char *identifier, char *alt, char *iupac, int length,
      long sites, long positives_with_motif, long negatives_with_motif,
      double log10_pvalue, double log10_evalue, double log10_unerased_evalue);

/*****************************************************************************
 * dreme > motifs > /motif
 * contains: pos+, match*
 ****************************************************************************/
  void (*end_motif)(void *ctx);

/*****************************************************************************
 * dreme > motifs > motif > pos
 *
 *  position          index of the motif position
 *  nfreqs            number of frequencies
 *  freqs             array of frequencies
 ****************************************************************************/
  void (*handle_pos)(void *ctx, int position, int nfreqs, double* freqs);

/*****************************************************************************
 * dreme > motifs > motif > match
 *
 *  seq               word which matches motif
 *  p                 number of positive sequences which contain this word
 *  n                 number of negative sequences which contain this word
 *  pvalue            the pvalue of this word (in log base 10)
 *  evalue            the evalue of this word (in log base 10)
 ****************************************************************************/
  void (*handle_match)(void *ctx, char *matching_word, long positives_with_word,
      long negatives_with_word, double log10_pvalue, double log10_evalue);

/*****************************************************************************
 * dreme > run_time
 *
 *  cpu             the time the program was running on a cpu.
 *  real            the real world time the program was running.
 *  stop            the reason dreme stopped.
 ****************************************************************************/
  void (*handle_run_time)(void *ctx, double cpu_time, double real_time, 
      DREME_STOP_EN stopping_reason);
};
typedef struct dreme_io_xml_callbacks DREME_IO_XML_CALLBACKS_T;

/*****************************************************************************
 * Register handlers on the xmlSAXHandler structure
 ****************************************************************************/
void register_dreme_io_xml_sax_handlers(xmlSAXHandler *handler);

/*****************************************************************************
 * Creates the data to be passed to the SAX parser
 ****************************************************************************/
void* create_dreme_io_xml_sax_context(void *user_data, DREME_IO_XML_CALLBACKS_T *callbacks);

/*****************************************************************************
 * Destroys the data used by the SAX parser
 ****************************************************************************/
void destroy_dreme_io_xml_sax_context(void *context);

#endif
