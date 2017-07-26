/********************************************************************
 * FILE: fimo.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy Bailey
 * CREATE DATE: 12/17/2004
 * PROJECT: MEME suite
 * COPYRIGHT: 2004-2007, UW
 ********************************************************************/

#define DEFINE_GLOBALS

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "alphabet.h"
#include "cisml.h"
#include "config.h"
#include "dir.h"
#include "fasta-io.h"
#include "fimo.h"
#include "fimo-output.h"
#include "heap.h"
#include "io.h"
#include "matrix.h"
#include "motif-in.h"
#include "pssm.h"
#include "prior-reader-from-psp.h"
#include "prior-reader-from-wig.h"
#include "reservoir.h"
#include "seq-reader-from-fasta.h"
#include "simple-getopt.h"
#include "string-list.h"
#include "utils.h"
#include "wiggle-reader.h"

char* program_name = "fimo";
VERBOSE_T verbosity = NORMAL_VERBOSE;

const char *threshold_type_to_string(THRESHOLD_TYPE type) {
  switch(type) {
    case PV_THRESH:
      return "p-value";
      break;
    case QV_THRESH:
      return "q-value";
      break;
    default:
      return "invalid";
      break;
  }
}

static DATA_BLOCK_READER_T *get_psp_reader(
  const char * prior_filename,
  double default_prior,
  BOOLEAN_T parse_genomic_coord
) {
  DATA_BLOCK_READER_T *psp_reader = NULL;

  // If suffix for file is '.wig' assume it's a wiggle file
  size_t len = strlen(prior_filename);
  if (strcmp(prior_filename + len - 4, ".wig") == 0) {
    psp_reader = new_prior_reader_from_wig(
      prior_filename, 
      default_prior
    );
  }
  else {
    psp_reader = new_prior_reader_from_psp(
        parse_genomic_coord,
        prior_filename
      );
  }
  
  return psp_reader;
}

/***********************************************************************
  Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options(FIMO_OPTIONS_T options) {
  myfree(options.command_line);
  myfree(options.html_path);
  myfree(options.text_path);
  myfree(options.gff_path);
  myfree(options.xml_path);
  myfree(options.cisml_path);
  free_string_list(options.selected_motifs);
  alph_release(options.alphabet);
}

/***********************************************************************
  Process command line options
 ***********************************************************************/
static FIMO_OPTIONS_T process_fimo_command_line(
  int argc,
  char* argv[]
) {

  FIMO_OPTIONS_T options;

  // Define command line options.
  cmdoption const fimo_options[] = {
    {"alpha", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"max-stored-scores", REQUIRED_VALUE},
    {"max-strand", NO_VALUE},
    {"motif", REQUIRED_VALUE},
    {"motif-pseudo", REQUIRED_VALUE},
    {"norc", NO_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"no-qvalue", NO_VALUE},
    {"parse-genomic-coord", NO_VALUE},
    {"psp", REQUIRED_VALUE},
    {"prior-dist", REQUIRED_VALUE},
    {"qv-thresh", NO_VALUE},
    {"text", NO_VALUE},
    {"skip-matched-sequence", NO_VALUE},
    {"thresh", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE},
    {"pval-lookup", REQUIRED_VALUE} // This option is hidden from users.
  };
  const int num_options = sizeof(fimo_options) / sizeof(cmdoption);

  // Define the usage message.
  options.usage =
    "Usage: fimo [options] <motif file> <sequence file>\n"
    "\n"
    "   Options:\n"
    "     --alpha <double> (default 1.0)\n"
    "     --bgfile <background file> (DNA and protein use NRDB by default)\n"
    "     --max-stored-scores <int> (default=100000)\n"
    "     --max-strand\n"
    "     --motif <id> (default=all)\n"
    "     --motif-pseudo <float> (default=0.1)\n"
    "     --no-qvalue\n"
    "     --norc\n"
    "     --o <output dir> (default=fimo_out)\n"
    "     --oc <output dir> (default=fimo_out)\n"
    "     --parse-genomic-coord\n"
    "     --psp <PSP filename> (default none)\n"
    "     --prior-dist <PSP distribution filename> (default none)\n"
    "     --qv-thresh\n"
    "     --skip-matched-sequence\n"
    "     --text\n"
    "     --thresh <float> (default = 1e-4)\n"
    "     --verbosity [1|2|3|4] (default 2)\n"
    "     --version (print the version and exit)\n"
    "\n"
    "   When scanning with a single motif use \'-\' for <sequence file> to\n"
    "     read the database from standard input.\n"
    "   Use \'--bgfile --motif--\' to read the background from the motif file.\n"
    "   Use \'--bgfile --uniform--\' to use a uniform background.\n"
    "\n";

  int option_index = 0;

  /* Make sure various options are set to NULL or defaults. */
  options.allow_clobber = TRUE;
  options.compute_qvalues = TRUE;
  options.max_strand = FALSE;
  options.parse_genomic_coord = FALSE;
  options.threshold_type = PV_THRESH;
  options.text_only = FALSE;
  options.scan_both_strands = TRUE;
  options.skip_matched_sequence = FALSE;

  options.bg_filename = NULL;
  options.command_line = NULL;
  options.meme_filename = NULL;
  options.output_dirname = "fimo_out";
  options.psp_filename = NULL;
  options.prior_distribution_filename = NULL;
  options.seq_filename = NULL;

  options.max_stored_scores = 100000;


  options.alpha = 1.0;
  options.pseudocount = 0.1;
  options.pseudocount = 0.1;
  options.output_threshold = 1e-4;

  options.selected_motifs = new_string_list();
  options.pval_lookup_filename = NULL;
  verbosity = 2;

  simple_setopt(argc, argv, num_options, fimo_options);

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
    if (strcmp(option_name, "bgfile") == 0){
      options.bg_filename = option_value;
    }
    else if (strcmp(option_name, "psp") == 0){
      options.psp_filename = option_value;
    }
    else if (strcmp(option_name, "prior-dist") == 0){
      options.prior_distribution_filename = option_value;
    }
    else if (strcmp(option_name, "alpha") == 0) {
      options.alpha = atof(option_value);
    }
    else if (strcmp(option_name, "max-stored-scores") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options.max_stored_scores = (int)atof(option_value);
    }
    else if (strcmp(option_name, "max-strand") == 0) {
      // Turn on the max-strand option
      options.max_strand = TRUE;
    }
    else if (strcmp(option_name, "motif") == 0){
      if (options.selected_motifs == NULL) {
        options.selected_motifs = new_string_list();
      }
      add_string(option_value, options.selected_motifs);
    }
    else if (strcmp(option_name, "motif-pseudo") == 0){
      options.pseudocount = atof(option_value);
    }
    else if (strcmp(option_name, "norc") == 0){
      options.scan_both_strands = FALSE;
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
    else if (strcmp(option_name, "parse-genomic-coord") == 0){
      options.parse_genomic_coord = TRUE;
    }
    else if (strcmp(option_name, "thresh") == 0){
      options.output_threshold = atof(option_value);
    }
    else if (strcmp(option_name, "qv-thresh") == 0){
      options.threshold_type = QV_THRESH;
    }
    else if (strcmp(option_name, "no-qvalue") == 0){
      options.compute_qvalues = FALSE;
    }
    else if (strcmp(option_name, "skip-matched-sequence") == 0){
      options.skip_matched_sequence = TRUE;
      options.text_only = TRUE;
    }
    else if (strcmp(option_name, "text") == 0){
      options.text_only = TRUE;
    }
    else if (strcmp(option_name, "verbosity") == 0){
      verbosity = atoi(option_value);
    }
    else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    }
    else if (strcmp(option_name, "pval-lookup") == 0) {
      options.pval_lookup_filename = option_value;
    }
  }

  // Check that positiion specific priors options are consistent
  if (options.psp_filename != NULL 
      && options.prior_distribution_filename == NULL) {
    die(
      "Setting the --psp option requires that the"
      " --prior-dist option be set as well.\n"
    );
  }
  if (options.psp_filename == NULL 
      && options.prior_distribution_filename != NULL) {
    die(
      "Setting the --prior-dist option requires that the"
      " --psp option be set as well.\n");
  }

  // Check that qvalue options are consistent
  if (options.compute_qvalues == FALSE && options.threshold_type == QV_THRESH) {
    die("The --no-qvalue option cannot be used with the --qv-thresh option");
  }

  // Turn off q-values if text only.
  if (options.text_only == TRUE) {
    if (options.threshold_type == QV_THRESH) {
      die("The --text option cannot be used with the --qv-thresh option");
    }
    if (options.compute_qvalues) {
      fprintf(stderr, "Warning: text mode turns off computation of q-values\n");
    }
    options.compute_qvalues = FALSE;
  }

  // Must have sequence and motif file names
  if (argc != option_index + 2) {
    fprintf(stderr, "%s", options.usage);
    exit(EXIT_FAILURE);
  }
  // Record the command line
  options.command_line = get_command_line(argc, argv);

  // Record the input file names
  options.meme_filename = argv[option_index];
  option_index++;
  options.seq_filename = argv[option_index];
  option_index++;

  // Set up path values for needed stylesheets and output files.
  options.HTML_FILENAME = "fimo.html";
  options.TEXT_FILENAME = "fimo.txt";
  options.GFF_FILENAME = "fimo.gff";
  options.XML_FILENAME = "fimo.xml";
  options.CISML_FILENAME = "cisml.xml";
  options.html_path = make_path_to_file(options.output_dirname, options.HTML_FILENAME);
  options.text_path = make_path_to_file(options.output_dirname, options.TEXT_FILENAME);
  options.gff_path = make_path_to_file(options.output_dirname, options.GFF_FILENAME);
  options.xml_path = make_path_to_file(options.output_dirname, options.XML_FILENAME);
  options.cisml_path = make_path_to_file(options.output_dirname, options.CISML_FILENAME);

  return options;

}

/**********************************************
* Read the motifs from the motif file.
**********************************************/
static void fimo_read_motifs(
  FIMO_OPTIONS_T *options, 
  ARRAYLST_T **motifs, 
  ARRAY_T **bg_freqs) {

  MREAD_T *mread;

  mread = mread_create(options->meme_filename, OPEN_MFILE);
  mread_set_bg_source(mread, options->bg_filename);
  mread_set_pseudocount(mread, options->pseudocount);

  *motifs = mread_load(mread, NULL);
  options->alphabet = alph_hold(mread_get_alphabet(mread));
  *bg_freqs = mread_get_background(mread);
  options->bg_filename = mread_get_other_bg_src(mread);
  mread_destroy(mread);

  // Check that we got back some motifs
  int num_motif_names = arraylst_size(*motifs);
  if (num_motif_names == 0) {
    die("No motifs could be read.\n");
  }

  // If motifs use protein alphabet we will not scan both strands
  if (!alph_has_complement(options->alphabet)) {
    options->scan_both_strands = FALSE;
  }

  if (options->scan_both_strands == TRUE) {
    // Correct background by averaging on freq. for both strands.
    average_freq_with_complement(options->alphabet, *bg_freqs);
    normalize_subarray(0, alph_size_core(options->alphabet), 0.0, *bg_freqs);
    // Make reverse complement motifs.
    add_reverse_complements(*motifs);
  }

}

/*************************************************************************
 * Write a motif match to the appropriate output files.
 *************************************************************************/
BOOLEAN_T fimo_record_score(
  const FIMO_OPTIONS_T options,
  PATTERN_T *pattern,
  SCANNED_SEQUENCE_T *scanned_seq,
  RESERVOIR_SAMPLER_T *reservoir,
  MATCHED_ELEMENT_T *match
) {

  double pvalue = get_matched_element_pvalue(match);
  if (!options.text_only) {
    // Keep count of how may positions we've scanned.
    add_scanned_sequence_scanned_position(scanned_seq);
    reservoir_sample(reservoir, pvalue);
  }

  BOOLEAN_T score_recorded = FALSE;
  if (pvalue <= options.output_threshold) {
    if (options.text_only) {
      print_site_as_text(match, scanned_seq);
    } else {
      score_recorded = add_pattern_matched_element(pattern, match);
    }
  }

  return score_recorded;
}

/*************************************************************************
 * Calculate the log odds score for a single motif-sized window.
 *************************************************************************/
static inline BOOLEAN_T fimo_score_site(
  const FIMO_OPTIONS_T options,
  char *seq,
  double prior,
  PSSM_T *pssm,
  double *pvalue, // OUT
  double *score // OUT
) {

  ARRAY_T* pv_lookup = pssm->pv;
  MATRIX_T* pssm_matrix = pssm->matrix;
  BOOLEAN_T scorable_site = TRUE;
  double scaled_log_odds = 0.0;

  // For each position in the site
  int motif_position;
  for (motif_position = 0; motif_position < pssm->w; motif_position++) {
    char c = seq[motif_position];
    int aindex = alph_indexc(options.alphabet, c);

    // Check for gaps and ambiguity codes at this site
    if (aindex == -1) {
      scorable_site = FALSE;
      break;
    }
    scaled_log_odds += get_matrix_cell(motif_position, aindex, pssm_matrix);
  }

  if (scorable_site == TRUE) {

    int w = pssm->w;
    *score = get_unscaled_pssm_score(scaled_log_odds, pssm);

    if (!isnan(prior)) {
    // Use the prior to adjust the motif site score.
    // Using the log-odds prior increases the width
      // of the scaled score table by 1.
      ++w;
      double prior_log_odds = (options.alpha) * prior;
      prior_log_odds = my_log2(prior_log_odds / (1.0 - prior_log_odds));
      *score += prior_log_odds;
      scaled_log_odds = raw_to_scaled(*score, w, pssm->scale, pssm->offset);
    }

    // Handle scores that are out of range
    // This should never happen and indicates a bug has been
    // introduced in the code if it does.
    int max_log_odds = get_array_length(pv_lookup) - 1;
    if (scaled_log_odds < 0.0 ) {
      fprintf(stderr, "Scaled log-odds score out of range: %d\n", (int) scaled_log_odds);
      fprintf(stderr, "Assigning 0 to scaled log-odds score.\n");
      scaled_log_odds = 0.0;
    }
    else if ((int) scaled_log_odds > max_log_odds) {
      fprintf(stderr, "Scaled log-odds score out of range: %d\n", (int) scaled_log_odds);
      fprintf(stderr, "Assigning %d to scaled log-odds score.\n", max_log_odds);
      scaled_log_odds = (float) max_log_odds;
    }
    *score = scaled_to_raw(scaled_log_odds, w, pssm->scale, pssm->offset);
    *pvalue = get_array_item((int) scaled_log_odds, pv_lookup);

  }

  return scorable_site;

}

/*************************************************************************
 * Calculate and record the log-odds score and p-value for each 
 * possible motif site in the sequence.
 *
 * Returns the length of the sequence.
 *************************************************************************/
static long fimo_score_sequence(
  const FIMO_OPTIONS_T options,
  RESERVOIR_SAMPLER_T *reservoir,
  DATA_BLOCK_READER_T *fasta_reader,
  DATA_BLOCK_READER_T *psp_reader,
  MOTIF_T* motif,
  MOTIF_T* rev_motif,
  ARRAY_T* bg_freqs,
  PSSM_T*  pssm,
  PSSM_T*  rev_pssm,
  PATTERN_T* pattern
)
{
  assert(motif != NULL);
  assert(bg_freqs != NULL);
  assert(pssm != NULL);

  long num_positions = 0L;

  char strand = (!alph_has_complement(options.alphabet) ? '.' : '+');
  char *fasta_seq_name = NULL;
  BOOLEAN_T fasta_result 
    = fasta_reader->get_seq_name(fasta_reader, &fasta_seq_name);

  if (psp_reader != NULL) {

    // Check current seq in psp_reader
    BOOLEAN_T psp_result = FALSE;
    char *psp_seq_name = NULL;
    psp_result = psp_reader->get_seq_name(psp_reader, &psp_seq_name);
    if (!psp_seq_name || strcmp(fasta_seq_name, psp_seq_name) != 0) {
      // The seq name for the seq reader doesn't match the seq name for the psp reader.
      // Try moving to the next sequence in the psp reader
      psp_result = psp_reader->go_to_next_sequence(psp_reader);
      if (psp_result) {
        // Replace old seq name with new seq name
        myfree(psp_seq_name); 
        psp_result = psp_reader->get_seq_name(psp_reader, &psp_seq_name);
        if (!psp_result) {
          die("Unable to read sequence name from psp reader");
        }
      }

      // Sequences must be in the same order in the FASTA file and PSP file.
      // So the seq names from the psp reader and seq reader had better match
      // now
      if (!psp_result || strcmp(fasta_seq_name, psp_seq_name) != 0) {
        die(
          "Sequence name %s from PSP file %s doesn't match\n"
          "sequence name %s from FASTA file %s.\n",
          psp_seq_name,
          options.psp_filename,
          fasta_seq_name,
          options.seq_filename
        );
      }
    }
    myfree(psp_seq_name);

  }

  // Create a scanned_sequence record and record it in pattern.
  SCANNED_SEQUENCE_T* scanned_seq = 
    allocate_scanned_sequence(fasta_seq_name, fasta_seq_name, pattern);

  // Score and record each possible motif site in the sequence

  // Reserve memory for incoming sequence and prior data
  DATA_BLOCK_T *seq_block = new_sequence_block(pssm->w);
  DATA_BLOCK_T  *prior_block = NULL;

  MATCHED_ELEMENT_T *fwd_match = NULL;
  MATCHED_ELEMENT_T *rev_match = NULL;

  if (options.text_only) {
    // In text only mode we only allocate one match and reuse it
    fwd_match = allocate_matched_element(0, 0, scanned_seq);
    if (rev_pssm) {
      rev_match = allocate_matched_element(0, 0, scanned_seq);
    }
  }

  while (fasta_reader->get_next_block(fasta_reader, seq_block) != FALSE) {

    double prior = NAN;
    ++num_positions;

    int fwd_start = get_start_pos_for_data_block(seq_block);
    int rev_stop = fwd_start;
    int fwd_stop = fwd_start + get_motif_length(motif) - 1;
    int rev_start = fwd_stop;
    if (options.text_only) {
      // Set the values for the current match
      set_matched_element_start(fwd_match, fwd_start);
      set_matched_element_stop(fwd_match, fwd_stop);
      if (rev_pssm) {
        set_matched_element_start(rev_match, rev_start);
        set_matched_element_stop(rev_match, rev_stop);
      }
    }
    else {
      fwd_match = allocate_matched_element(fwd_start, fwd_stop, scanned_seq);
      if (rev_pssm) {
        rev_match = allocate_matched_element(rev_start, rev_stop, scanned_seq);
      }
    }

    // Get sequence data
    char *raw_seq = get_sequence_from_data_block(seq_block);
    if (!options.skip_matched_sequence) {
      set_matched_element_sequence(fwd_match, raw_seq);
    }
    set_matched_element_strand(fwd_match, '+');
    char *rc_seq = NULL;
    if (rev_pssm) {
      // Since we're using the reverse complemment motif
      // convert sequence to reverse complment for output.
      if (!options.skip_matched_sequence) {
        rc_seq = strdup(raw_seq);
        if (!rc_seq) {
          die("Unable to allocate memory for RC sequence string.\n");
        }
        invcomp_seq(options.alphabet, rc_seq, get_motif_length(motif));
        set_matched_element_sequence(rev_match, rc_seq);
      }
      set_matched_element_strand(rev_match, '-');
    }

    // Get corresponding priors.
    if (psp_reader) { 
      get_prior_from_reader(
        psp_reader, 
        fasta_seq_name, 
        fwd_start,
        &prior_block,
        &prior
      ); 
    }

    BOOLEAN_T scoreable_site = FALSE;
    double fwd_score = NAN;
    double fwd_pvalue = NAN;
    double rev_score = NAN;
    double rev_pvalue = NAN;

    // Always score forward strand
    scoreable_site = fimo_score_site(options, raw_seq, prior, pssm, &fwd_pvalue, &fwd_score);
    if (scoreable_site) {
      set_matched_element_score(fwd_match, fwd_score);
      set_matched_element_pvalue(fwd_match, fwd_pvalue);
    }
    if (scoreable_site && rev_pssm != NULL) {
      // Score reverse strand if reverse PSSM available.
      scoreable_site = fimo_score_site(options, raw_seq, prior, rev_pssm, &rev_pvalue, &rev_score);
      if (scoreable_site) {
        set_matched_element_score(rev_match, rev_score);
        set_matched_element_pvalue(rev_match, rev_pvalue);
      }
    }

    BOOLEAN_T fwd_match_recorded = FALSE;
    BOOLEAN_T rev_match_recorded = FALSE;
    if (scoreable_site) {
      if (rev_match && options.max_strand) {
        // Report the larger of the two scores for this site or choose 
        // randomly if tied..
        if (rev_score > fwd_score || (rev_score == fwd_score && drand_mt()>0.5)) {
          rev_match_recorded = 
            fimo_record_score(options, pattern, scanned_seq, reservoir, rev_match);
        }
        else {
          fwd_match_recorded = 
            fimo_record_score(options, pattern, scanned_seq, reservoir, fwd_match);
        }
      }
      else {
        // Certainly record the forward score
        fwd_match_recorded = 
          fimo_record_score(options, pattern, scanned_seq, reservoir, fwd_match);
        if (rev_match) {
          // Record the rev. score too
          rev_match_recorded =
            fimo_record_score(options, pattern, scanned_seq, reservoir, rev_match);
        }
      }
    }


    if (!options.text_only && !fwd_match_recorded) {
      free_matched_element(fwd_match);
    }
    if (rev_pssm != NULL) {
      if (!options.text_only && !rev_match_recorded) {
          free_matched_element(rev_match);
      }
      myfree(rc_seq);
    }
  }

  if (options.text_only) {
    free_matched_element(fwd_match);
    if (rev_match) {
        free_matched_element(rev_match);
    }
  }

  // Count reminaing positions in the sequence.
  num_positions += get_num_read_into_data_block(seq_block);

  free_data_block(seq_block);
  if (prior_block) free_data_block(prior_block);
  myfree(fasta_seq_name);

  return  num_positions;

}

static void fimo_build_pssms(
  MOTIF_T *motif,
  MOTIF_T *rev_motif,
  const FIMO_OPTIONS_T options,
  ARRAY_T *bg_freqs,
  PRIOR_DIST_T *prior_dist,
  PSSM_T **pos_pssm,
  PSSM_T **rev_pssm
) {
  // Build PSSM for motif and tables for p-value calculation.
  // FIXME: the non-averaged freqs should be used for p-values
  *pos_pssm = build_motif_pssm(
     motif, 
     bg_freqs, 
     bg_freqs, 
     prior_dist, 
     options.alpha,
     PSSM_RANGE, 
     0,    // no GC bins
     FALSE // make log-likelihood pssm
  );

  // Open the output file for the p-value lookup table, if requested.
  FILE* pval_lookup_file = NULL;
  if (options.pval_lookup_filename != NULL) {
  }

  if (options.pval_lookup_filename != NULL) {

    FILE *pval_lookup_file = NULL;
    if (open_file(options.pval_lookup_filename, "w", FALSE,
      "p-value lookup table",
      "p-value lookup table", &pval_lookup_file) == 0) {
      exit(EXIT_FAILURE);
    }
    // Print the cumulative lookup table.
    print_array((*pos_pssm)->pv, 8, 6, TRUE, pval_lookup_file);
    // Also print the non-cumulative version.
    int i;
    for (i = 1; i < get_array_length((*pos_pssm)->pv); i++) {
      fprintf(pval_lookup_file, "%g ",
      get_array_item(i-1, (*pos_pssm)->pv) - get_array_item(i, (*pos_pssm)->pv));
    }
    fprintf(pval_lookup_file, "\n");
  }

  // If needed, build the PSSM for the reverse complement motif.
  if (options.scan_both_strands) {
    // FIXME: the non-averaged freqs should be used for p-values
    *rev_pssm = build_motif_pssm(
      rev_motif, 
      bg_freqs, 
      bg_freqs, 
      prior_dist, 
      options.alpha,
      PSSM_RANGE, 
      0, // GC bins
      FALSE
    );
  }
}

/**************************************************************
 * Score each of the sites for each of the selected motifs.
 **************************************************************/
static void fimo_score_each_motif(
  const FIMO_OPTIONS_T options,
  ARRAY_T *bg_freqs,
  ARRAYLST_T *motifs,
  CISML_T *cisml,
  PRIOR_DIST_T *prior_dist,
  DATA_BLOCK_READER_T *fasta_reader,
  DATA_BLOCK_READER_T *psp_reader,
  int *num_scanned_sequences,
  long *num_scanned_positions
) {

  // Create p-value sampling reservoir
  RESERVOIR_SAMPLER_T *reservoir = NULL;
  if (!options.text_only) {
    reservoir = new_reservoir_sampler(10000, NULL);
  }

  // Iterate over all motifs (including rev. comp).
  bool first_motif = true;
  bool need_reset = false;
  int num_motifs = arraylst_size(motifs);
  int num_selected_motifs = get_num_strings(options.selected_motifs);
  int motif_index = 0;
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {

    MOTIF_T* motif = (MOTIF_T *) arraylst_get(motif_index, motifs);
    char* motif_id = get_motif_st_id(motif);
    char* bare_motif_id = get_motif_id(motif);
    char* bare_motif_id2 = get_motif_id2(motif);
    int motif_length = get_motif_length(motif);

    // Is this a selected motif?
    if (num_selected_motifs > 0 
         && have_string(bare_motif_id, options.selected_motifs) == FALSE ) {
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Skipping motif %s.\n", motif_id);
      }
      continue;
    }

    // Only count sequences and positions once.
    if (first_motif) {
      *num_scanned_sequences = 0;
      *num_scanned_positions = 0;
    } else {
      first_motif = false;
    }

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using motif %s of width %d.\n", motif_id, motif_length);
    }

    MOTIF_T* rev_motif = NULL;
    if (options.scan_both_strands) {
      ++motif_index;
      rev_motif = (MOTIF_T *) arraylst_get(motif_index, motifs);
      char * rev_motif_id = get_motif_st_id(rev_motif);
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Using motif %s of width %d.\n", rev_motif_id, motif_length);
      }
    }

    // Only do this if we have more than one motif so users can
    // actually use "-" for sequences with one motif
    if (need_reset) {
      // Reset to start of sequence data. 
      fasta_reader->reset(fasta_reader);
      if (psp_reader) {
         psp_reader->reset(psp_reader);
      }

      if (!options.text_only) {
        clear_reservoir(reservoir);
      }
    }

    PSSM_T* pos_pssm = NULL;
    PSSM_T* rev_pssm = NULL;
    fimo_build_pssms(motif, rev_motif, options, bg_freqs, prior_dist, &pos_pssm, &rev_pssm);

    //if (first_motif) *num_scanned_sequences = 0;

    // Create cisml pattern for this motif.
    PATTERN_T* pattern = allocate_pattern(bare_motif_id, bare_motif_id2);
    if (!options.text_only) {
      add_cisml_pattern(cisml, pattern);
      set_pattern_max_stored_matches(pattern, options.max_stored_scores);
      set_pattern_max_pvalue_retained(pattern, options.output_threshold);
    }

    // Read the FASTA file one sequence at a time.
    while (
      fasta_reader->go_to_next_sequence(fasta_reader) != FALSE
    ) {

      long scanned_positions = fimo_score_sequence(
        options,
        reservoir,
        fasta_reader,
        psp_reader,
        motif,
        rev_motif,
        bg_freqs,
        pos_pssm,
        rev_pssm,
        pattern
      );
      if (first_motif) {
        ++(*num_scanned_sequences);
        *num_scanned_positions += scanned_positions;
      }

    }  // All sequences parsed

    // The pattern is complete.
    if (!options.text_only) {
      set_pattern_is_complete(pattern);
      // Compute q-values, if requested.
      if (options.compute_qvalues) {
        int num_samples = get_reservoir_num_samples_retained(reservoir);
        ARRAY_T *sampled_pvalues = allocate_array(num_samples);
        fill_array(get_reservoir_samples(reservoir), sampled_pvalues);
        pattern_calculate_qvalues(pattern, sampled_pvalues);
        free_array(sampled_pvalues);
      }
    }
    else {
      free_pattern(pattern);
    }

    // Free memory associated with this motif.
    free_pssm(pos_pssm);
    free_pssm(rev_pssm);

    need_reset = true;

  } // All motifs parsed

  if (reservoir != NULL) {
    free_reservoir(reservoir);
  }

}

/*************************************************************************
 * Entry point for fimo
 *************************************************************************/
int main(int argc, char *argv[]) {

  // Get command line arguments
  FIMO_OPTIONS_T options = process_fimo_command_line(argc, argv);

  // Set up motif input
  ARRAYLST_T *motifs = NULL;
  ARRAY_T *bg_freqs = NULL;
  fimo_read_motifs(&options, &motifs, &bg_freqs);

  // Set up sequence input
  DATA_BLOCK_READER_T *fasta_reader = new_seq_reader_from_fasta(
    options.parse_genomic_coord,
    options.alphabet, 
    options.seq_filename
  );

  // Set up priors input
  DATA_BLOCK_READER_T *psp_reader = NULL;
  PRIOR_DIST_T *prior_dist = NULL;
  if (options.psp_filename && options.prior_distribution_filename) {
    prior_dist = new_prior_dist(options.prior_distribution_filename);
    double default_prior = get_prior_dist_median(prior_dist);
    psp_reader =  get_psp_reader(options.psp_filename, default_prior, options.parse_genomic_coord);
  }

  // Create cisml data structure for recording results
  CISML_T *cisml = allocate_cisml("fimo", options.command_line, options.meme_filename, options.seq_filename);
  if (options.threshold_type == PV_THRESH) {
    set_cisml_site_pvalue_cutoff(cisml, options.output_threshold);
  }
  else if (options.threshold_type == QV_THRESH) {
    set_cisml_site_qvalue_cutoff(cisml, options.output_threshold);
  }

  // Initialize tracking variables
  int num_scanned_sequences = 0;
  long num_scanned_positions = 0;

  // Scan all sequences using each motif
  fimo_score_each_motif(
    options, 
    bg_freqs, 
    motifs, 
    cisml,
    prior_dist,
    fasta_reader,
    psp_reader,
    &num_scanned_sequences,
    &num_scanned_positions
  );

  if (options.text_only != TRUE) {
    print_fimo_results(
      cisml, 
      options, 
      bg_freqs,
      motifs,
      num_scanned_sequences,
      num_scanned_positions
    );
  }

  // Clean up.
  free_cisml(cisml);

  fasta_reader->close(fasta_reader);
  free_data_block_reader(fasta_reader);

  if (psp_reader != NULL) {
    psp_reader->close(psp_reader);
    free_data_block_reader(psp_reader);
  }
  if (prior_dist != NULL) {
    free_prior_dist(prior_dist);
  }

  free_motifs(motifs);
  free_array(bg_freqs);
  cleanup_options(options);

  return 0;

}

