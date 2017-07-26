#include "array-list.h"
#include "fimo.h"

/***********************************************************************
  Print plain text record for motif site to standard output.
 ***********************************************************************/
void print_site_as_text(
  MATCHED_ELEMENT_T *match,
  SCANNED_SEQUENCE_T *scanned_seq
);

/***********************************************************************
 * Print FIMO settings information to an XML file
 ***********************************************************************/
void print_settings_xml(FILE *out, FIMO_OPTIONS_T options);

/***********************************************************************
 * Print FIMO specific information to an XML file
 ***********************************************************************/
void print_fimo_xml_file(
  CISML_T *cisml,
  FILE *out,
  FIMO_OPTIONS_T options,
  ARRAYLST_T *motifs,
  ARRAY_T *bgfreq,
  char  *stylesheet,
  int num_seqs,
  long num_residues
);

/**********************************************************************
 * This function saves FIMO results as a tab-delimited text file
 *********************************************************************/
void print_fimo_text_file(FILE *fimo_file, CISML_T *cisml, FIMO_OPTIONS_T options);

/**********************************************************************
 * This function saves FIMO results as a tab-delimited text file
 *********************************************************************/
void print_fimo_gff_file(FILE *fimo_file, CISML_T *cisml, FIMO_OPTIONS_T options);

void fimo_print_version(FILE *fimo_file);

void fimo_print_database_and_motifs(
  FILE *fimo_file,
  FIMO_OPTIONS_T options,
  ARRAYLST_T *motifs,
  ARRAY_T *bgfreq,
  int num_seqs,
  long num_residues
);

void fimo_print_counts(FILE *fimo_file, CISML_T *cisml, FIMO_OPTIONS_T options);

void fimo_print_match_table(
  FILE *fimo_file,
  CISML_T *cisml,
  FIMO_OPTIONS_T options
);

void fimo_print_command_line(FILE *fimo_file, FIMO_OPTIONS_T options);

void fimo_print_parameters(FILE *fimo_file, FIMO_OPTIONS_T options);

/**********************************************************************
 * This function saves FIMO results as an HTML file
 *********************************************************************/
void print_fimo_html_file(
  FILE *fimo_file,
  CISML_T *cisml,
  FIMO_OPTIONS_T options,
  ARRAYLST_T *motifs,
  ARRAY_T *bg_freq,
  int num_seqs,
  long num_residues
);

/**********************************************************************
 * This function saves the FIMO results as a set of files in a
 * directory:
 *
 *   xml_filename will be the name of the FIMO output
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *   gff_filename will be the name of the GFF output
 *
 * allow_clobber will determine whether or not existing files will
 * be overwritten.
 *********************************************************************/
void print_fimo_results(
  CISML_T *cisml, 
  FIMO_OPTIONS_T options,
  ARRAY_T *bg_freqs,
  ARRAYLST_T *motifs,
  int num_seqs,
  long num_residues
);
