#include "motif-in.h"
#include "cisml.h"
#include "config.h"
#include "fimo.h"
#include "fimo-html-string.h"
#include "io.h"

const int MAX_HTML_MATCHES = 1000;

/***********************************************************************
  Print plain text record for motif site to standard output.
 ***********************************************************************/
void print_site_as_text(
  MATCHED_ELEMENT_T *match,
  SCANNED_SEQUENCE_T *scanned_seq
) {

  static BOOLEAN_T print_header = TRUE;

  if (print_header) {
    fprintf(
      stdout, 
      "# motif_id"
      "\tmotif_alt_id"
      "\tsequence_name"
      "\tstart"
      "\tstop"
      "\tstrand"
      "\tscore"
      "\tp-value"
      "\tq-value"
      "\tmatched_sequence\n"
      );
    print_header = FALSE;
  }

  PATTERN_T *pattern = get_scanned_sequence_parent(scanned_seq);
  char *motif_id = get_pattern_accession(pattern);
  char *motif_id2 = get_pattern_name(pattern);
  char *seq_name = get_scanned_sequence_name(scanned_seq);
  char *seq = (char *) get_matched_element_sequence(match);
  int start = get_matched_element_start(match);
  int stop = get_matched_element_stop(match);
  if (stop < start) {
    SWAP(int, start, stop);
  }
  fprintf(
    stdout, 
    "%s\t%s\t%s\t%d\t%d\t%c\t%g\t%.3g\t\t%s\n", 
    motif_id,
    motif_id2,
    seq_name,
    start,
    stop,
    get_matched_element_strand(match),
    get_matched_element_score(match),
    get_matched_element_pvalue(match),
    seq ? seq : ""
  );
}

/***********************************************************************
 * Print FIMO settings information to an XML file
 ***********************************************************************/
void print_settings_xml( FILE *out, FIMO_OPTIONS_T options) {

  fputs("<settings>\n",  out);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "output directory", options.output_dirname);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "MEME file name", options.meme_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "sequence file name", options.seq_filename);
  if (options.bg_filename) {
    fprintf(out, "<setting name=\"%s\">%s</setting>\n", "background file name", options.bg_filename);
  }
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "allow clobber",
    boolean_to_string(options.allow_clobber)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "compute q-values",
    boolean_to_string(options.compute_qvalues)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "parse genomic coord.",
    boolean_to_string(options.parse_genomic_coord)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "text only",
    boolean_to_string(options.text_only)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "scan both strands",
    boolean_to_string(options.scan_both_strands)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%3.2g</setting>\n",
    "output threshold",
    options.output_threshold);
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "threshold type",
    threshold_type_to_string(options.threshold_type)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%d</setting>\n",
    "max stored scores",
    options.max_stored_scores
  );
  fprintf(
    out,
    "<setting name=\"%s\">%3.2g</setting>\n",
    "pseudocount",
    options.pseudocount
  );
  fprintf(
    out,
    "<setting name=\"%s\">%d</setting>\n",
    "verbosity",
    verbosity
  );
  int i = 0;
  int num_strings = get_num_strings(options.selected_motifs);
  for(i = 0; i < num_strings; i++) {
    fprintf(
      out,
      "<setting name=\"%s\">%s</setting>\n", "selected motif",
      get_nth_string(i, options.selected_motifs)
    );
  }

  fputs("</settings>\n",  out);

}

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
) {

  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
  if (stylesheet != NULL) {
    fprintf(out, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n", stylesheet);
  }
  fputs("<!-- Begin document body -->\n", out);
  fputs("<fimo version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n", out);
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("  xmlns:fimo=\"http://noble.gs.washington.edu/schema/fimo\"\n>\n", out);
  fprintf(out, "<command-line>%s</command-line>\n", options.command_line);
  print_settings_xml(out, options);
  fprintf(
    out, 
    "<sequence-data num-sequences=\"%d\" num-residues=\"%ld\" />\n", 
    num_seqs,
    num_residues
  );
  alph_print_xml(options.alphabet, "alphabet", "", "", out);
  int i;
  int num_motifs = arraylst_size(motifs);
  for (i = 0; i < num_motifs; i++) {
    MOTIF_T *motif = arraylst_get(i, motifs);
    char *bare_motif_id = get_motif_id(motif);
    char *bare_motif_id2 = get_motif_id2(motif);
    char *best_possible_match = get_best_possible_match(motif);
    fprintf(
      out,
      "<motif name=\"%s\" alt=\"%s\" width=\"%d\" best-possible-match=\"%s\"/>\n",
      bare_motif_id,
      bare_motif_id2,
      get_motif_length(motif),
      best_possible_match
    );
    if (options.scan_both_strands == TRUE) {
      // Skip RC motif
      ++i;
    }
    myfree(best_possible_match);
  }
  fprintf(
    out,
    "<background source=\"%s\">\n",
    options.bg_filename
  );
  for (i = 0; i < alph_size_core(options.alphabet); i++) {
    fprintf(
      out,
      "<value letter=\"%c\">%1.3f</value>\n",
      alph_char(options.alphabet, i),
      get_array_item(i, bgfreq)
    );
  }
  fputs("</background>\n", out);
  /*
   * FIXME CEG 
   * Use Patterns in cisml to get this information
  int num_incomplete_motifs = get_num_strings(incomplete_motifs);
    if (num_incomplete_motifs > 0) {
    fputs("<incomplete-motifs>\n", out);
    int motif_idx;
    for(motif_idx = 0; motif_idx < num_incomplete_motifs; ++motif_idx) {
      char *name = get_nth_string(motif_idx, incomplete_motifs);
      double max_pvalue = get_nth_score(motif_idx, incomplete_motifs);
      fprintf(out, "<incomplete-motif name=\"%s\" max-pvalue=\"%3.2g\"/>\n", name, max_pvalue);
    }
    fputs("</incomplete-motifs>\n", out);
  }
  */
  fputs("<cisml-file>cisml.xml</cisml-file>\n", out);
  fputs("</fimo>\n", out);
}

/**********************************************************************
 * This function saves FIMO results as a tab-delimited text file
 *********************************************************************/
void print_fimo_text_file(
  FILE *fimo_file, 
  CISML_T *cisml, 
  FIMO_OPTIONS_T options
) {

  char *header = 
    "# motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand"
    "\tscore\tp-value\tq-value\tmatched_sequence\n";
  fputs(header, fimo_file);

  CISML_MATCH_IT_T *it = allocate_cisml_match_iterator(cisml);
  MATCHED_ELEMENT_T *match = NULL;
  while ((match = cisml_match_iterator_next(it)) !=NULL) {
    double pvalue = get_matched_element_pvalue(match);
    double qvalue = get_matched_element_qvalue(match);
    if (options.threshold_type == PV_THRESH) {
      if (pvalue > options.output_threshold) continue;
    }
    else {
      if (qvalue > options.output_threshold) continue;
    }
    int start =  get_matched_element_start(match);
    int stop = get_matched_element_stop(match);
    if (stop < start) {
      SWAP(int, start, stop);
    }
    SCANNED_SEQUENCE_T *scanned_seq = get_matched_element_scanned_seq(match);
    PATTERN_T *pattern = get_scanned_sequence_parent(scanned_seq);
    fprintf(
      fimo_file,
      "%s\t%s\t%s\t%d\t%d\t%c\t%g\t%.3g\t%.3g\t%s\n",
      get_pattern_accession(pattern),
      get_pattern_name(pattern),
      get_scanned_sequence_name(scanned_seq),
      start,
      stop,
      get_matched_element_strand(match),
      get_matched_element_score(match),
      pvalue,
      qvalue,
      get_matched_element_sequence(match)
    );
  }

  free_cisml_match_iterator(it);
}

/**********************************************************************
 * This function saves FIMO results as a tab-delimited text file
 *********************************************************************/
void print_fimo_gff_file(
  FILE *fimo_file,
  CISML_T *cisml,
  FIMO_OPTIONS_T options
) {

  char *header = "##gff-version 3\n";
  fputs(header, fimo_file);

  int num_patterns = get_cisml_num_patterns(cisml);
  PATTERN_T **patterns = get_cisml_patterns(cisml);
  const char *seq_ontology_term = NULL;
  if (alph_extends_dna(options.alphabet) || alph_extends_rna(options.alphabet)) {
    seq_ontology_term = "nucleotide_motif";
  }
  else {
    seq_ontology_term = "sequence_motif";
  }
  int i = 0;
  for (i = 0; i < num_patterns; ++i) {
    int num_seqs = get_pattern_num_scanned_sequences(patterns[i]);
    SCANNED_SEQUENCE_T **seqs = get_pattern_scanned_sequences(patterns[i]);
    int j = 0;
    for (j = 0; j < num_seqs; ++j) {
      int pattern_id = 0;
      int num_matches = get_scanned_sequence_num_matched_elements(seqs[j]);
      MATCHED_ELEMENT_T **matches = get_scanned_sequence_matched_elements(seqs[j]);
      int k = 0;
      for (k = 0; k < num_matches; ++k) {
        ++pattern_id;
        int start =  get_matched_element_start(matches[k]);
        int stop = get_matched_element_stop(matches[k]);
        if (stop < start) {
          SWAP(int, start, stop);
        }
        char *motif_id = get_pattern_accession(patterns[i]);
        char *motif_id2 = get_pattern_name(patterns[i]);
        char qvalue_field[100];
        if (options.compute_qvalues) {
          sprintf(qvalue_field, "qvalue= %.3g;", get_matched_element_qvalue(matches[k]));
        } else {
          qvalue_field[0] = '\0';
        }
        fprintf(
          fimo_file,
          "%s\tfimo\t%s\t%d\t%d\t%3.3g\t%c\t.\t"
          "Name=%s_%s%c;Alias=%s;ID=%s%s%s-%d-%s;pvalue=%.3g;%ssequence=%s;\n",
          get_scanned_sequence_name(seqs[j]),
          seq_ontology_term,
          start,
          stop,
          MIN(1000.0, -4.342 * my_log(get_matched_element_pvalue(matches[k]))), // -10 * log10(x)
          get_matched_element_strand(matches[k]),
          motif_id, get_scanned_sequence_name(seqs[j]), get_matched_element_strand(matches[k]),
          motif_id2,
          motif_id,
          motif_id2 && *motif_id2 ? "-" : "",
          motif_id2 ? motif_id2 : "",
          pattern_id,
          get_scanned_sequence_name(seqs[j]),
          get_matched_element_pvalue(matches[k]),
          qvalue_field,
          get_matched_element_sequence(matches[k])
        );
      }
    }
  }
}

void fimo_print_version(FILE *fimo_file) {
  fprintf(fimo_file, "FIMO version %s, (Release date: %s)", VERSION, ARCHIVE_DATE);
};

void fimo_print_database_and_motifs(
  FILE *fimo_file,
  FIMO_OPTIONS_T options,
  ARRAYLST_T *motifs,
  ARRAY_T *bgfreq,
  int num_seqs,
  long num_residues
) {

  fprintf(
    fimo_file,
    "<p>\n"
    "  DATABASE %s\n"
    "  <br />\n"
    "  Database contains %d sequences, %ld residues\n"
    "</p>\n",
    options.seq_filename,
    num_seqs,
    num_residues
   );
  fprintf(
    fimo_file,
    "<p>\n"
    "  MOTIFS %s (%s)\n"
    "  <table>\n"
    "    <thead>\n"
    "      <tr>\n"
    "        <th style=\"border-bottom: 1px dashed;\">MOTIF</th>\n"
    "        <th style=\"border-bottom: 1px dashed; padding-left: 1em;\">WIDTH</th>\n"
    "        <th style=\"border-bottom: 1px dashed; padding-left: 1em;text-align:left;\" >\n"
    "         BEST POSSIBLE MATCH\n"
    "        </th>\n"
    "      </tr>\n"
    "    </thead>\n"
    "    <tbody>\n",
    options.meme_filename,
    alph_name(options.alphabet)
  );
  int i;
  int num_motifs = arraylst_size(motifs);
  for (i = 0; i < num_motifs; i++) {
    MOTIF_T *motif = arraylst_get(i, motifs);
    char motif_strand = get_motif_strand(motif);
    if (motif_strand == '-') {
      // Motifs on the reverse strand are just
      // rev. complement of the motifs on the forward strand.
      continue;
    }
    char *bare_motif_id = get_motif_id(motif);
    char *best_possible_match = get_best_possible_match(motif);
    fprintf(
      fimo_file,
      "      <tr>\n"
      "        <td style=\"text-align:right;\">%s</td>\n"
      "        <td style=\"text-align:right;padding-left: 1em;\">%d</td>\n"
      "        <td style=\"text-align:left;padding-left: 1em;\">%s</td>\n"
      "       </tr>\n",
      bare_motif_id,
      get_motif_length(motif),
      best_possible_match
    );
    myfree(best_possible_match);
  }
  fprintf(
    fimo_file,
    "    </tbody>\n"
    "  </table>\n"
    "</p>\n"
    "<p>\n"
    "Random model letter frequencies (%s):\n"
    "<br/>\n",
    options.bg_filename
  );
  for (i = 0; i < alph_size_core(options.alphabet); i++) {
    if (i % 9 == 0) {
      fputc('\n', fimo_file);
    }
    fprintf(
      fimo_file,
      "%c %1.3f ",
      alph_char(options.alphabet, i),
      get_array_item(i, bgfreq)
    );
  }
  fputs("</p>", fimo_file);
}

void fimo_print_counts(FILE *fimo_file, CISML_T *cisml, FIMO_OPTIONS_T options) {
  int num_stored_matches = get_cisml_num_stored_matches(cisml);
  fprintf(
    fimo_file,
    "There were %d motif occurences with a %c-value less than %.3g.\n",
    num_stored_matches,
    options.threshold_type == PV_THRESH ? 'p' : 'q',
    options.output_threshold
  );
  if (num_stored_matches >= MAX_HTML_MATCHES) {
    fprintf(
      fimo_file,
      "<b>Only the most significant %d matches are shown here.</b>\n",
      MAX_HTML_MATCHES
    );
  }
}

void fimo_print_match_table(
  FILE *fimo_file,
  CISML_T *cisml,
  FIMO_OPTIONS_T options
) {

  // FIXME Handle output threshold

  // print table header
  char *table_header = NULL;
  if (options.compute_qvalues) {
    table_header =
      "<table border=\"1\">\n"
      "<thead>\n"
      "<tr>\n"
      "<th>Motif ID</th>\n"
      "<th>Alt ID</th>\n"
      "<th>Sequence Name</th>\n"
      "<th>Strand</th>\n"
      "<th>Start</th>\n"
      "<th>End</th>\n"
      "<th>p-value</th>\n"
      "<th>q-value</th>\n"
      "<th>Matched Sequence</th>\n"
      "</tr>\n"
      "</thead>\n"
      "<tbody>\n";
  }
  else {
    table_header =
      "<table border=\"1\">\n"
      "<thead>\n"
      "<tr>\n"
      "<th>Motif</th>\n"
      "<th>Sequence Name</th>\n"
      "<th>Strand</th>\n"
      "<th>Start</th>\n"
      "<th>End</th>\n"
      "<th>p-value</th>\n"
      "<th>Matched Sequence</th>\n"
      "</tr>\n"
      "</thead>\n"
      "<tbody>\n";
  }
  fputs(table_header, fimo_file);

  // print table body
  CISML_MATCH_IT_T *it = allocate_cisml_match_iterator(cisml);
  MATCHED_ELEMENT_T *match = NULL;
  int num_matches = 0;

  while ((match = cisml_match_iterator_next(it)) !=NULL) {
    if (num_matches >= MAX_HTML_MATCHES) {
      break;
    }
    double pvalue = get_matched_element_pvalue(match);
    double qvalue = get_matched_element_qvalue(match);
    if (options.threshold_type == PV_THRESH) {
      if (pvalue > options.output_threshold) continue;
    }
    else {
      if (qvalue > options.output_threshold) continue;
    }
    SCANNED_SEQUENCE_T *scanned_seq = get_matched_element_scanned_seq(match);
    PATTERN_T *pattern = get_scanned_sequence_parent(scanned_seq);
    int start =  get_matched_element_start(match);
    int stop = get_matched_element_stop(match);
    if (stop < start) {
      SWAP(int, start, stop);
    }
    if (options.compute_qvalues) {
      fprintf(
        fimo_file,
        "    <tr>\n"
        "      <td style=\"text-align:left;\">%s</td>\n"
        "      <td style=\"text-align:left;\">%s</td>\n"
        "      <td style=\"text-align:left;\">%s</td>\n"
        "      <td style=\"text-align:center;\">%c</td>\n"
        "      <td style=\"text-align:left;\">%d</td>\n"
        "      <td style=\"text-align:left;\">%d</td>\n"
        "      <td style=\"text-align:left;\">%.3g</td>\n"
        "      <td style=\"text-align:left;\">%.3g</td>\n"
        "      <td style=\"text-align:left;font-size:x-large;font-family:monospace;\">%s</td>\n"
        "   </tr>\n",
        get_pattern_accession(pattern),
        get_pattern_name(pattern),
        get_scanned_sequence_name(scanned_seq),
        get_matched_element_strand(match),
        start,
        stop,
        pvalue,
        qvalue,
        get_matched_element_sequence(match)
      );
    }
    else {
      fprintf(
        fimo_file,
        "    <tr>\n"
        "      <td style=\"text-align:left;\">%s</td>\n"
        "      <td style=\"text-align:left;\">%s</td>\n"
        "      <td style=\"text-align:center;\">%c</td>\n"
        "      <td style=\"text-align:left;\">%d</td>\n"
        "      <td style=\"text-align:left;\">%d</td>\n"
        "      <td style=\"text-align:left;\">%.3g</td>\n"
        "      <td style=\"text-align:left;font-size: x-large;font-family:monospace;\">%s</td>\n"
        "   </tr>\n",
        get_pattern_name(pattern),
        get_scanned_sequence_name(scanned_seq),
        get_matched_element_strand(match),
        start,
        stop,
        pvalue,
        get_matched_element_sequence(match)
      );
    }
    ++num_matches;
  }

  // print table close
  fputs("</tbody>\n</table>\n", fimo_file);

  free_cisml_match_iterator(it);
}

void fimo_print_command_line(FILE *fimo_file, FIMO_OPTIONS_T options) {
  fputs(options.command_line, fimo_file);
}

void fimo_print_parameters(FILE *fimo_file, FIMO_OPTIONS_T options) {

  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">output_directory = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">MEME file name = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">sequence file name = %s</td>\n"
    "  </tr>",
    options.output_dirname,
    options.meme_filename,
    options.seq_filename
  );
  if (options.psp_filename != NULL) {
    fprintf(
      fimo_file,
      "  <tr>\n"
      "    <td style=\"padding-right: 2em\">PSP filename = %s</td>\n"
      "    <td style=\"padding-left: 5em; padding-right: 2em\">prior dist. filename = %s</td>\n"
      "    <td style=\"padding-left: 5em; padding-right: 2em\"></td>\n"
      "  </tr>",
      options.psp_filename,
      options.prior_distribution_filename
    );
  }
  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">background file name = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">alphabet = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">max stored scores = %d</td>\n"
    "  </tr>",
    options.bg_filename,
    alph_name(options.alphabet),
    options.max_stored_scores
  );
  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">allow clobber = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">compute q-values = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">parse genomic coord. = %s</td>\n"
    "  </tr>\n",
    options.allow_clobber ? "true" : "false",
    options.compute_qvalues ? "true" : "false",
    options.parse_genomic_coord ? "true" : "false"
  );
  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">text only = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">scan both strands = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">max strand = %s</td>\n"
    "  </tr>\n",
    options.text_only ? "true" : "false",
    options.scan_both_strands ? "true" : "false",
    options.max_strand ? "true" : "false"
  );
  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">threshold type = %s</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">output theshold = %g</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">pseudocount = %g</td>\n"
    "  </tr>\n",
    options.threshold_type == PV_THRESH ? "p-value" : "q-value",
    options.output_threshold,
    options.pseudocount
  );
  fprintf(
    fimo_file,
    "  <tr>\n"
    "    <td style=\"padding-right: 2em\">alpha = %g</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\">verbosity = %d</td>\n"
    "    <td style=\"padding-left: 5em; padding-right: 2em\"></td>\n"
    "  </tr>\n",
    options.alpha,
    verbosity
  );

};

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
) {
  const int MAX_TAG_SIZE = 1000;
  int html_string_size = strlen(fimo_html_string);
  int i = 0;
  for (i = 0; i < html_string_size; ++i) {
    if (fimo_html_string[i] != '@') {
      fputc(fimo_html_string[i], fimo_file);
    }
    else {
      char buffer[MAX_TAG_SIZE];
      ++i;
      int j = 0;
      while (fimo_html_string[i] != '@' && j < (MAX_TAG_SIZE - 1)) {
        buffer[j] = fimo_html_string[i];
        ++j;
        ++i;
      }
      if (fimo_html_string[i] != '@') {
        die("FIMO tag buffer length exceeded\n");
      }
      buffer[j] = '\0';
      if (strcmp("version", buffer) == 0) {
        fimo_print_version(fimo_file);
      }
      else if (strcmp("database_and_motifs", buffer) == 0) {
        fimo_print_database_and_motifs(fimo_file, options, motifs, bg_freq, num_seqs, num_residues);
      }
      else if (strcmp("counts", buffer) == 0) {
        fimo_print_counts(fimo_file, cisml, options);
      }
      else if (strcmp("match_table", buffer) == 0) {
        fimo_print_match_table(fimo_file, cisml, options);
      }
      else if (strcmp("command_line", buffer) == 0) {
        fimo_print_command_line(fimo_file, options);
      }
      else if (strcmp("parameters", buffer) == 0) {
        fimo_print_parameters(fimo_file, options);
      }
    }
  }

}

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
) {

  const BOOLEAN_T PRINT_WARNINGS = FALSE;
  if (create_output_directory(
       options.output_dirname,
       options.allow_clobber,
       PRINT_WARNINGS
      )
    ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }

  // Print CisML.
  FILE *cisml_file = fopen(options.cisml_path, "w");
  if (!cisml_file) {
    die("Couldn't open file %s for output.\n", options.cisml_path);
  }
  print_cisml(cisml_file, cisml, TRUE, NULL, TRUE);
  fclose(cisml_file);

  // Print XML.
  FILE *fimo_file = fopen(options.xml_path, "w");
  if (!fimo_file) {
    die("Couldn't open file %s for output.\n", options.xml_path);
  }
  print_fimo_xml_file(
    cisml,
    fimo_file,
    options,
    motifs,
    bg_freqs,
    NULL,
    num_seqs,
    num_residues
  );
  fclose(fimo_file);

  // Print plain text.
  fimo_file = fopen(options.text_path, "w");
  if (!fimo_file) {
    die("Couldn't open file %s for output.\n", options.text_path);
  }
  print_fimo_text_file(fimo_file, cisml, options);
  fclose(fimo_file);

  // Print GFF.
  fimo_file = fopen(options.gff_path, "w");
  if (!fimo_file) {
    die("Couldn't open file %s for output.\n", options.gff_path);
  }
  print_fimo_gff_file(fimo_file, cisml, options);
  fclose(fimo_file);

  // Print HTML
  fimo_file = fopen(options.html_path, "w");
  if (!fimo_file) {
    die("Couldn't open file %s for output.\n", options.html_path);
  }
  print_fimo_html_file(
    fimo_file, 
    cisml, 
    options, 
    motifs,
    bg_freqs,
    num_seqs, 
    num_residues
  );
  fclose(fimo_file);
}

