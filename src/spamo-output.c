#define _GNU_SOURCE
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <sys/wait.h>

#include "config.h"
#include "spamo-output.h"
#include "spamo-matches.h"
#include "utils.h"
#include "seq-reader-from-fasta.h"

/**************************************************************************
 * Outputs txt for the spacings
 **************************************************************************/
void output_secondary_motif_txt(FILE *txt_output, 
    MOTIF_DB_T *primary_db, MOTIF_T *primary_motif, 
    SECONDARY_MOTIF_T *parent, SECONDARY_MOTIF_T *smotif, 
    int n_secondary_motifs, LINKLST_T *rmotifs) {
  LINK_T *node;
  SIGSPACE_T sig;
  int i, quad;
  char* orient_names[NORIENTS];
  // this maintains the previous definition of orientation but actual names would be better
  orient_names[LEFT | SAME] = "0";
  orient_names[LEFT | OPPO] = "1";
  orient_names[RIGHT | SAME] = "2";
  orient_names[RIGHT | OPPO] = "3";
  orient_names[UP_SEC_PAL] = "4";
  orient_names[UP_PRI_PAL] = "5";
  orient_names[DOWN_PRI_PAL] = "6";
  orient_names[DOWN_SEC_PAL] = "7";
  orient_names[BOTH_PAL] = "8";

  for (i = 0; i < smotif->sig_count; ++i) {
    sig = smotif->sigs[i];
    // Output txt line.
    fprintf(txt_output, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.2e\t%d\t%s\t%d\t%d\t%.2e\t%.2e\n",
      primary_db->name,
      get_motif_id(primary_motif),
      strcmp(get_motif_id2(primary_motif), "") ? get_motif_id2(primary_motif) : ".",
      get_motif_consensus(primary_motif),
      smotif->db->name,
      get_motif_id(smotif->motif),
      strcmp(get_motif_id2(smotif->motif), "") ? get_motif_id2(smotif->motif) : ".",
      get_motif_consensus(smotif->motif),
      get_motif_trim_left(smotif->motif),
      get_motif_trim_right(smotif->motif),
      parent ? smotif->db->name : ".",
      parent ? get_motif_id(parent->motif) : ".",
      parent && strcmp(get_motif_id2(parent->motif), "") ? get_motif_id2(parent->motif) : ".",
      smotif->min_pvalue * n_secondary_motifs, // E-value
      sig.bin, // gap
      orient_names[sig.orient], // orientation
      smotif->spacings[sig.orient].count[sig.bin], // count
      smotif->total_spacings, // total occurrences of secondary in margins
      sig.pvalue < 1e-10 ? sig.pvalue/NORIENTS : 1-pow((1-sig.pvalue),1.0/NORIENTS),
      sig.pvalue
    );
  }

  if (rmotifs != NULL && linklst_size(rmotifs) > 0) {
    for (node = linklst_first(rmotifs); node != NULL; node = linklst_next(node)) {
      SECONDARY_MOTIF_T *rmotif = linklst_get(node);
      output_secondary_motif_txt(txt_output, primary_db, primary_motif, smotif, rmotif, n_secondary_motifs, NULL);
    }
  }
}


/**************************************************************************
 * Dump sequence matches sorted by the name of the sequence.
 *
 * Outputs Columns:
 *   1) Trimmed lowercase sequence with uppercase matches.
 *   2) Position of the secondary match within the whole sequence.
 *   3) Sequence fragment that the primary matched.
 *   4) Strand of the primary match (+|-)
 *   5) Sequence fragment that the secondary matched.
 *   6) Strand of the secondary match (+|-)
 *   7) Is the primary match on the same strand as the secondary (s|o)
 *   8) Is the secondary match downstream or upstream (d|u)
 *   9) The gap between the primary and secondary matches
 *  10) The name of the sequence
 *  11) The p-value of the bin containing the match (adjusted for # of bins)
 *  ---if the FASTA input file sequence names are in Genome Browser format:
 *  12-14) Position of primary match in BED coordinates
 *  15) Position of primary match in Genome Browser coordinates
 *  16-18) Position of secondary match in BED coordinates
 *  19) Position of secondary match in Genome Browser coordinates
 *
 * If you wish to sort based on the gap column:
 * Sort individual output:
 *  sort -n -k 9,9 -o seqs_primary_secondary.txt seqs_primary_secondary.txt
 * Or sort all outputs:
 *  for f in seqs_*.txt; do sort -n -k 9,9 -o $f $f; done
 * Or to get just locations of primary motif in BED coordinates
 * where the secondary is on the opposite strand, upstream with a gap of 118bp:
 *   awk '$7=="o" && $8=="u" && $9==118 {print $12"\t"$13"\t"$14;}' seqs_primary_secondary.txt 
 *
 **************************************************************************/
static void dump_sequence_matches(FILE *out, int margin, int bin, 
    double sigthresh, BOOLEAN_T sig_only, RBTREE_T *sequences,
    MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif,
    ARRAY_T **matches) {
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int idx, seqlen, i, j, start, end, secondary, secondary_pos, primary_len, secondary_len, distance;
  BOOLEAN_T primary_rc, secondary_rc, downstream; 
  char *buffer, *seq, *primary_match, *secondary_match;
  ARRAY_T *secondary_array;
  ALPH_T *alph;
  // get the alphabet
  alph = get_motif_alph(primary_motif);
  // allocate a buffer for copying the trimmed sequence into and modify it
  seqlen = margin * 2 + get_motif_trimmed_length(primary_motif);
  buffer = (char*)mm_malloc(sizeof(char) * (seqlen + 1));
  // get the lengths of the motifs
  primary_len = get_motif_trimmed_length(primary_motif);
  secondary_len = get_motif_trimmed_length(secondary_motif->motif); 
  // allocate some strings for storing the matches
  primary_match = (char*)mm_malloc(sizeof(char) * (primary_len + 1));
  secondary_match = (char*)mm_malloc(sizeof(char) * (secondary_len + 1));
  // add null byte at the end of the match strings
  primary_match[primary_len] = '\0';
  secondary_match[secondary_len] = '\0';

  // iterate over all the sequences
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    primary_rc = get_array_item(0, sequence->primary_matches) < 0;

    //secondary = matches[sequence->index];
    secondary_array = matches[sequence->index];
    if (! secondary_array) continue;
    int n_secondary_matches = get_array_length(secondary_array);
    for (idx=0; idx<n_secondary_matches; idx++) {
      secondary = get_array_item(idx, secondary_array);
      secondary_rc = secondary < 0;
      secondary_pos = abs(secondary);

      // calculate the distance
      if (secondary_pos <= margin) {
        distance = margin - secondary_pos - secondary_len + 1;
        downstream = primary_rc;
      } else {
        distance = secondary_pos - margin - primary_len - 1;
        downstream = !primary_rc;
      }

      // copy the trimmed sequence
      seq = sequence->data;
      for (i = 0; i < seqlen; ++i) {
        buffer[i] = (alph_is_case_insensitive(alph) ? tolower(seq[i]) : seq[i]);
      }
      buffer[seqlen] = '\0';

      // uppercase primary
      start = margin;
      end = margin + primary_len;
      for (i = start, j = 0; i < end; ++i, ++j) {
        buffer[i] = (alph_is_case_insensitive(alph) ? toupper(buffer[i]) : buffer[i]);
        primary_match[j] = buffer[i];
      }

      // uppercase secondary
      // note orign was one, subtract 1 to make origin zero as required for arrays
      start = secondary_pos -1;
      end = start + secondary_len;
      for (i = start, j = 0; i < end; ++i, ++j) {
        buffer[i] = (alph_is_case_insensitive(alph) ? toupper(buffer[i]) : buffer[i]);
        secondary_match[j] = buffer[i];
      }

      // get the p-value of the seconndary match
      SPACING_T *spacings;
      if (secondary_rc == primary_rc) {
        spacings = downstream ? secondary_motif->spacings+(SAME+RIGHT) : secondary_motif->spacings+(SAME+LEFT); 
      } else {
        spacings = downstream ? secondary_motif->spacings+(OPPO+RIGHT) : secondary_motif->spacings+(OPPO+LEFT); 
      }
      double p_value = spacings->pvalue[distance/bin];

      // skip match if not significant and only reporting significant matches
      if (sig_only && (p_value > sigthresh)) continue;

      // output line to file
      fprintf(out, "%s    %3d    %s    %s    %s    %s    %s    %s    %3d    %s    %.1e", 
          buffer, 
          secondary_pos, 
          primary_match, 
          (primary_rc ? "-" : "+"), 
          secondary_match, 
          (secondary_rc ? "-" : "+"), 
          (secondary_rc == primary_rc ? "s" : "o"),
          (downstream ? "d" : "u"), 
          distance, 
          sequence->name,
          p_value
      );

      // Parse the sequence name to see if we can get genomic coordinates
      // and print additional columns with primary and secondary matches
      // in both BED and Genome Browser coordinates.
      char *chr_name;
      size_t chr_name_len;
      int start_pos, end_pos;
      if (parse_genomic_coordinates_helper(
          sequence->name,
          &chr_name,
          &chr_name_len,
          &start_pos,
          &end_pos))
      {
        // Get the start and end of the primary match in 
        // 0-relative, half-open genomic coordinates.
        int p_start = start_pos + fabs(get_array_item(0, sequence->primary_matches)) - 1;
        int p_end = p_start + primary_len;
        // Get the start and end of the secondary match in 
        // 0-relative, half-open genomic coordinates.
        int s_start, s_end;
        if ( (!primary_rc && downstream) || (primary_rc && !downstream) ) {
          s_start = p_end + distance;
          s_end = s_start + secondary_len;
        } else {
          s_end = p_start - distance;
          s_start = s_end - secondary_len;
        }
        fprintf(out, "    %s    %d    %d    %s:%d-%d", 
          chr_name, p_start, p_end, chr_name, p_start+1, p_end);
        fprintf(out, "    %s    %d    %d    %s:%d-%d\n", 
          chr_name, s_start, s_end, chr_name, s_start+1, s_end);
      } else {
        fprintf(out, "\n");
      }

    } // secondary match
  } // primary match

  free(buffer);
  free(primary_match);
  free(secondary_match);
}

/**************************************************************************
 * Makes the histogram file name from the details in the motifs and the
 * file extension.
 * Caller is responsible for freeing memory
 **************************************************************************/
static char* make_pattern_file_name(char *prefix, char *ext, MOTIF_T* primary, SECONDARY_MOTIF_T *secondary) {
  const char *fmt = "%s_%s_db%d_%s.%s";
  char dummy[1];
  char *ret;
  int len;
  len = snprintf(dummy, 1, fmt, prefix, get_motif_id(primary), secondary->db->id, get_motif_id(secondary->motif), ext);
  ret = mm_malloc(sizeof(char) * (len+1));
  snprintf(ret, (len+1), fmt, prefix, get_motif_id(primary), secondary->db->id, get_motif_id(secondary->motif), ext);
  return ret;
}

/**************************************************************************
 * Create an output file and dump the sequence matches to file.
 **************************************************************************/
void output_sequence_matches(char *dir, int margin, int bin, double sigthresh,
    BOOLEAN_T sig_only, RBTREE_T *sequences, MOTIF_T *primary_motif,
    SECONDARY_MOTIF_T *secondary_motif, ARRAY_T **matches) {
  FILE *out;
  int file_name_len;
  char *file_path, *file_name;
  file_name = make_pattern_file_name("seqs", "txt", primary_motif, secondary_motif);
  file_path = make_path_to_file(dir, file_name);
  out = fopen(file_path, "w");
  dump_sequence_matches(out, margin, bin, sigthresh, sig_only, sequences, primary_motif, secondary_motif, matches);
  fclose(out);
  free(file_path);
  free(file_name);
}

