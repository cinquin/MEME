#include <assert.h>
#include <math.h>
#include <string.h>
#include "alphabet.h"
#include "config.h"
#include "html-monolith.h"
#include "json-writer.h"
#include "mcast-match.h"
#include "object-list.h"
#include "projrel.h"
#include "utils.h"
#include "cisml.h"

struct motif_hit_t {
  char *motif_id;
  int motif_index;
  char *seq;
  size_t start;
  size_t stop;
  char strand;
  double pvalue;
};

struct mcast_match_t {
  int cluster_id;
  OBJECT_LIST_T *motif_hits;
  char *seq_name;
  size_t seq_length; // length of the input sequence
  size_t seq_start; // starting position of the sequence (normally 0 unless parse-genomic-coords used).
  char *sequence; // the smallest sequence segment that includes all the motif hits
  char *lflank; // a small piece of sequence on the left flank of the matching sequence
  char *rflank; // a small piece of sequence on the right flank of the matching sequence
  int gc_bin;
  double score;
  double gc;
  double evalue;
  double pvalue;
  double qvalue;
  size_t start;
  size_t stop;
};

/***********************************************************************
 * allocate_mcast_match
 *
 * This function generates a new mcast match object.
 ***********************************************************************/
MCAST_MATCH_T *allocate_mcast_match() {
  static int cluster_id = 0;
  MCAST_MATCH_T *mcast_match = mm_malloc(sizeof(MCAST_MATCH_T));
  ++cluster_id;
  mcast_match->cluster_id = cluster_id;

  mcast_match->motif_hits
    = new_object_list(NULL, NULL, NULL, free_motif_hit);
  mcast_match->seq_name = NULL;
  mcast_match->seq_length = 0;
  mcast_match->seq_start = 0;
  mcast_match->sequence = NULL;
  mcast_match->lflank = NULL;
  mcast_match->rflank = NULL;
  mcast_match->gc_bin = -1;
  mcast_match->score = -1.0;
  mcast_match->gc = NAN;
  mcast_match->evalue = NAN;
  mcast_match->pvalue = NAN;
  mcast_match->qvalue = NAN;
  mcast_match->start = -1;
  mcast_match->stop = -1;

  return mcast_match;
}

/***********************************************************************
 * copy_mcast_match
 *
 * This function deep copies an mcast match object.
 *
 ***********************************************************************/
void *copy_mcast_match(void *m) {
  MCAST_MATCH_T *mcast_match = m;
  MCAST_MATCH_T *new_mcast_match = mm_malloc(sizeof(MCAST_MATCH_T));

  new_mcast_match->motif_hits
    = new_object_list(NULL, NULL, NULL, free_motif_hit);
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while (motif_hit != NULL) {
    // FIXME actually do the copying
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }

  new_mcast_match->seq_name = strdup(mcast_match->seq_name);
  new_mcast_match->seq_length = mcast_match->seq_length;
  new_mcast_match->seq_start = mcast_match->seq_start;
  new_mcast_match->sequence = strdup(mcast_match->sequence);
  new_mcast_match->lflank = strdup(mcast_match->lflank);
  new_mcast_match->rflank = strdup(mcast_match->rflank);
  new_mcast_match->gc_bin = mcast_match->gc_bin;
  new_mcast_match->score = mcast_match->score;
  new_mcast_match->gc = mcast_match->gc;
  new_mcast_match->evalue = mcast_match->evalue;
  new_mcast_match->pvalue = mcast_match->pvalue;
  new_mcast_match->qvalue = mcast_match->qvalue;
  new_mcast_match->start = mcast_match->start;
  new_mcast_match->stop = mcast_match->stop;

  return new_mcast_match;
}

/***********************************************************************
 * free_mcast_match
 *
 * This function frees the memory associated with a mcast_matc object.
 ***********************************************************************/
void free_mcast_match(void* m) {
  MCAST_MATCH_T *mcast_match = m;
  free_object_list(mcast_match->motif_hits);
  myfree(mcast_match->seq_name);
  myfree(mcast_match->sequence);
  myfree(mcast_match->lflank);
  myfree(mcast_match->rflank);
  myfree(mcast_match);
}

/***********************************************************************
 * compare_mcast_matches
 *
 * This function compares two mcast match objects.
 * mcast match objects are ordered by pvalue from lower to higher.
 * This function returns 
 *    1 if m1->pvalue > m2->pvalue,
 *    0 if m1->pvalue = m2->pvalue,
 *    -1 if m1->pvalue < m2->pvalue,
 *
 ***********************************************************************/
int compare_mcast_matches(void *m1, void *m2) {
  MCAST_MATCH_T *match1 = m1;
  MCAST_MATCH_T *match2 = m2;

  if (match1->pvalue > match2->pvalue) {
    return -1;
  }
  else if (match1->pvalue == match2->pvalue) {
    return 0;
  }
  else {
    return 1;
  }
}

/***********************************************************************
 * compare_mcast_match_pvalues
 *
 * This function compares the p-values of two MCAST matches
 * and returns -1, 0, or 1 as the p-value of the 1st argument is
 * less than, equal to, or greater than the p-value of the 2nd argument.
 *
 ***********************************************************************/
int compare_mcast_match_pvalues(const void *param1, const void *param2) {
  MCAST_MATCH_T *match1 = *((MCAST_MATCH_T **) param1);
  MCAST_MATCH_T *match2 = *((MCAST_MATCH_T **) param2);

  if (match1->pvalue > match2->pvalue) {
    return 1;
  }
  else if (match1->pvalue < match2->pvalue) {
    return -1;
  }
  else {
    return 0;
  }
}

/***********************************************************************
 * rev_compare_mcast_match_pvalues
 *
 * This function compares the p-values of two MCAST matches
 * and returns 1, 0, or -1 as the p-value of the 1st argument is
 * less than, equal to, or greater than the p-value of the 2nd argument.
 *
 ***********************************************************************/
int rev_compare_mcast_match_pvalues(const void *param1, const void *param2) {
  MCAST_MATCH_T *match1 = *((MCAST_MATCH_T **) param1);
  MCAST_MATCH_T *match2 = *((MCAST_MATCH_T **) param2);

  if (match1->pvalue > match2->pvalue) {
    return -1;
  }
  else if (match1->pvalue < match2->pvalue) {
    return 1;
  }
  else {
    return 0;
  }
}

/***********************************************************************
 * get_mcast_match_cluster_id
 *
 * This function returns the cluster id from an MCAST match object.
 *
 ***********************************************************************/
int get_mcast_match_cluster_id(MCAST_MATCH_T *match) {
  return match->cluster_id;
}

/***********************************************************************
 * set_mcast_match_seq_name
 *
 * This function sets the seq_name for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_seq_name(MCAST_MATCH_T *match, char *name) {
  match->seq_name = strdup(name);
}

/***********************************************************************
 * get_mcast_match_seq_name
 *
 * This function returns the seq_name from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_seq_name(MCAST_MATCH_T *match) {
  return match->seq_name;
}

/***********************************************************************
 * set_mcast_match_seq_length
 *
 * This function sets the seq_length for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_seq_length(MCAST_MATCH_T *match, size_t length) {
  match->seq_length = length;
}

/***********************************************************************
 * get_mcast_match_seq_length
 *
 * This function returns the seq_length from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_seq_length(MCAST_MATCH_T *match) {
  return match->seq_length;
}

/***********************************************************************
 * set_mcast_match_seq_start
 *
 * This function sets the seq_start for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_seq_start(MCAST_MATCH_T *match, size_t start) {
  match->seq_start = start;
}

/***********************************************************************
 * get_mcast_match_seq_start
 *
 * This function returns the seq_start from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_seq_start(MCAST_MATCH_T *match) {
  return match->seq_start;
}

/***********************************************************************
 * set_mcast_match_sequence
 *
 * This function sets the sequence for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_sequence(MCAST_MATCH_T *match, char *sequence) {
  match->sequence = strdup(sequence);
}

/***********************************************************************
 * get_mcast_match_sequence
 *
 * This function returns the sequence from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_sequence(MCAST_MATCH_T *match) {
  return match->sequence;
}

/***********************************************************************
 * set_mcast_match_lflank
 *
 * This function sets the lflank for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_lflank(MCAST_MATCH_T *match, char *lflank) {
  match->lflank = strdup(lflank);
}

/***********************************************************************
 * get_mcast_match_lflank
 *
 * This function returns the lflank from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_lflank(MCAST_MATCH_T *match) {
  return match->lflank;
}

/***********************************************************************
 * set_mcast_match_rflank
 *
 * This function sets the rflank for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_rflank(MCAST_MATCH_T *match, char *rflank) {
  match->rflank = strdup(rflank);
}

/***********************************************************************
 * get_mcast_match_rflank
 *
 * This function returns the rflank from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_rflank(MCAST_MATCH_T *match) {
  return match->rflank;
}

/***********************************************************************
 * set_mcast_match_gc_bin
 *
 * This function sets the gc_bin for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_gc_bin(MCAST_MATCH_T *match, int gc_bin) {
  match->gc_bin = gc_bin;
}

/***********************************************************************
 * get_mcast_match_gc_bin
 *
 * This function gets the gc_bin for an MCAST match object.
 *
 ***********************************************************************/
int get_mcast_match_gc_bin(MCAST_MATCH_T *match) {
  return match->gc_bin;
}

/***********************************************************************
 * set_mcast_match_score
 *
 * This function sets the score for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_score(MCAST_MATCH_T *match, double score) {
  match->score = score;
}

/***********************************************************************
 * get_mcast_match_score
 *
 * This function returns the score from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_score(MCAST_MATCH_T *match) {
  return match->score;
}

/***********************************************************************
 * set_mcast_match_gc
 *
 * This function sets the gc value for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_gc(MCAST_MATCH_T *match, double gc) {
  match->gc = gc;
}

/***********************************************************************
 * get_mcast_match_gc
 *
 * This function returns the gc value from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_gc(MCAST_MATCH_T *match) {
  return match->gc;
}

/***********************************************************************
 * set_mcast_match_evalue
 *
 * This function sets the evalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_evalue(MCAST_MATCH_T *match, double evalue) {
  match->evalue = evalue;
}

/***********************************************************************
 * get_mcast_match_evalue
 *
 * This function returns the evalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_evalue(MCAST_MATCH_T *match) {
  return match->evalue;
}

/***********************************************************************
 * set_mcast_match_pvalue
 *
 * This function sets the pvalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_pvalue(MCAST_MATCH_T *match, double pvalue) {
  match->pvalue = pvalue;
}

/***********************************************************************
 * get_mcast_match_pvalue
 *
 * This function returns the pvalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_pvalue(MCAST_MATCH_T *match) {
  return match->pvalue;
}

/***********************************************************************
 * set_mcast_match_qvalue
 *
 * This function sets the qvalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_qvalue(MCAST_MATCH_T *match, double qvalue) {
  match->qvalue = qvalue;
}

/***********************************************************************
 * get_mcast_match_qvalue
 *
 * This function returns the qvalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_qvalue(MCAST_MATCH_T *match) {
  return match->qvalue;
}

/***********************************************************************
 * set_mcast_match_start
 *
 * This function sets the match start position for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_start(MCAST_MATCH_T *match, size_t start) {
  match->start = start;
}

/***********************************************************************
 * get_mcast_match_start
 *
 * This function returns the match start position from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_start(MCAST_MATCH_T *match) {
  return match->start;
}

/***********************************************************************
 * set_mcast_match_stop
 *
 * This function sets the match stop position for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_stop(MCAST_MATCH_T *match, size_t stop) {
  match->stop = stop;
}

/***********************************************************************
 * get_mcast_match_stop
 *
 * This function returns the match start position from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_stop(MCAST_MATCH_T *match) {
  return match->stop;
}

/***********************************************************************
 * allocate_motif_hit
 *
 * This function generates a new motif hit object.
 ***********************************************************************/
MOTIF_HIT_T *allocate_motif_hit(
  char *motif_id, 
  int motif_index,
  char *seq,
  char strand, 
  size_t start, 
  size_t stop, 
  double pvalue
) {
  MOTIF_HIT_T *motif_hit = mm_malloc(sizeof(MOTIF_HIT_T));

  motif_hit->motif_id = strdup(motif_id);
  motif_hit->motif_index = motif_index;
  motif_hit->seq = strdup(seq);
  motif_hit->strand = strand;
  motif_hit->start = start;
  motif_hit->stop = stop;
  motif_hit->pvalue = pvalue;

  return motif_hit;
}

/***********************************************************************
 * free_motif_hit
 *
 * This function frees the memory associated with a motif hit object.
 ***********************************************************************/
void free_motif_hit(MOTIF_HIT_T* motif_hit) {
  myfree(motif_hit->motif_id);
  myfree(motif_hit->seq);
  myfree(motif_hit);
}

/***********************************************************************
 * get_motif_hit_motif_id
 *
 * This function returns a pointer to a string containing the motif id
 * associated with a motif hit object. 
 *
 * The caller is NOT responsible for freeing the returnted string.
 * It will be freed when the motif hit is freed.
 *
 ***********************************************************************/
const char* get_motif_hit_motif_id(MOTIF_HIT_T* motif_hit) {
  return motif_hit->motif_id;
}

/***********************************************************************
 * get_motif_hit_seq
 *
 * This function returns a pointer to a string containing the motif id
 * associated with a motif hit object. 
 *
 * The caller is NOT responsible for freeing the returnted string.
 * It will be freed when the motif hit is freed.
 *
 ***********************************************************************/
const char* get_motif_hit_seq(MOTIF_HIT_T* motif_hit) {
  return motif_hit->seq;
}

/***********************************************************************
 * get_motif_hit_strand
 *
 * This function returns the strand assocaited with a motif hit object.
 *
 ***********************************************************************/
char get_motif_hit_strand(MOTIF_HIT_T* motif_hit) {
  return motif_hit->strand;
}

/***********************************************************************
 * get_motif_hit_start
 *
 * This function returns the starting position assocaited with 
 * a motif hit object.
 *
 ***********************************************************************/
size_t get_motif_hit_start(MOTIF_HIT_T* motif_hit) {
  return motif_hit->start;
}

/***********************************************************************
 * get_motif_hit_stop
 *
 * This function returns the stop position assocaited with 
 * a motif hit object.
 *
 ***********************************************************************/
size_t get_motif_hit_stop(MOTIF_HIT_T* motif_hit) {
  return motif_hit->stop;
}

/***********************************************************************
 * get_motif_hit_pvalue
 *
 * This function returns the pvalue assocaited with a motif hit object.
 *
 ***********************************************************************/
double get_motif_hit_pvalue(MOTIF_HIT_T* motif_hit) {
  return motif_hit->pvalue;
}

/***********************************************************************
 * add_mcast_match_motif_hit
 *
 * This function adds a motif_hit to an mcast_match.
 * Freeing the mcast_match will free the motif_hits that have been added.
 ***********************************************************************/
void add_mcast_match_motif_hit(
    MCAST_MATCH_T* mcast_match,
    MOTIF_HIT_T *motif_hit
) {
  store_object(motif_hit, NULL, 0.0, mcast_match->motif_hits);
}

/***********************************************************************
 * print_mcast_match
 *
 * This function prints the data from mcast_match object.
 *
 * This is used for debugging.
 *
 ***********************************************************************/
void print_mcast_match(FILE *output, void *m) {
  MCAST_MATCH_T *mcast_match = m;
  fprintf(
    output,
    "%s\t%5.1g\t%3.1g\t%3.1g\t%zd\t%zd\n",
    mcast_match->seq_name,
    mcast_match->score,
    mcast_match->pvalue,
    mcast_match->evalue,
    mcast_match->start + 1, // convert to 1-indexed
    mcast_match->stop + 1 // convert to 1-indexed
  );
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while(motif_hit != NULL) {
    fprintf(
      output,
      "\t%s\t%3.1g\t%c\t%zd\t%zd\n",
      motif_hit->motif_id,
      motif_hit->pvalue,
      motif_hit->strand,
      motif_hit->start + 1, // convert to 1-indexed
      motif_hit->stop + 1 // convert to 1-indexed
    );
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }
}
/**********************************************************************
  print_mcast_match_as_cisml

  Print a heap of mcast_matches asCisML XML
**********************************************************************/
void print_mcast_match_as_cisml(
  FILE* out, 
  MCAST_MATCH_T *mcast_match
) {

  assert(out != NULL);
  assert(mcast_match != NULL);

  int cluster_id = get_mcast_match_cluster_id(mcast_match);
  char *seq_name = get_mcast_match_seq_name(mcast_match);
  char *match_seq = get_mcast_match_sequence(mcast_match);
  double match_score = get_mcast_match_score(mcast_match);
  double match_evalue = get_mcast_match_evalue(mcast_match);
  double match_pvalue = get_mcast_match_pvalue(mcast_match);
  double match_qvalue = get_mcast_match_qvalue(mcast_match);
  size_t match_start = get_mcast_match_start(mcast_match);
  size_t match_stop = get_mcast_match_stop(mcast_match);

  fprintf(out, "<multi-pattern-scan");
  fprintf(out, " score=\"%g\"", match_score);
  fprintf(out, " pvalue=\"%.5g\"", match_pvalue);
  fprintf(out, ">\n");
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while(motif_hit != NULL) {
    const char *motif_id = get_motif_hit_motif_id(motif_hit);
    const char *hit_seq = get_motif_hit_seq(motif_hit);
    char strand = get_motif_hit_strand(motif_hit);
    double hit_pvalue = get_motif_hit_pvalue(motif_hit);
    size_t hit_start = get_motif_hit_start(motif_hit);
    size_t hit_stop = get_motif_hit_stop(motif_hit);
    if (strand == '-') {
      // Reverse strand, swap start and stop
      size_t tmp = hit_stop;
      hit_stop = hit_start;
      hit_start = tmp;
    }
    fprintf(
      out,
      "<pattern accession=\"%s\" name=\"%s\">\n",
      motif_id,
      motif_id
    );
    fprintf(
      out,
      "<scanned-sequence accession=\"%s\" name=\"%s\">\n",
      seq_name,
      seq_name
    );
    fprintf(
      out,
      "<matched-element start=\"%zd\" stop=\"%zd\" pvalue=\"%.5g\">\n",
      hit_start + 1, // convert to 1-indexed
      hit_stop + 1, // convert to 1-indexed
      hit_pvalue
    );
    fprintf(out, "<sequence>%s</sequence>\n", hit_seq);
    fputs("</matched-element>\n", out);
    fputs("</scanned-sequence>\n", out);
    fputs("</pattern>\n", out);
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }
  fprintf(
    out, 
    "<mem:match cluster-id=\"cluster-%d\" "
    "seq-name=\"%s\" start=\"%zd\" stop=\"%zd\" "
    "evalue=\"%.5g\" qvalue=\"%.5g\""
    ">",
    cluster_id,
    seq_name,
    match_start + 1, // convert to 1-indexed
    match_stop + 1, // convert to 1-indexed
    match_evalue,
    match_qvalue
  );
  fprintf(out, "%s\n", match_seq);
  fputs("</mem:match>\n", out);
  fputs("</multi-pattern-scan>\n", out);
}

/**********************************************************************
  mcast_print_results_as_cisml

  Print a heap of mcast_matches as CisML XML
**********************************************************************/
void mcast_print_results_as_cisml(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {

  assert(options != NULL);
  assert(mcast_matches != NULL);

  // Open the CisML XML file for output
  FILE *out = fopen(options->cisml_path, "w");
  if (!out) {
    die("Couldn't open file %s for output.\n", options->cisml_path);
  }

  print_cisml_start(out, options->program, TRUE, NULL, TRUE);

  fprintf(out, "<parameters>\n");
  fprintf(
    out,
    "<pattern-file>%s</pattern-file>\n",
    options->motif_filename 
  );
  fprintf(
    out,
    "<sequence-file>%s</sequence-file>\n",
    options->seq_filename
  );
  if (options->bg_filename != NULL) {
    fprintf(
      out,
      "<background-seq-file>%s</background-seq-file>\n",
      options->bg_filename
    );
  }
  fprintf(
    out,
    "<pattern-pvalue-cutoff>%g</pattern-pvalue-cutoff>\n",
    options->motif_pthresh
  );
  fprintf(
    out,
    "<sequence-pvalue-cutoff>%g</sequence-pvalue-cutoff>\n",
    options->p_thresh
  );
  fprintf(out, "</parameters>\n" );

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      print_mcast_match_as_cisml(out, match);
    }
  }

  print_cisml_end(out);

  fclose(out);
}

/**********************************************************************
  mcast_print_results_as_gff

  Print an array of mcast_matches as gff.
**********************************************************************/
void mcast_print_results_as_gff(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {
  assert(options != NULL);
  assert(mcast_matches != NULL);

  FILE *out = fopen(options->gff_path, "w");
  if (!out) {
    die("Couldn't open file %s for output.\n", options->gff_path);
  }

  // Print header line
  char *header = "##gff-version 3\n";
  fprintf(out, "%s", header);

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      int id = get_mcast_match_cluster_id(match);
      char *seq_name = get_mcast_match_seq_name(match);
      double score = MIN(1000.0, -4.342 * my_log(pvalue)); // -10 * log10(pvalue)
      size_t start = get_mcast_match_start(match);
      size_t stop = get_mcast_match_stop(match);
      fprintf(
        out,
        "%s\tMCAST\ttranscriptional_cis_regulatory_region\t%zd\t%zd\t%3.3g\t.\t.\t"
        "ID=cluster-%d;pvalue=%.3g;evalue=%3.3g;qvalue=%.3g\n",
        seq_name,
        start + 1, // convert to 1-indexed
        stop + 1, // convert to 1-indexed
        score,
        id,
        pvalue,
        evalue,
        qvalue
      );
    }
  }

  fclose(out);
}

static void mcast_print_match_as_json(
  JSONWR_T *json,
  MCAST_MATCH_T *match
) {
  MOTIF_HIT_T *hit;
  size_t lflank_len, rflank_len, match_len;
  char* match_with_flanks;
  // concat the match with flanks
  lflank_len = strlen(match->lflank);
  rflank_len = strlen(match->rflank);
  match_len = strlen(match->sequence);
  match_with_flanks = mm_malloc(sizeof(char) * (lflank_len + match_len + rflank_len + 1));
  strncpy(match_with_flanks, match->lflank, lflank_len);
  strncpy(match_with_flanks+lflank_len, match->sequence, match_len);
  strncpy(match_with_flanks+(lflank_len + match_len), match->rflank, rflank_len);
  match_with_flanks[lflank_len + match_len + rflank_len] = '\0';
  // output
  jsonwr_start_object_value(json);
  jsonwr_lng_prop(json, "db", 0); // currently only one DB
  jsonwr_str_prop(json, "name", match->seq_name);
  jsonwr_lng_prop(json, "start", match->seq_start);
  jsonwr_lng_prop(json, "length", match->seq_length);
  jsonwr_dbl_prop(json, "score", match->score);
  // as MCAST allows NAN values which can't be represented in JSON we check that the number is not NaN first
  if (!isnan(match->pvalue)) jsonwr_dbl_prop(json, "pvalue", match->pvalue);
  if (!isnan(match->evalue)) jsonwr_dbl_prop(json, "evalue", match->evalue);
  if (!isnan(match->qvalue)) jsonwr_dbl_prop(json, "qvalue", match->qvalue);
  // output segments
  jsonwr_property(json, "segs");
  jsonwr_start_array_value(json);
  jsonwr_start_object_value(json); //start seg - currently MCAST only outputs one segment
  jsonwr_lng_prop(json, "pos", match->start - match->seq_start - lflank_len); // zero-indexed relative to scanned sequence
  jsonwr_str_prop(json, "data", match_with_flanks); // a sequence segment

  jsonwr_property(json, "hits");
  jsonwr_start_array_value(json);
  // WARNING, object-list uses internal state to do the iteration and there is
  // no way to reset the damn thing so if this starts misbehaving check that
  // the internal state isn't being set elsewhere... Sigh.
  while ((hit = retrieve_next_object(match->motif_hits)) != NULL) {
    jsonwr_start_object_value(json); // start hit
    jsonwr_lng_prop(json, "idx", hit->motif_index - 1); // convert to zero-indexed
    jsonwr_lng_prop(json, "pos", hit->start - match->seq_start); // zero-indexed relative to scanned sequence
    jsonwr_bool_prop(json, "rc", hit->strand == '-'); // strand can be '+', '-' or '.'
    // as MCAST allows NAN values which can't be represented in JSON we check that the number is not NaN first
    if (!isnan(hit->pvalue)) jsonwr_dbl_prop(json, "pvalue", hit->pvalue);
    jsonwr_end_object_value(json); // end hit
  }
  jsonwr_end_array_value(json); // end hits

  jsonwr_end_object_value(json); // end seg
  jsonwr_end_array_value(json); // end segs

  jsonwr_end_object_value(json); // end sequence
  free(match_with_flanks);
}

/**********************************************************************
  mcast_print_results_as_html

  Print an array of mcast_matches as HTML.
**********************************************************************/
#define TOP_RESULTS 100
void mcast_print_results_as_html(
  int argc,
  char **argv,
  time_t *start,
  double duration,
  ARRAY_T *bgfreqs,
  char *psp_file,
  char *prior_dist_file,
  int num_motifs,
  MOTIF_T *motifs,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  int num_seqs, 
  long num_residues,
  long max_stored_scores,
  double motif_pthresh,
  double alpha,
  bool synth,
  uint32_t seed,
  bool genome_coords,
  bool stats_available,
  MHMMSCAN_OPTIONS_T *options
) {
  HTMLWR_T *html;
  JSONWR_T *json;
  MOTIF_T *motif;
  int i, j, k;
  bool norc;
  assert(options != NULL);
  assert(mcast_matches != NULL);
  // Check if the alphabet has a complement and the motifs were created with both strands
  // if both conditions are true then MAST searched with the reverse complement, otherwise norc is true.
  // TODO maybe add a norc option?
  norc = !(alph_has_complement(options->alphabet) && get_motif_strands(motif_at(motifs, 0)));
  // setup html monolith writer
  json = NULL;
  if ((html = htmlwr_create(get_meme_etc_dir(), "mcast_template.html", true))) {
    htmlwr_set_dest_name(html, options->output_dirname, "mcast.html");
    htmlwr_replace(html, "mcast_data.js", "data");
    htmlwr_replace(html, "mcast_data_aux.js", "data_aux");
    json = htmlwr_output(html);
    if (json == NULL || strcmp("data", htmlwr_replacing(html)) != 0) die("Template does not contain data section.\n");
  } else {
    if (verbosity >= QUIET_VERBOSE) fprintf(stderr, "Failed to open html template file.\n");
    return;
  }
  // now output some json
  // output some top level variables
  jsonwr_str_prop(json, "program", "MCAST");
  jsonwr_str_prop(json, "version", VERSION);
  jsonwr_str_prop(json, "release", ARCHIVE_DATE);
  jsonwr_str_prop(json, "revision", REVISION);
  
  // output command line
  jsonwr_str_array_prop(json, "cmd", argv, argc);

  // output settings
  jsonwr_property(json, "settings");
  jsonwr_start_object_value(json);
  // general settings
  jsonwr_lng_prop(json, "max_stored_scores", max_stored_scores);
  jsonwr_dbl_prop(json, "motif_p_thresh", motif_pthresh);
  switch (options->output_thresh_type) {
    case PVALUE:
      jsonwr_dbl_prop(json, "p_thresh", options->p_thresh);
      break;
    case EVALUE:
      jsonwr_dbl_prop(json, "e_thresh", options->e_thresh);
      break;
    case QVALUE:
      jsonwr_dbl_prop(json, "q_thresh", options->q_thresh);
      break;
  }
  jsonwr_dbl_prop(json, "alpha", alpha);
  jsonwr_dbl_prop(json, "min_match_score", options->dp_thresh);
  jsonwr_lng_prop(json, "max_gap", options->max_gap);
  jsonwr_dbl_prop(json, "cost_factor", options->egcost);
  jsonwr_dbl_prop(json, "gap_open_cost", options->gap_open);
  jsonwr_dbl_prop(json, "gap_extend_cost", options->gap_extend);
  jsonwr_bool_prop(json, "hard_mask", options->hard_mask);
  jsonwr_bool_prop(json, "norc", norc);
  jsonwr_bool_prop(json, "synth", synth);
  jsonwr_bool_prop(json, "genomic_coord", genome_coords);
  jsonwr_bool_prop(json, "stats", stats_available);
  jsonwr_lng_prop(json, "seed", seed);
  // end of settings
  jsonwr_end_object_value(json);

  // alphabet
  jsonwr_property(json, "alphabet");
  alph_print_json(options->alphabet, json);

  // determine background source
  char *bg_source;
  bool bg_file = false;
  if (options->bg_filename == NULL || strcmp(options->bg_filename, "--nrdb--") == 0) {
    bg_source = "nrdb";
  } else if (strcmp(options->bg_filename, "--uniform--") == 0) {
    bg_source = "uniform";
  } else if (strcmp(options->bg_filename, "--motif--") == 0) {
    bg_source = "motif";
  } else {
    bg_source = "file";
    bg_file = true;
  }
  // background
  jsonwr_property(json, "background");
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "source", bg_source);
  if (bg_file) jsonwr_str_prop(json, "file", options->bg_filename);
  jsonwr_property(json, "freqs");
  jsonwr_start_array_value(json);
  for (i = 0; i < alph_size_core(options->alphabet); i++) {
    jsonwr_dbl_value(json, get_array_item(i, bgfreqs));
  }
  jsonwr_end_array_value(json);
  jsonwr_end_object_value(json);

  // motif dbs
  jsonwr_property(json, "motif_dbs");
  jsonwr_start_array_value(json);
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "file", options->motif_filename);
  jsonwr_end_object_value(json);
  jsonwr_end_array_value(json);
  
  // sequence dbs (currently there is only one but I use an array for future expansion)
  jsonwr_property(json, "sequence_dbs");
  jsonwr_start_array_value(json);
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "file", options->seq_filename);
  if (psp_file != NULL && prior_dist_file != NULL) {
    jsonwr_str_prop(json, "psp_file", psp_file);
    jsonwr_str_prop(json, "dist_file", prior_dist_file);
  }
  jsonwr_lng_prop(json, "sequence_count", num_seqs);
  jsonwr_lng_prop(json, "residue_count", num_residues);
  jsonwr_end_object_value(json);
  jsonwr_end_array_value(json);

  // motifs
  jsonwr_property(json, "motifs");
  jsonwr_start_array_value(json);
  for (i = 0; i < num_motifs; i += (norc ? 1 : 2)) {
    motif = motif_at(motifs, i);
    jsonwr_start_object_value(json);
    jsonwr_lng_prop(json, "db", 0);
    jsonwr_str_prop(json, "id", get_motif_id(motif));
    jsonwr_str_prop(json, "alt", get_motif_id2(motif));
    jsonwr_lng_prop(json, "len", get_motif_length(motif));
    jsonwr_dbl_prop(json, "nsites", get_motif_nsites(motif));
    jsonwr_log10num_prop(json, "evalue", get_motif_log_evalue(motif), 1);
    jsonwr_str_prop(json, "url", get_motif_url(motif));
    jsonwr_property(json, "pwm");
    jsonwr_start_array_value(json);
    for (j = 0; j < get_motif_length(motif); j++) {
      jsonwr_start_array_value(json);
      for (k = 0; k < get_motif_alph_size(motif); k++) {
        jsonwr_dbl_value(json, get_matrix_cell(j, k, get_motif_freqs(motif)));
      }
      jsonwr_end_array_value(json);
    }
    jsonwr_end_array_value(json);
    jsonwr_end_object_value(json);
  }
  jsonwr_end_array_value(json);

  // do some calculations needed by the HTML
  // this enables us to split the data section in the HTML into 2.
  long max_cluster_len = 0;
  long max_display_num = 0;
  for (i = 0; i < num_matches; i++) {
    MCAST_MATCH_T *match;
    match = mcast_matches[i];
    // skip any matches that don't pass the threshold
    if (!stats_available || (
          match->evalue <= options->e_thresh &&
          match->pvalue <= options->p_thresh &&
          match->qvalue <= options->q_thresh)) {
      size_t lflank_len, rflank_len, cluster_len, display_num;
      lflank_len = strlen(match->lflank);
      rflank_len = strlen(match->rflank);
      cluster_len = lflank_len + (match->stop - match->start + 1) + rflank_len;
      if (max_cluster_len < cluster_len) max_cluster_len = cluster_len;
      display_num = match->stop + rflank_len;
      if (max_display_num < display_num) max_display_num = display_num;
    }
  }
  jsonwr_property(json, "calc");
  jsonwr_start_object_value(json);
  // used to scale the clusters
  jsonwr_lng_prop(json, "max_cluster_len", max_cluster_len);
  // used to provide accurate edge padding for the needles
  jsonwr_lng_prop(json, "max_display_num", max_display_num);
  jsonwr_end_object_value(json); // end calc


  // write out the top matches
  int top = (num_matches < TOP_RESULTS ? num_matches : TOP_RESULTS);
  jsonwr_property(json, "matches");
  jsonwr_start_array_value(json);
  for (i = 0; i < top; i++) {
    MCAST_MATCH_T *match;
    match = mcast_matches[i];
    // skip any matches that don't pass the threshold
    if (!stats_available || (
          match->evalue <= options->e_thresh &&
          match->pvalue <= options->p_thresh &&
          match->qvalue <= options->q_thresh)) {
      mcast_print_match_as_json(json, match);
    }
  }
  jsonwr_end_array_value(json); // end matches
  
  // output runtime
  jsonwr_property(json, "runtime");
  jsonwr_start_object_value(json);
  jsonwr_str_prop(json, "host", hostname());
  jsonwr_time_prop(json, "when", start);
  jsonwr_dbl_prop(json, "seconds", duration);
  jsonwr_end_object_value(json); // runtime

  // now write out aux matches
  json = htmlwr_output(html);
  if (json == NULL || strcmp("data_aux", htmlwr_replacing(html)) != 0) die("Template does not contain aux data section.\n");
  jsonwr_property(json, "matches");
  jsonwr_start_array_value(json);
  for (i = top; i < num_matches; i++) {
    MCAST_MATCH_T *match;
    match = mcast_matches[i];
    // skip any matches that don't pass the threshold
    if (!stats_available || (
          match->evalue <= options->e_thresh &&
          match->pvalue <= options->p_thresh &&
          match->qvalue <= options->q_thresh)) {
      mcast_print_match_as_json(json, match);
    }
  }
  jsonwr_end_array_value(json); // end matches

  // finish writing html file
  if (htmlwr_output(html) != NULL) {
    die("Found another JSON replacement!\n");
  }
  htmlwr_destroy(html);
}

/**********************************************************************
  mcast_print_results_as_text

  Print a heap of mcast_matches as plain text.
**********************************************************************/
void mcast_print_results_as_text(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {

  assert(options != NULL);
  assert(mcast_matches != NULL);

  FILE *out = NULL;
  if (options->text_only) {
    out = stdout;
  }
  else {
    out = fopen(options->text_path, "w");
    if (!out) {
      die("Couldn't open file %s for output.\n", options->text_path);
    }
  }

  // Print header line
  fprintf(
    out,
    "# pattern_name\t"
    "sequence_name\t"
    "start\t"
    "stop\t"
    "score\t"
    "p-value\t"
    "E-value\t"
    "q-value\t"
    "matched_sequence\n"
  );

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      int id = get_mcast_match_cluster_id(match);
      char *seq_name = get_mcast_match_seq_name(match);
      char *seq = get_mcast_match_sequence(match);
      // translate to primary symbols in matching sequence
      // this has the effect of uppercasing the standard alphabets
      char *seq_prime = strdup(seq);
      translate_seq(options->alphabet, seq_prime, ALPH_NO_UNKNOWN | ALPH_NO_ALIASES);
      double score = get_mcast_match_score(match);
      size_t start = get_mcast_match_start(match);
      size_t stop = get_mcast_match_stop(match);
      fprintf(
        out,
        "cluster-%d\t%s\t%zd\t%zd\t%g\t%.5g\t%.5g\t%.5g\t%s\n",
        id,
        seq_name,
        start + 1, // convert to 1-indexed
        stop + 1, // convert to 1-indexed
        score,
        pvalue,
        evalue,
        qvalue,
        seq_prime
      );
      free(seq_prime);
    }
  }

  if (out != stdout) {
    fclose(out);
  }
}

