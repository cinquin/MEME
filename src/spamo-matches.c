
#include <errno.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "binomial.h"
#include "spamo-matches.h"
#include "linked-list.h"
#include "seq.h"
#include "alignment.h"
#include "motif.h"
#include "motif-spec.h"

#define min(a,b)      ((a)<(b))?(a):(b)
#define max(a,b)      ((a)>(b))?(a):(b)

const uint8_t orient2set[] = {
  Q_TOP_LEFT, // upstream same strand
  Q_TOP_RIGHT, // downstream same strand
  Q_BOTTOM_LEFT, // upstream opposite strand
  Q_BOTTOM_RIGHT,  // downstream opposite strand
  Q_TOP_LEFT | Q_BOTTOM_LEFT, // upstream with a palindromic secondary motif
  Q_TOP_LEFT | Q_BOTTOM_RIGHT, // upstream with a palindromic primary motif
  Q_TOP_RIGHT | Q_BOTTOM_LEFT, // downstream with a palindromic primary motif
  Q_TOP_RIGHT | Q_BOTTOM_RIGHT, // downstream with a palindromic secondary motif
  Q_TOP_LEFT | Q_TOP_RIGHT | Q_BOTTOM_LEFT | Q_BOTTOM_RIGHT // both motifs are palindromic
};
// Initialize logo motif names.
static char* LOGO_NAME[NQUADS] = {
  " Upstream, Same Strand",
  " Upstream, Opposite Strand",
  " Downstream, Same Strand",
  " Downstream, Opposite Strand",
};

/**************************************************************************
 * Compares two keys of a secondary motif for equality.
 * First compares based on the motif database and second compares
 * based on the name of the motif.
 **************************************************************************/
int secondary_key_compare(const void *p1, const void *p2) {
  const SECONDARY_KEY_T *k1, *k2;
  k1 = (SECONDARY_KEY_T*)p1;
  k2 = (SECONDARY_KEY_T*)p2;
  if (k1->db_id == k2->db_id) {
    return strcmp(k1->motif_id, k2->motif_id);
  } else if (k1->db_id < k2->db_id) {
    return -1;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Copies a key of a secondary motif
 **************************************************************************/
void* secondary_key_copy(void *p) {
  SECONDARY_KEY_T *source, *dest;
  source = (SECONDARY_KEY_T*)p;
  dest = mm_malloc(sizeof(SECONDARY_KEY_T));
  dest->db_id = source->db_id;
  dest->motif_id = source->motif_id;
  return dest;
}

/**************************************************************************
 * Creates a structure to hold the sequence database information.
 * All strings are copied.
 **************************************************************************/
SEQUENCE_DB_T* create_sequence_db(char *file) {
  SEQUENCE_DB_T* db;
  struct stat stbuf;
  char *c;
  int len, i, stat_result;
  //allocate memory
  db = (SEQUENCE_DB_T*)mm_malloc(sizeof(SEQUENCE_DB_T));
  //copy the source
  copy_string(&(db->source), file);
  //copy the name, assumes unix style file separator
  copy_string(&(db->name), ((c = strrchr(file, '/')) ? c+1 : file));
  //clean up the name
  len = strlen(db->name);
  //if it ends with .fasta or .fa then truncate the string
  if (len > 6 && strcmp(db->name+(len - 6), ".fasta") == 0) db->name[len-6] = '\0';
  else if (len > 3 && strcmp(db->name+(len - 3), ".fa") == 0) db->name[len-3] = '\0';
  //replace underscores; removed by TLB
  //for (c = db->name; *c != '\0'; c++) if (*c == '_') *c = ' ';
  //stat the file and get the last modified date
  if (stat(file, &stbuf) < 0) {
    //an error occured
    die("Failed to stat file \"%s\", error given as: %s\n", file, strerror(errno));
  }
  db->last_mod = stbuf.st_mtime;
  //initilize defaults
  db->loaded = 0;
  db->excluded_nomatch = 0;
  db->excluded_similar = 0;
  return db;
}

/**************************************************************************
 * Destroy the sequence db
 **************************************************************************/
void destroy_sequence_db(SEQUENCE_DB_T *db) {
  free(db->source);
  free(db->name);
  free(db);
}

/**************************************************************************
 * Creates a structure to hold the motif database information.
 * All strings are copied.
 **************************************************************************/
MOTIF_DB_T* create_motif_db(int id, char *file) {
  MOTIF_DB_T* db;
  struct stat stbuf;
  char *c;
  int len, i, stat_result;
  //allocate memory
  db = (MOTIF_DB_T*)mm_malloc(sizeof(MOTIF_DB_T));
  //copy the id
  db->id = id;
  //copy the source
  copy_string(&(db->source), file);
  //copy the name, assumes unix style file separator
  copy_string(&(db->name), ((c = strrchr(file, '/')) ? c+1 : file));
  //clean up the name
  len = strlen(db->name);
  //if it ends with .meme then truncate the string
  if (len > 5 && strcmp(db->name+(len - 5), ".meme") == 0) db->name[len-5] = '\0';
  //replace underscores; Removed by TLB.
  //for (c = db->name; *c != '\0'; c++) if (*c == '_') *c = ' ';
  //stat the file and get the last modified date
  if (stat(file, &stbuf) < 0) {
    //an error occured
    die("Failed to stat file \"%s\", error given as: %s\n", file, strerror(errno));
  }
  db->last_mod = stbuf.st_mtime;
  //initilize defaults
  db->loaded = 0;
  db->excluded = 0;
  return db;
}

/**************************************************************************
 * Destroy the motif db.
 **************************************************************************/
void destroy_motif_db(void *secondary_db) {
  MOTIF_DB_T *db;

  db = (MOTIF_DB_T*)secondary_db;
  free(db->source);
  free(db->name);
  free(db);
}

/**************************************************************************
 * Destroy the sequence (including the name attribute)
 * Note that the data attribute is destroyed seperately
 **************************************************************************/
void destroy_sequence(void *sequence) {
  SEQUENCE_T *seq;

  seq = (SEQUENCE_T*)sequence;
  free(seq->name);
  free(seq->primary_matches);
  if (seq->trimmed_primary_matches) free(seq->trimmed_primary_matches);
  free(seq);
}

/**************************************************************************
 * Used by create_secondary_motif to simplify the initilization of the
 * spacings array for each of the NQUAD spacing types.
 **************************************************************************/
static void init_spacings(SPACING_T *frequencies, int bin_count) {
  int bin;
  frequencies->bin_count = bin_count;
  frequencies->count = (int*)mm_malloc(sizeof(int) * bin_count);
  frequencies->pvalue = (double*)mm_malloc(sizeof(double) * bin_count);
  frequencies->sequences = (LINKLST_T**)mm_malloc(sizeof(LINKLST_T*) * bin_count);
  for (bin = 0; bin < bin_count; ++bin) {
    frequencies->count[bin] = 0;
    frequencies->pvalue[bin] = 1;
    frequencies->sequences[bin] = linklst_create();
  }
}

/**************************************************************************
 * Used by destroy_secondary_motif to deallocate memory for the count,
 * pvalue and sequences fields of the spacing type.
 **************************************************************************/
static void destroy_spacings(SPACING_T *frequencies) {
  int bin;
  free(frequencies->count);
  free(frequencies->pvalue);
  for (bin = 0; bin < frequencies->bin_count; ++bin) {
    linklst_destroy(frequencies->sequences[bin]);
  }
}

/**************************************************************************
 * Create a secondary motif. As the number of sequences is unknown at this
 * point the sequence_matches array is left unallocated. All pvalues are
 * initilized to 1.
 **************************************************************************/
SECONDARY_MOTIF_T* create_secondary_motif(int margin, int bin, 
    MOTIF_DB_T *db, MOTIF_T *motif) {
  int bin_count, i, norients;
  SECONDARY_MOTIF_T *smotif;
  bool revcomp;
  // check if reverse complement motifs are supported
  revcomp = alph_has_complement(get_motif_alph(motif));
  norients = (revcomp ? NORIENTS : 2);
  // create the secondary motif
  smotif = mm_malloc(sizeof(SECONDARY_MOTIF_T));
  smotif->idx = -1; // not assigned until written out
  smotif->db = db;
  smotif->motif = motif;
  //calculate the number of bins needed for this motif
  bin_count = (int)((margin - get_motif_trimmed_length(motif) + 1) / bin) + 1;
  //allocate spacings
  for (i = 0; i < norients; ++i) init_spacings((smotif->spacings)+i, bin_count);
  for (; i < NORIENTS; ++i) memset((smotif->spacings)+i, 0, sizeof(SPACING_T));
  smotif->total_spacings = 0;
  smotif->max_in_one_bin = 0;
  //these will be allocated after we've filled the spacings tables
  //and calculated the most significant spacings
  smotif->sigs = NULL;
  smotif->sig_count = 0;
  smotif->min_pvalue = DBL_MAX;
  smotif->best_orient = 0;
  smotif->best_bin = 0;
  smotif->seqs = NULL;
  smotif->seq_count = 0;
  return smotif;
}

/**************************************************************************
 * Destroys a secondary motif. 
 * It takes a void * pointer so it can be used in the collection objects.
 **************************************************************************/
void destroy_secondary_motif(void *p) {
  int i;
  SECONDARY_MOTIF_T *smotif = (SECONDARY_MOTIF_T*)p;
  free_motif(smotif->motif);
  free(smotif->motif);
  for (i = 0; i < NORIENTS; ++i) destroy_spacings((smotif->spacings)+i);
  if (smotif->sigs) free(smotif->sigs);
  if (smotif->seqs) free(smotif->seqs);
  free(smotif);
}

/**************************************************************************
 * Create a group with a good secondary motif and a group
 * of redundant secondary motifs.
 **************************************************************************/
GROUPED_MOTIF_T* create_grouped_motif(SECONDARY_MOTIF_T* best) {
  GROUPED_MOTIF_T *group;
  group = (GROUPED_MOTIF_T*)mm_malloc(sizeof(GROUPED_MOTIF_T));
  group->best = best;
  group->others = linklst_create();
  return group;
}

/**************************************************************************
 * Destroys a grouped motif
 * 
 * Relies on the secondary motifs begin destroyed elsewhere.
 **************************************************************************/
void destroy_grouped_motif(void *p) {
  GROUPED_MOTIF_T *group = (GROUPED_MOTIF_T*)p;
  linklst_destroy(group->others);
  free(group);
}

/**************************************************************************
 * Copy a (segment of a) string, optionally with reverse complementation.
 *
 * The destination string must be allocated and freed by the caller
 * and must be at least size "length+1".
 * Copying begins at the address "from_str" and the segment of 
 * length bytes is copied to "to_str", with reverse complementation
 * if requested.  Then a null is added.
 **************************************************************************/
void copy_string_with_rc(
  ALPH_T *alph,
  char *from_str, // Start of segment to copy.
  char *to_str, // Destination for copy.
  int length, // The number of characters to copy.
  bool rc // Copy the reverse complement of the segment.  
) 
{
  // copy sequence segment, reverse-complementing if necessary
  if (rc) {
    int i;
    char *p;
    p = from_str+(length - 1);
    for (i = 0; i < length; i++, p--) to_str[i] = comp_sym(alph, *p);
  } else {
    strncpy(to_str, from_str, length);
  }
  to_str[length] = '\0';
} // copy_string_with_rc

/**************************************************************************
 * Puts counts into the spacing bins.
 **************************************************************************/
void bin_matches(int margin, int bin_size, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, ARRAY_T **matches) {
  int i, bin, primary_len, secondary_len, secondary, secondary_pos, primary_rc, secondary_rc, quad, distance, max_distance;
  RBNODE_T *node;
  SECONDARY_MOTIF_T *smotif;
  SEQUENCE_T *sequence;
  SPACING_T *spacing;
  ARRAY_T *secondary_array;

  primary_len = get_motif_trimmed_length(primary_motif);

  smotif = secondary_motif;
  secondary_len = get_motif_trimmed_length(smotif->motif);

  // Note that distance counts from zero
  max_distance = margin - secondary_len;

  // for each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    // get the strand of the primary match
    primary_rc = get_array_item(0, sequence->primary_matches) < 0;
    // get the list of secondary matches
    secondary_array = matches[sequence->index];
    if (! secondary_array) continue;
    // process each secondary match
    int n_secondary_matches = get_array_length(secondary_array);
    for (i=0; i<n_secondary_matches; i++) {
      secondary = get_array_item(i, secondary_array);
      // convert the encoded form into easier to use form
      secondary_rc = secondary < 0;
      secondary_pos = (secondary_rc ? -secondary : secondary);
      // calculate the distance (counts from zero) and side
      if (secondary_pos <= margin) {
        distance = margin - secondary_pos - secondary_len + 1;
        quad = (primary_rc ? RIGHT : LEFT); //rotate reference direction
      } else {
        distance = secondary_pos - margin - primary_len - 1;
        quad = (primary_rc ? LEFT : RIGHT); //rotate reference direction
      }
      // check that we're within the acceptable range
      if (distance < 0 || distance > max_distance) {
        die("Secondary motif match not within margin as it should be due to prior checks!");
      }
      // calculate the strand
      quad |= (secondary_rc == primary_rc ? SAME : OPPO);
      // add a count to the frequencies
      spacing = smotif->spacings+(quad);
      bin = distance / bin_size;
      spacing->count[bin] += 1;
      linklst_add(sequence, spacing->sequences[bin]);
      smotif->total_spacings += 1;
    } // secondary match
  } // primary match
} // bin_matches

/**************************************************************************
 * Mark sequences that contribute to a significant spacing.
 *************************************************************************/
void mark_contributing_sequences(
  int bin,
  int orient,
  SECONDARY_MOTIF_T *smotif
) {
  int quad;
  LINKLST_T *list;
  LINK_T *link;
  SEQUENCE_T *sequence;

  // Mark contributing sequences in each quadrant
  for (quad = 0; quad < NQUADS; quad++) { 
    // Skip non-contributing quadrants.
    if (!QUADRANT_CONTRIBUTES(quad, orient)) continue;
    // now mark the sequences
    list = smotif->spacings[quad].sequences[bin];
    for (link = linklst_first(list); link != NULL; link = linklst_next(link)) {
      sequence = (SEQUENCE_T*)linklst_get(link);
      sequence->contributes = true;
    }
  }
}

/**************************************************************************
 * Create alignments of the contributing sequences in each
 * contributing quadrant in the orientation.
 *************************************************************************/
void create_alignment_motifs(
  int bin,
  int orient,
  MOTIF_T *pmotif,
  SECONDARY_MOTIF_T *smotif, 
  MOTIF_T **alignment_motif
)
{
  int quad, nquads = 0;
  ALPH_T *alph;
  alph = get_motif_alph(pmotif);

  // Make an alignment for each contributing quadrant.
  for (quad = 0; quad < NQUADS; quad++) { 

    // Skip non-contributing quadrants.
    if (!QUADRANT_CONTRIBUTES(quad, orient)) {
       alignment_motif[nquads++] = NULL;
       continue;
    }

    // Skip quadrant if it has no contributing sequences.
    LINKLST_T *contr_seq_list = smotif->spacings[quad].sequences[bin];
    LINK_T* link = linklst_first(contr_seq_list);
    if (link == NULL || linklst_size(contr_seq_list) == 0) {
      alignment_motif[nquads++] = NULL;
      continue;
    }

    // Make an alignment for the current quadrant.
    int num_seqs = linklst_size(contr_seq_list);
    SEQ_T** contr_seqs = (SEQ_T**)mm_malloc(num_seqs * sizeof(SEQ_T*));
    SEQUENCE_T *sequence = (SEQUENCE_T*) linklst_get(link);
    char *tmp_seq = (char*) mm_malloc((sequence->length+1) * sizeof(char));
    int plength = get_motif_trimmed_length(pmotif);
    int slength = get_motif_trimmed_length(smotif->motif);
    int i = 0;
    // Add each contributing sequence to the alignment.
    for ( ; link != NULL; link=linklst_next(link)) {
      sequence = (SEQUENCE_T*) linklst_get(link);
      int length = sequence->length;
      int center = length/2.0;
      int primary_start = center - (plength/2);
      int primary_end = primary_start + plength - 1;
      int first = primary_start - bin - slength;
      int last = primary_end + bin + slength;
      int pad = min(10, min(first, length-last-1));
      first -= pad;
      last += pad;
      int new_len = last - first + 1;
      // copy sequence segment, reverse-complementing if necessary
      int primary_match = (int) get_array_item(0, sequence->primary_matches);
      copy_string_with_rc(alph, &(sequence->data[first]), tmp_seq, new_len, (primary_match < 0));
      contr_seqs[i++] = allocate_seq(NULL, NULL, 0, tmp_seq);
    } // contributing sequence

    // Free space for temporary sequence.
    myfree(tmp_seq);
    myassert(TRUE, i == num_seqs, "The number of sequences (%d) should be (%d).\n", i, num_seqs);

    // Create the alignment.
    ALIGNMENT_T *alignment = allocate_alignment("", "", num_seqs, contr_seqs, "");

    // Convert alignment to frequency matrix.
    MATRIX_T *freqs = get_freq_matrix_from_alignment(alph, alignment);
    
    // Convert frequency matrix to motif.
    alignment_motif[nquads] = allocate_motif(LOGO_NAME[nquads], "", alph, freqs, NULL);
    set_motif_nsites(alignment_motif[nquads++], num_seqs);

    // Free the space for the alignment.
    for (i=0; i<num_seqs; i++) {
      free_seq(contr_seqs[i]);
    } 
    myfree(contr_seqs);
    free_alignment(alignment);

  } // quad

} //create_alignment_motifs

/**************************************************************************
* Create an alignment of the "inferred" secondary motif by aligning
* the contributing regions from the sequences with the given spacing.
* Returns the inferred motif.
quad:   0       1       2       3       4 (0,1) 5 (0,3) 6 (1,2) 7 (2,3) 8 (all)
orient: IV      III     I       II      III-IV  II-IV   I-III   I-II    all
 *************************************************************************/
MOTIF_T *create_inferred_secondary_motif(
  int bin,
  int orient,
  MOTIF_T *pmotif,
  SECONDARY_MOTIF_T *smotif
)
{
  int quad, i;
  int nseqs = 0;
  //int trim_left = get_motif_trim_left(smotif->motif);
  //int trim_right = get_motif_trim_right(smotif->motif);
  ALPH_T* alph = get_motif_alph(pmotif);
  

  // Get the total number of contributing sequences for the spacing.
  LINKLST_T *contr_seq_list = smotif->spacings[orient].sequences[bin];
  int num_seqs = linklst_size(contr_seq_list);

  // No sequences?
  // Create a dummy sequence containing all wildcards for the alignment
  // so every motif will have a valid secondary motif entry in HTML.
  SEQ_T** contr_seqs = (SEQ_T**)mm_malloc((num_seqs+1) * sizeof(SEQ_T*));
  if (num_seqs == 0) {
    int plength = get_motif_trimmed_length(pmotif);
    int slength = get_motif_trimmed_length(smotif->motif);
    int length = plength + slength;
    char *tmp_seq = (char*) mm_malloc((length+1) * sizeof(char));
    for (i = 0; i < length; i++) tmp_seq[i] = alph_wildcard(alph);
    contr_seqs[nseqs++] = allocate_seq(NULL, NULL, 0, tmp_seq);
    num_seqs = 1;
  }

  // Include segments from all contribution quadrants.
  for (quad=0; quad < NQUADS; quad++) { 

    // Skip non-contributing quadrants.
    if (!QUADRANT_CONTRIBUTES(quad, orient)) continue;

    // Skip quadrant if it has no sequences.
    contr_seq_list = smotif->spacings[quad].sequences[bin];
    LINK_T* link = linklst_first(contr_seq_list);
    if (link == NULL) continue;

    // Add the contributing sequence segments.
    SEQUENCE_T *sequence = (SEQUENCE_T*) linklst_get(link);
    char *tmp_seq = (char*) mm_malloc((sequence->length+1) * sizeof(char));
    int plength = get_motif_trimmed_length(pmotif);
    int slength = get_motif_trimmed_length(smotif->motif);

    // Add each contributing sequence segment to the alignment.
    for ( ; link != NULL; link=linklst_next(link)) {
      sequence = (SEQUENCE_T*) linklst_get(link);
      int length = sequence->length;
      int center = length/2.0;
      int primary_start = center - (plength/2);
      int primary_end = primary_start + plength - 1;
      int primary_match = (int) get_array_item(0, sequence->primary_matches);
      BOOLEAN_T rc = ( (primary_match>0 && STRAND(quad)==OPPO) || 
        (primary_match<0 && STRAND(quad)==SAME) );
      BOOLEAN_T up = ( (primary_match>0 && SIDE(quad)==LEFT) ||
        (primary_match<0 && SIDE(quad)==RIGHT) );
      //int trim = ((up && rc) || (!up && !rc)) ? trim_left : trim_right;
      int trim = 0;
      int first = up ? primary_start-bin-slength+trim : primary_end+bin+1-trim;
      // copy sequence segment, reverse-complementing if necessary
      copy_string_with_rc(alph, &(sequence->data[first]), tmp_seq, slength, rc);
      contr_seqs[nseqs++] = allocate_seq(NULL, NULL, 0, tmp_seq);
     
    } // contributing sequence

    // Free space for temporary sequence.
    myfree(tmp_seq);

  } // quad

  // Create the alignment.
  //myassert(TRUE, num_seqs == nseqs, "The total number of sequences (%d) should be (%d).\n", nseqs, num_seqs);
  ALIGNMENT_T *alignment = allocate_alignment("", "", nseqs, contr_seqs, "");

  // Convert alignment to frequency matrix.
  MATRIX_T *freqs = get_freq_matrix_from_alignment(alph, alignment);
  
  // Convert frequency matrix to motif.
  MOTIF_T *inferred_motif = allocate_motif(" inferrred secondary motif", "", alph, freqs, NULL);
  set_motif_nsites(inferred_motif, nseqs);

  // Free the space for the alignment.
  for (i=0; i<nseqs; i++) {
    free_seq(contr_seqs[i]);
  } 
  myfree(contr_seqs);
  free_alignment(alignment);

  return(inferred_motif);
} // create_inferred_secondary_motif

static void add_secondary_spacing(SECONDARY_MOTIF_T *smotif, MOTIF_T *pmotif, int orient, int bin, double pvalue) {
  SIGSPACE_T *sig;

  // create a new spacing object
  smotif->sig_count += 1;
  smotif->sigs = mm_realloc(smotif->sigs, sizeof(SIGSPACE_T) * smotif->sig_count);
  sig = smotif->sigs+(smotif->sig_count - 1);
  // set the values
  sig->pvalue = pvalue;
  sig->bin = bin;
  sig->orient = orient;
  // track which sequence names we must output
  mark_contributing_sequences(bin, orient, smotif);
  // create an alignment pwm for each of the quadrants
  create_alignment_motifs(bin, orient, pmotif, smotif, sig->alignment_motif);
  // combine the alignment into an inferred secondary motif
  sig->inferred_motif = create_inferred_secondary_motif(bin, orient, pmotif, smotif);
}

/**************************************************************************
* Compute a pvalue for a bin for the NORIENTS combinations of 
* 1, 2, or 4 quadrants that we condider biologically interesting.
*   |-------------|
 *   | IV   |   I  |
 *   | ------------|
 *   | III  |   II |
 *   |-------------|
 * quad:                0       1       2       3       4 (0,1) 5 (0,3) 6 (1,2) 7 (2,3) 8 (all)
 * orientation:         IV      III     I       II      III-IV  II-IV   I-III   I-II    all
 * combined quadrant names:
 *                                              	4 up-secondary-pal
 *                                                      	5 up-primary-pal
 *                                                             		6 down-primary-pal
 *                                                                      	7 down-secondary-pal
 *                                                                              	8 primary-secondary-pal
 * Store the all the p-values.
 **************************************************************************/
static void compute_spacing_pvalue(
    int tests,          // number of distances times number of combinations of quadrants
    double threshold,   // unadjusted p-value significance threshold
    int bin,            // the gap (or range of gaps) between primary and secondary
    int test_max,       // number of gaps (or range of gaps) being tested
    double prob,        // prior probability of the gap (or range of gaps) in a single quadrant
    MOTIF_T *pmotif,    // primary motif
    SECONDARY_MOTIF_T *smotif   // secondary motif
  ) {
  double pvalue;
  int counts[NORIENTS];
  double probs[NORIENTS];
  int quad, orient, i, total_spacings, nquads, norients;
  double log_tests = log(tests);
  bool revcomp;

  // check if the motifs have reverse complements
  revcomp = alph_has_complement(get_motif_alph(smotif->motif));
  nquads = (revcomp ? NQUADS : 2);
  norients = (revcomp ? NORIENTS : 2);
 
  // Update the maximum counts in one bin
  for (quad = 0; quad < nquads; quad++) {
    if (smotif->spacings[quad].count[bin] > smotif->max_in_one_bin) {
      smotif->max_in_one_bin = smotif->spacings[quad].count[bin];
    }
  }

  // Skip this bin (range of gaps) if beyond the requested range
  if (bin < test_max) {

    // Initialize the counts to 0.
    for(i = 0; i < NORIENTS; i++) {
      counts[i] = 0;
      probs[i] = 0.0;
    }

    // Number of observed occurrences of the given gap or range of gaps
    total_spacings = smotif->total_spacings;

    // Single quadrant counts and priors.
    for (quad = 0; quad < nquads; quad++) {
      counts[quad] = smotif->spacings[quad].count[bin];
      probs[quad] = prob;
    }

    if (revcomp) {
      // Double quadrant counts and priors
      smotif->spacings[UP_SEC_PAL].count[bin] = counts[UP_SEC_PAL] = counts[LEFT | SAME] + counts[LEFT | OPPO];
      smotif->spacings[UP_SEC_PAL].sequences[bin] = linklst_plus(smotif->spacings[LEFT | SAME].sequences[bin], smotif->spacings[LEFT | OPPO].sequences[bin]);

      smotif->spacings[UP_PRI_PAL].count[bin] = counts[UP_PRI_PAL] = counts[LEFT | SAME] + counts[RIGHT | OPPO];
      smotif->spacings[UP_PRI_PAL].sequences[bin] = linklst_plus(smotif->spacings[LEFT | SAME].sequences[bin], smotif->spacings[RIGHT | OPPO].sequences[bin]);

      smotif->spacings[DOWN_PRI_PAL].count[bin] = counts[DOWN_PRI_PAL] = counts[RIGHT | SAME] + counts[LEFT | OPPO];
      smotif->spacings[DOWN_PRI_PAL].sequences[bin] = linklst_plus(smotif->spacings[RIGHT | SAME].sequences[bin], smotif->spacings[LEFT | OPPO].sequences[bin]);

      smotif->spacings[DOWN_SEC_PAL].count[bin] = counts[DOWN_SEC_PAL] = counts[RIGHT | SAME] + counts[RIGHT | OPPO];
      smotif->spacings[DOWN_SEC_PAL].sequences[bin] = linklst_plus(smotif->spacings[RIGHT | SAME].sequences[bin], smotif->spacings[RIGHT | OPPO].sequences[bin]);

      probs[UP_SEC_PAL] = probs[UP_PRI_PAL] = probs[DOWN_PRI_PAL] = probs[DOWN_SEC_PAL] = 2 * prob;

      // "All" quadrant counts and priors
      smotif->spacings[BOTH_PAL].count[bin] = counts[BOTH_PAL] = counts[UP_SEC_PAL] + counts[DOWN_SEC_PAL];
      smotif->spacings[BOTH_PAL].sequences[bin] = linklst_plus(smotif->spacings[UP_SEC_PAL].sequences[bin], smotif->spacings[DOWN_SEC_PAL].sequences[bin]);
      probs[BOTH_PAL] = 4 * prob;
    }

    // Compute the unadjusted and adusted p-value for each orientation.
    for (orient = 0; orient < norients; orient++) {
      // unadjusted p-value for a given gap and orientation
      pvalue = one_minus_binomial_cdf(probs[orient], counts[orient], total_spacings);
      // adjust p-value for the number of distances (gaps) being tested 
      // times the number of combinations of quadrants 
      double adj_pvalue = exp(LOGEV(log_tests, log(pvalue)));
      smotif->spacings[orient].pvalue[bin] = adj_pvalue;
      // Save the best adjusted p-value over all tested spacings/orientations of secondary motif
      if (adj_pvalue < smotif->min_pvalue) {
        smotif->min_pvalue = adj_pvalue;
        smotif->best_orient = orient;
        smotif->best_bin = bin;
      }
      // Save the spacing if significant.
      if (adj_pvalue <= threshold) add_secondary_spacing(smotif, pmotif, orient, bin, pvalue);
    } // orientation
  } else {
    //don't test this bin, set pvalue in orient 0 to special value to indicate it is untested
    pvalue = 2;
    smotif->spacings[0].pvalue[bin] = pvalue;
  }

} // compute_spacing_pvalue

/**************************************************************************
 * compare the pvalues for two significant spacings
 **************************************************************************/
int compare_sigs(const void *p1, const void *p2) {
  SIGSPACE_T *s1 = (SIGSPACE_T*)p1;
  SIGSPACE_T *s2 = (SIGSPACE_T*)p2;
  if (s1->pvalue < s2->pvalue) {
    return -1;
  } else if (s1->pvalue == s2->pvalue) {
    if (s1->orient < s2->orient) {
      return -1;
    } else if (s1->orient == s2->orient) {
      if (s1->bin < s2->bin) {
        return -1;
      } else if (s1->bin == s2->bin) {
        return 0; // shouldn't happen as they should be unique
      } else {
        return 1;
      }
    } else {
      return 1;
    }
  } else {
    return 1;
  }
}

/**************************************************************************
 * compute the pvalues for the frequencies of each spacing in each quadrant
 **************************************************************************/
void compute_spacing_stats(int margin, int bin_size, int n_secondary_motifs, int test_max, 
  double threshold, double motif_evalue_cutoff, MOTIF_T *pmotif, SECONDARY_MOTIF_T *smotif) {
  int quad_opt_count, quad_bin_count, quad_leftover, total_opt_count, i, j;
  double general_prob, leftover_prob;
  bool revcomp;
  // check if reverse complement motifs are available
  revcomp = alph_has_complement(get_motif_alph(pmotif));
  //the number of possible values for spacings in one quadrant
  quad_opt_count = margin - get_motif_trimmed_length(smotif->motif) + 1;
  //the number of bins in one quadrant (excluding a possible leftover bin)
  quad_bin_count = (int)(quad_opt_count / bin_size);
  //the number of spacings that don't fit in the full bins (the number that would go into the leftover bin)
  quad_leftover = quad_opt_count % bin_size;
  //the total number of possible values for spacings
  total_opt_count = quad_opt_count * (revcomp ? 4 : 2);
  //prior probability of a bin that has bin_size possible spacings that could go into it
  general_prob = (double)bin_size / total_opt_count;
  //prior probability of the final bin that has less than bin_size possible spacings that could go into it
  leftover_prob = (double)quad_leftover / total_opt_count;
  //calculate the number of independent tests for NORIENTS combinations of 4 quadrants
  int tests = (revcomp ? NORIENTS : 2) * (min(quad_bin_count,test_max) + (quad_leftover == 0 ? 0 : 1));
  //calculate the significance of each bin
  for (j = 0; j < quad_bin_count; ++j) {
    compute_spacing_pvalue(tests, threshold, j, test_max, general_prob, pmotif, smotif);
  }
  if (quad_leftover) { // bin only exists if quad_leftover is non-zero
    compute_spacing_pvalue(tests, threshold, j, test_max, leftover_prob, pmotif, smotif);
  }
  // judge if the motif passes the threshold
  smotif->passes_evalue_cutoff = ((smotif->min_pvalue * n_secondary_motifs) <= motif_evalue_cutoff);
  // when the motif is significant but no spacings pass the significance test track the best one
  if (smotif->passes_evalue_cutoff && smotif->sigs == 0) {
    add_secondary_spacing(smotif, pmotif, smotif->best_orient, smotif->best_bin, smotif->min_pvalue);
  }
  //sort the significant finds
  qsort(smotif->sigs, smotif->sig_count, sizeof(SIGSPACE_T), compare_sigs);
}

/**************************************************************************
 * compute the list of ids for the most significant spacing
 **************************************************************************/
void compute_idset(int margin, int bin_size, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, ARRAY_T **matches) {
  int i, primary_len, secondary_len, secondary, secondary_pos, primary_rc, secondary_rc, quad, distance;
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  ARRAY_T *secondary_array;

  if (secondary_motif->sig_count == 0) return;

  primary_len = get_motif_trimmed_length(primary_motif);
  secondary_len = get_motif_trimmed_length(secondary_motif->motif);

  // for each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    secondary_array = matches[sequence->index];
    if (! secondary_array) continue;
    // process each secondary match
    int n_secondary_matches = get_array_length(secondary_array);
    for (i=0; i<n_secondary_matches; i++) {
      secondary = get_array_item(i, secondary_array);
      // convert the encoded form into easier to use form
      primary_rc = get_array_item(0, sequence->primary_matches) < 0;
      secondary_rc = secondary < 0;
      secondary_pos = (secondary_rc ? -secondary : secondary);
      // calculate the distance and side
      // note that distance can be zero meaning the primary is next to the secondary
      if (secondary_pos <= margin) {
        distance = margin - secondary_pos - secondary_len + 1;
        quad = LEFT;
      } else {
        distance = secondary_pos - margin - primary_len;
        quad = RIGHT;
      }
      // calculate the strand
      quad |= (secondary_rc == primary_rc ? SAME : OPPO);
      // add the sequence id to the set if the bin and quadrant matches    
      // the most significant bin and orientation
      int bin = secondary_motif->sigs[0].bin;
      int orient = secondary_motif->sigs[0].orient;
      if ( ((distance / bin_size) == bin) && QUADRANT_CONTRIBUTES(quad, orient) ) {
        secondary_motif->seq_count += 1;
        secondary_motif->seqs = (int*)mm_realloc(secondary_motif->seqs, sizeof(int) * secondary_motif->seq_count);
        secondary_motif->seqs[secondary_motif->seq_count-1] = sequence->index;
      }
    } // secondary match
  } // primary match
} // compute_idset

/**************************************************************************
 * This does most of the calculation steps after the matching
 * positions of the motifs have been derived from a scan
 * file. The steps undertaken are:
 * 1) bin the matches
 * 2) compute the pvalues of the tested region around the primary
 * 3) compute the sequence id set of the most significant peak
 *
 **************************************************************************/
void process_matches(
  int margin, 
  int bin, 
  double significance_threshold, 
  double motif_evalue_cutoff,
  int test_max,
  MOTIF_T *primary_motif, 
  RBTREE_T *sequences,
  SECONDARY_MOTIF_T *secondary_motif,
  int n_secondary_motifs,
  ARRAY_T **secondary_matches
) {
  bin_matches(margin, bin, sequences, primary_motif, secondary_motif, secondary_matches);
  compute_spacing_stats(margin, bin, n_secondary_motifs, test_max, significance_threshold,
    motif_evalue_cutoff, primary_motif, secondary_motif);
  compute_idset(margin, bin, sequences, primary_motif, secondary_motif, secondary_matches);
}
