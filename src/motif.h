/***********************************************************************
 * FILE: motif.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-13-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Data structure for representing one motif.
 ***********************************************************************/
#ifndef MOTIF_H
#define MOTIF_H

#include "matrix.h" // include before array.h

#include "alphabet.h"
#include "array-list.h"
#include "red-black-tree.h"
#include "utils.h"
#include "mtwist.h"

// For reading and writing, store null ID chars as this.
#define NON_MOTIF_ID_CHAR '.'

// Maximum number of characters in a motif ID.
#define MAX_MOTIF_ID_LENGTH 100
#define MAX_MOTIF_URL_LENGTH 500

typedef struct motif_t MOTIF_T;

/***********************************************************************
 * Set a boolean to true on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
void set_motif_mark
  (MOTIF_T* motif, int mark_no);

/***********************************************************************
 * Set a boolean to false on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
void clear_motif_mark
  (MOTIF_T* motif, int mark_no);

/***********************************************************************
 * Test a boolean on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
bool test_motif_mark
  (MOTIF_T* motif, int mark_no);

/***********************************************************************
 * Get the index of the motif in the file. The first motif is indexed 1.
 ***********************************************************************/
int get_motif_idx
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif id.
 ***********************************************************************/
char* get_motif_id
  (MOTIF_T* motif);

char* get_motif_id2
  (MOTIF_T* motif);

/***********************************************************************
 * Get set the consensus of a motif.
 ***********************************************************************/
char *get_motif_consensus
  (MOTIF_T* motif);

/***********************************************************************
 * Get the strand of a motif.
 ***********************************************************************/
char get_motif_strand
  (MOTIF_T *motif);

/***********************************************************************
 * Get motif id with leading +/- indicating strand.
 * For protein motifs this acts like get_motif_id.
 ***********************************************************************/
char* get_motif_st_id
  (MOTIF_T *motif);

/***********************************************************************
 * Get the frequencies
 ***********************************************************************/
MATRIX_T* get_motif_freqs
  (MOTIF_T* motif);

/***********************************************************************
 * Get the scores
 ***********************************************************************/
MATRIX_T* get_motif_scores
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif length
 ***********************************************************************/
int get_motif_length
  (const MOTIF_T* motif);

/***********************************************************************
 * Get the motif length after trimming
 ***********************************************************************/
int get_motif_trimmed_length
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif alphabet
 ***********************************************************************/
ALPH_T* get_motif_alph
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif strands
 * TRUE if the motif was created looking at both strands
 ***********************************************************************/
BOOLEAN_T get_motif_strands
  (MOTIF_T *motif);

/***********************************************************************
 * Get the motif alphabet size (non-ambiguous letters)
 ***********************************************************************/
int get_motif_alph_size
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif ambiguous alphabet size (only ambiguous letters)
 ***********************************************************************/
int get_motif_ambiguous_size
  (MOTIF_T* motif);

/***********************************************************************
 * Get the E-value of a motif.
 ***********************************************************************/
double get_motif_evalue
  (MOTIF_T* motif);

/***********************************************************************
 * Get the log E-value of a motif.
 ***********************************************************************/
double get_motif_log_evalue
  (MOTIF_T* motif);

/***********************************************************************
 * Get the complexity of a motif.
 ***********************************************************************/
double get_motif_complexity
  (MOTIF_T *motif);

/***********************************************************************
 * Get the number of sites of a motif.
 ***********************************************************************/
double get_motif_nsites
  (MOTIF_T* motif);

/***********************************************************************
 * Get the position specific score for a letter
 ***********************************************************************/
double get_motif_score
  (MOTIF_T* motif, int position, int i_alph);

/***********************************************************************
 * Return one column of a motif, as a newly allocated array of counts.
 ***********************************************************************/
ARRAY_T* get_motif_counts
  (int      position,
   MOTIF_T* motif);

/***********************************************************************
 * Get the url of a motif
 ***********************************************************************/
char* get_motif_url
  (MOTIF_T* motif);

/***********************************************************************
 * Check if the motif has a URL
 ***********************************************************************/
BOOLEAN_T has_motif_url
  (MOTIF_T *motif);

/***********************************************************************
 * Get the number of positions to trim from the left of the motif
 ***********************************************************************/
int get_motif_trim_left
  (MOTIF_T *motif);

/***********************************************************************
 * Get the number of positions to trim from the right of the motif
 ***********************************************************************/
int get_motif_trim_right
  (MOTIF_T *motif);

/***********************************************************************
 * Clear the motif trim
 ***********************************************************************/
void clear_motif_trim
  (MOTIF_T *motif);

/***********************************************************************
 * Check the motif to see it it has any probabilities that are zero.
 * Some algorithms can not handle motifs with probabilities that are zero.
 ***********************************************************************/
bool has_motif_zeros
  (MOTIF_T *motif);

/***********************************************************************
 * Determine whether a given motif is in a given list of motifs.
 ***********************************************************************/
BOOLEAN_T have_motif
  (char*    motif_id,
   int      num_motifs,
   MOTIF_T* motifs);

/***********************************************************************
 * Copy a motif from one place to another.
 ***********************************************************************/
void copy_motif
  (MOTIF_T* source,
   MOTIF_T* dest);

/***********************************************************************
 * Allocates a new matrix with the columns rearranged to suite the target
 * alphabet. Any missing columns are filled in with the specified value.
 *
 * Note that the target alphabet must have all the primary core symbols
 * of the source alphabet defined as a core symbol.
 ***********************************************************************/
MATRIX_T* convert_matrix_alphabet
  (MATRIX_T *in, 
   MTYPE value, 
   ALPH_T *source_alph, 
   ALPH_T *target_alph);

/***********************************************************************
 * Shuffle the positions of the motif
 ***********************************************************************/
void shuffle_motif
  (MOTIF_T* motif, mt_state* prng);

/***********************************************************************
 * Takes a matrix of letter probabilities and converts them into meme
 * score. This adds the supplied pseudo count to avoid log of zero.
 *
 * Assuming the probability is nonzero the score is just: 
 * s = log2(p / bg) * 100
 *
 ***********************************************************************/
MATRIX_T* convert_freqs_into_scores
  (ALPH_T *alph,
   MATRIX_T *freqs,
   ARRAY_T *bg,
   int site_count,
   double pseudo_count);

/***********************************************************************
 * Takes a matrix of meme scores and converts them into letter 
 * probabilities. This assumes that the scores had pseudo counts applied
 * when they were originally created and attempts to reverse that.
 *
 * The probablility can be got by:
 * p = (2 ^ (s / 100)) * bg
 *
 ***********************************************************************/
MATRIX_T* convert_scores_into_freqs
  (ALPH_T *alph,
   MATRIX_T *scores,
   ARRAY_T *bg,
   int site_count,
   double pseudo_count);

/***********************************************************************
 * Turn a given motif into its own reverse complement.
 ***********************************************************************/
void reverse_complement_motif
  (MOTIF_T* a_motif);

/***********************************************************************
 * Apply a pseudocount to the motif pspm.
 ***********************************************************************/
void apply_pseudocount_to_motif
  (MOTIF_T* motif, ARRAY_T *background, double pseudocount);

/***********************************************************************
 * Calculate the ambiguous letters from the concrete ones.
 ***********************************************************************/
void calc_motif_ambigs
  (MOTIF_T *motif);

/***********************************************************************
 * Normalize the motif's pspm
 ***********************************************************************/
void normalize_motif
  (MOTIF_T *motif, double tolerance);

/***********************************************************************
 * Set the trimming bounds on the motif.
 *
 * Reads from the left and right until it finds a motif position with
 * an information content larger or equal to the threshold in bits.
 * 
 ***********************************************************************/
void trim_motif_by_bit_threshold(
  MOTIF_T *a_motif, 
  double threshold_bits
);

/***********************************************************************
 * Compute the complexity of a motif as a number between 0 and 1.
 ***********************************************************************/
double compute_motif_complexity
  (MOTIF_T* a_motif);

/***********************************************************************
 * Compute the number of positions from the start or end of a motif
 * that contain a given percentage of the information content.
 *
 * Information content is the same as relative entropy, and is computed
 * as
 *
 *  \sum_i p_i log(p_i/f_i)
 *
 ***********************************************************************/
int get_info_content_position
  (BOOLEAN_T from_start, // Count from start?  Otherwise, count from end.
   float     threshold,  // Information content threshold (in 0-100).
   ARRAY_T*  background, // Background distribution.
   MOTIF_T*  a_motif);

/***********************************************************************
 * Returns the string that is the best possible match to the given motif.
 ***********************************************************************/
char *get_best_possible_match(MOTIF_T *motif);

/***********************************************************************
 * Duplicates and reverse complements the motif
 ***********************************************************************/
MOTIF_T* dup_rc_motif
  (MOTIF_T *motif);

/***********************************************************************
 * Duplicates the motif
 ***********************************************************************/
MOTIF_T* duplicate_motif
  (MOTIF_T *motif);

/***********************************************************************
 * Free dynamic memory used by a given motif assuming that the
 * structure itself does not need to be freed.
 ***********************************************************************/
void free_motif
  (MOTIF_T * a_motif);

/***********************************************************************
 * Free dynamic memory used by a given motif and free the structure.
 * To be useable by collections it takes a void * but expects
 * a MOTIF_T *.
 ***********************************************************************/
void destroy_motif
  (void * a_motif);

/***********************************************************************
 * Convert a list of motifs into an array of motifs with a count.
 ***********************************************************************/
void motif_list_to_array(ARRAYLST_T *motif_list, MOTIF_T **motif_array, int *num);

/***********************************************************************
 * Convert a tree of motifs into an array of motifs with a count.
 ***********************************************************************/
void motif_tree_to_array(RBTREE_T *motif_tree, MOTIF_T **motif_array, int *num);

/***********************************************************************
 * Get the motif at the selected index in the array 
 ***********************************************************************/
MOTIF_T* motif_at(MOTIF_T *array_of_motifs, int index);

/***********************************************************************
 * Free an array of motifs
 ***********************************************************************/
void free_motif_array(MOTIF_T *motif_array, int num);

/***********************************************************************
 * Free an array list and the contained motifs
 ***********************************************************************/
void free_motifs(ARRAYLST_T *motif_list);

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 *
 * The reverse complement is always listed directly after the original.
 ***********************************************************************/
void add_reverse_complements(ARRAYLST_T* motifs);

/***********************************************************************
*
* Dump the frequencies of the motif to the output.
*
 ***********************************************************************/
void dump_motif_freqs(FILE *out, MOTIF_T* m);

/***********************************************************************
 *
 * Check that all alphabets are the same or convertible in a set of
 * motif files
 *
 ***********************************************************************/
void read_motif_alphabets(ARRAYLST_T* motif_sources, bool xalph, ALPH_T** alph);

/**************************************************************************
 * Compares two AP_T and orders largest probability to smallest probability.
 **************************************************************************/
int ap_cmp(const void *p1, const void *p2);

/**************************************************************************
 * Compares two uint8_t and orders smallest to largest.
 **************************************************************************/
int idx_cmp(const void *p1, const void *p2);

/**************************************************************************
 * Attempts to make a reasonable consensus representation of a motif.
 **************************************************************************/
void motif2consensus(MOTIF_T* motif, STR_T* consensus, bool single_letter);
#endif

