/***********************************************************************
 * FILE: motif.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-13-97
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Data structure for representing one motif.
 ***********************************************************************/
#include <assert.h>
#include <math.h>
#include <stdlib.h> /* For definition of NULL .*/
#include <string.h>

#include "macros.h"
#include "motif.h"
#include "motif-spec.h"
#include "motif-in.h"

// same default site count as Tomtom's website
#define DEFAULT_SITE_COUNT 20

/**************************************************************************
 * Holds a probability for an alphabet index. This is useful because a 
 * list of these can be sorted largest to smallest and the indexes can be
 * extracted.
 **************************************************************************/
typedef struct alph_prob {
  uint8_t idx;
  double prob;
} AP_T;


/**************************************************************************
 * Create a number in the range 0 to n (inclusive).
 **************************************************************************/
static inline uint32_t rand_uint32(mt_state* prng, const uint32_t n) {
  uint64_t d, threshold;
  uint32_t r;
  assert(n > 0);
  d = (uint64_t)n + 1;
  threshold = (uint64_t)UINT32_MAX + 1;
  threshold = threshold - (threshold % d);
  assert(threshold > 0);
  do {
    r = mts_lrand(prng);
  } while (r >= threshold);
  return r % d;
}

/***********************************************************************
 * Allocate memory for a MEME motif and initialize from score and
 * frequency matrices.
 ***********************************************************************/
MOTIF_T* allocate_motif(
  char *id,
  char *id2,
  ALPH_T *alph,
  MATRIX_T* freqs,
  MATRIX_T* scores
){
  MOTIF_T* motif = mm_malloc(sizeof(MOTIF_T));

  assert(id != NULL);
  assert(id2 != NULL);
  if (freqs == NULL && scores == NULL) {
    die(
      "A matrix of scores, or frequencies, or both, "
      "must be provided when allocating a motif.\n"
    );
  }
  set_motif_strand('+', motif);
  set_motif_id(id, strlen(id), motif);
  set_motif_id2(id2, strlen(id2), motif);
  //motif->consensus = NULL;
  motif->length = freqs ? get_num_rows(freqs) : get_num_rows(scores);
  motif->alph = alph_hold(alph);
  motif->flags = 0;
  motif->evalue = 0.0;
  motif->log_evalue = -HUGE_VAL;
  motif->num_sites = 0.0;
  motif->complexity = 0.0;
  motif->freqs = freqs ? duplicate_matrix(freqs) : NULL;
  motif->scores = scores ? duplicate_matrix(scores) : NULL;
  motif->url = NULL;
  motif->trim_left = 0;
  motif->trim_right = 0;

  // Compute a single-letter consensus for the motif. 
  STR_T *cons_buf = str_create(50);
  str_clear(cons_buf);
  motif2consensus(motif, cons_buf, true);
  motif->consensus = str_internal(cons_buf);

  return motif;
}

/***********************************************************************
 * Calculates the information content of a position of the motif.
 ***********************************************************************/
static inline double position_information_content(
  MOTIF_T *a_motif,
  int position
) {
  int i, asize;
  double H, item;
  ARRAY_T *freqs;

  asize = alph_size_core(a_motif->alph);
  H = 0;
  freqs = get_matrix_row(position, a_motif->freqs);
  for (i = 0; i < asize; ++i) {
    item = get_array_item(i, freqs);
    H -= item*my_log2(item);
  }
  return my_log2(asize) - H;
}

#define INT_BIT (((sizeof(int) / sizeof(char)) * CHAR_BIT) - 1)

/***********************************************************************
 * Set a boolean to true on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
void set_motif_mark
  (MOTIF_T* motif, int mark_no) 
{
  int bit_pos = mark_no + MOTIF_RESERVED_FLAGS;
  if (mark_no < 0 || bit_pos > INT_BIT) die("Insufficent motif flag bits to mark %d", mark_no);
  motif->flags |= (1 << bit_pos);
}

/***********************************************************************
 * Set a boolean to false on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
void clear_motif_mark
  (MOTIF_T* motif, int mark_no)
{
  int bit_pos = mark_no + MOTIF_RESERVED_FLAGS;
  if (mark_no < 0 || bit_pos > INT_BIT) die("Insufficent motif flag bits to mark %d", mark_no);
  motif->flags &= ~(1 << bit_pos);
}

/***********************************************************************
 * Test a boolean on the motif object.
 * This exists to allow applications to store simple state on a motif.
 ***********************************************************************/
bool test_motif_mark
  (MOTIF_T* motif, int mark_no) 
{
  int bit_pos = mark_no + MOTIF_RESERVED_FLAGS;
  if (mark_no < 0 || bit_pos > INT_BIT) die("Insufficent motif flag bits to mark %d", mark_no);
  return ((motif->flags & (1 << bit_pos)) != 0);
}

/***********************************************************************
 * Get the index of the motif in the file. The first motif is indexed 1.
 ***********************************************************************/
int get_motif_idx
  (MOTIF_T* motif)
{
  return motif->idx;
}

/***********************************************************************
 * Get or set the identifier of a motif.
 ***********************************************************************/
char* get_motif_id
  (MOTIF_T* motif)
{
  return motif->id + 1; // without strand indicator
}

char* get_motif_st_id
  (MOTIF_T* motif) 
{
  // skip strand indicator for alphabets that don't have a complement
  if (!alph_has_complement(motif->alph)) return motif->id+1;
  return motif->id; // with strand indicator
}

void set_motif_id
  (const char* id,
   int len,
   MOTIF_T* motif)
{
  // the first character is the strand indicator which defaults to '?'
  if (motif->id[0] != '+' && motif->id[0] != '-') motif->id[0] = '?';
  len = (len < MAX_MOTIF_ID_LENGTH ? len : MAX_MOTIF_ID_LENGTH);
  strncpy(motif->id+1, id, len);
  motif->id[len+1] = '\0'; // ensure null termination
}

char* get_motif_id2
  (MOTIF_T* motif)
{
  return motif->id2;
}

void set_motif_id2
  (const char* id2,
   int len,
   MOTIF_T* motif)
{
  len = (len < MAX_MOTIF_ID_LENGTH ? len : MAX_MOTIF_ID_LENGTH);
  strncpy(motif->id2, id2, len);
  motif->id2[len] = '\0'; // ensure null termination
}

/***********************************************************************
 * Get set the consensus of a motif.
 ***********************************************************************/
char *get_motif_consensus
  (MOTIF_T* motif)
{
  return motif->consensus;
}

/***********************************************************************
 * Get or set the strand of a motif.
 ***********************************************************************/
void set_motif_strand
  (char strand,
   MOTIF_T *motif)
{
  assert(strand == '?' || strand == '-' || strand == '+');
  motif->id[0] = strand;
}

char get_motif_strand
  (MOTIF_T *motif)
{
  return motif->id[0];
}

/***********************************************************************
 * Set the number of sites of a motif.
 ***********************************************************************/
void set_motif_nsites
  (MOTIF_T* motif, int num_sites)
{
  motif->num_sites = num_sites;
}

/***********************************************************************
 * Get the frequencies
 ***********************************************************************/
MATRIX_T* get_motif_freqs(MOTIF_T* motif){
  return motif->freqs;
}

/***********************************************************************
 * Get the scores
 ***********************************************************************/
MATRIX_T* get_motif_scores(MOTIF_T* motif){
  return motif->scores;
}

/***********************************************************************
 * Get the E-value of a motif.
 ***********************************************************************/
double get_motif_evalue
  (MOTIF_T* motif)
{
  return motif->evalue;
}

/***********************************************************************
 * Get the log (base 10) E-value of a motif.
 ***********************************************************************/
double get_motif_log_evalue
  (MOTIF_T* motif)
{
  return motif->log_evalue;
}

/***********************************************************************
 * Get the complexity of a motif.
 ***********************************************************************/
double get_motif_complexity
  (MOTIF_T *motif)
{
  return motif->complexity;
}

/***********************************************************************
 * Get the number of sites of a motif.
 ***********************************************************************/
double get_motif_nsites
  (MOTIF_T* motif)
{
  return motif->num_sites;
}

/***********************************************************************
 * Get the motif length
 ***********************************************************************/
int get_motif_length
  (const MOTIF_T* motif) 
{
  return motif->length;
}

/***********************************************************************
 * Get the motif alphabet
 ***********************************************************************/
ALPH_T* get_motif_alph
  (MOTIF_T* motif) 
{
  return motif->alph;
}

/***********************************************************************
 * Get the motif strands
 * TRUE if the motif was created looking at both strands
 ***********************************************************************/
BOOLEAN_T get_motif_strands
  (MOTIF_T *motif)
{
  return ((motif->flags & MOTIF_BOTH_STRANDS) != 0);
}

/***********************************************************************
 * Get the motif length after trimming
 ***********************************************************************/
int get_motif_trimmed_length
  (MOTIF_T* motif)
{
  return motif->length - motif->trim_left - motif->trim_right;
}

/***********************************************************************
 * Get the motif alphabet size (non-ambiguous letters)
 ***********************************************************************/
int get_motif_alph_size
  (MOTIF_T* motif) 
{
  return alph_size_core(motif->alph);
}


/***********************************************************************
 * Get the motif ambiguous alphabet size (only ambiguous letters)
 ***********************************************************************/
int get_motif_ambiguous_size
  (MOTIF_T* motif)
{
  return (motif->flags & MOTIF_HAS_AMBIGS ? alph_size_ambig(motif->alph) : 0);
}

/***********************************************************************
 * Get the position specific score for a letter
 ***********************************************************************/
double get_motif_score
  (MOTIF_T* motif, int position, int i_alph)
{
  return get_matrix_cell(position, i_alph, motif->scores);
}

/***********************************************************************
 * Return one column of a motif, as a newly allocated array of counts.
 * This assumes that num_sites is a reasonable value and not zero...
 ***********************************************************************/
ARRAY_T* get_motif_counts
  (int      position,
   MOTIF_T* motif)
{
  int i_alph, asize;
  ARRAY_T* return_value;
  
  asize = alph_size_core(motif->alph);
  return_value = allocate_array(asize);

  for (i_alph = 0; i_alph < asize; i_alph++) {
    set_array_item(i_alph, motif->num_sites * 
        get_matrix_cell(position, i_alph, motif->freqs), return_value);
  }
  return(return_value);
}

/***********************************************************************
 * Get the url of a motif
 ***********************************************************************/
char* get_motif_url
  (MOTIF_T *motif)
{
  return motif->url;
}

/***********************************************************************
 * Set the url of a motif
 ***********************************************************************/
void set_motif_url
  (char    *url,
   MOTIF_T *motif)
{
  if (motif->url) {
    free(motif->url);
    motif->url = NULL;
  }
  copy_string(&(motif->url), url);
}

/***********************************************************************
 * Check if the motif has a URL
 ***********************************************************************/
BOOLEAN_T has_motif_url
  (MOTIF_T *motif)
{
  if (!motif->url) return FALSE;
  if (motif->url[0] == '\0') return FALSE;
  return TRUE;
}

/***********************************************************************
 * Get the number of positions to trim from the left of the motif
 ***********************************************************************/
int get_motif_trim_left
  (MOTIF_T *motif)
{
  return motif->trim_left;
}

/***********************************************************************
 * Get the number of positions to trim from the right of the motif
 ***********************************************************************/
int get_motif_trim_right
  (MOTIF_T *motif)
{
  return motif->trim_right;
}

/***********************************************************************
 * Clear the motif trim
 ***********************************************************************/
void clear_motif_trim
  (MOTIF_T *motif)
{
  motif->trim_left = 0;
  motif->trim_right = 0;
}

/***********************************************************************
 * Check the motif to see it it has any probabilities that are zero.
 * Some algorithms can not handle motifs with probabilities that are zero.
 ***********************************************************************/
bool has_motif_zeros
  (MOTIF_T *motif)
{
  int row, col;
  for (row = 0; row < get_num_rows(motif->freqs); row++) {
    for (col = 0; col < get_num_cols(motif->freqs); col++) {
      if (get_matrix_cell(row, col, motif->freqs) == 0) {
        return true;
      }
    }
  }
  return false;
}


/***********************************************************************
 * Determine whether a given motif is in a given list of motifs.
 ***********************************************************************/
BOOLEAN_T have_motif
  (char*    motif_id,
   int      num_motifs,
   MOTIF_T* motifs)
{
  int i_motif;

  for (i_motif = 0; i_motif < num_motifs; i_motif++) {
    if (strcmp(motifs[i_motif].id, motif_id) == 0) {
      return(TRUE);
    }
  }
  
  return(FALSE);

}

/***********************************************************************
 * Copy a motif from one place to another.
 ***********************************************************************/
void copy_motif
  (MOTIF_T* source,
   MOTIF_T* dest)
{
  int size;
  memset(dest, 0, sizeof(MOTIF_T));
  strcpy(dest->id, source->id);
  strcpy(dest->id2, source->id2);
  if (source->consensus) dest->consensus = strdup(source->consensus);
  dest->length = source->length;
  dest->alph = alph_hold(source->alph);
  dest->flags = source->flags;
  dest->evalue = source->evalue;
  dest->log_evalue = source->log_evalue;
  dest->num_sites = source->num_sites;
  dest->complexity = source->complexity;
  if (source->freqs) {
    size = (dest->flags & MOTIF_HAS_AMBIGS ?
        alph_size_full(dest->alph) : 
        alph_size_core(dest->alph));
    // Allocate memory for the matrix.
    dest->freqs = allocate_matrix(dest->length, size);
    // Copy the matrix.
    copy_matrix(source->freqs, dest->freqs);
  } else {
    dest->freqs = NULL;
  }
  if (source->scores) {
    // Allocate memory for the matrix. Note that scores don't contain ambigs.
    dest->scores = allocate_matrix(dest->length, alph_size_core(dest->alph));
    // Copy the matrix.
    copy_matrix(source->scores, dest->scores);
  } else {
    dest->scores = NULL;
  }
  if (dest->url != NULL) {
    free(dest->url);
    dest->url = NULL;
  }
  copy_string(&(dest->url), source->url);
  dest->trim_left = source->trim_left;
  dest->trim_right = source->trim_right;
}

/***********************************************************************
 * Allocates a new matrix with the columns rearranged to suite the target
 * alphabet. Any missing columns are filled in with the specified value.
 *
 * Note that the target alphabet must have all the primary core symbols
 * of the source alphabet defined as a core symbol.
 ***********************************************************************/
MATRIX_T* convert_matrix_alphabet(MATRIX_T *in, MTYPE value, ALPH_T *source_alph, ALPH_T *target_alph) {
  MATRIX_T *out;
  int i_col_from, i_col_to, num_rows, i_row;
  uint32_t bitset[4] = {0, 0, 0, 0};
  num_rows = get_num_rows(in);
  // create a new matrix the size of the extended alphabet
  out = allocate_matrix(num_rows, alph_size_core(target_alph));
  init_matrix(value, out);
  // copy the old matrix to the new one
  for (i_col_from = 0; i_col_from < alph_size_core(source_alph); i_col_from++) {
    // find the new index
    i_col_to = alph_indexc(target_alph, alph_char(source_alph, i_col_from));
    if (i_col_to < 0) {
      die("Failed to promote matrix from %s to %s because %c is missing.",
          alph_name(source_alph), alph_name(target_alph),
          alph_char(source_alph, i_col_from));
      return NULL; // placate compiler
    } else if (i_col_to >= 128) {
      die("Alphabet index is too large! This should not be possible");
      return NULL; // placate compiler
    }
    // track which columns we've used
    if ((bitset[i_col_to / 32] & (1 << (i_col_to % 32))) != 0) {
      die("Failed to promote matrix from %s to %s because %c is becomes the "
          "same column as another core symbol.",
          alph_name(source_alph), alph_name(target_alph),
          alph_char(source_alph, i_col_from));
      return NULL; // placate compiler
    }
    bitset[i_col_to / 32] |= (1 << (i_col_to % 32));
    // copy the column over
    for (i_row = 0; i_row < num_rows; i_row++) {
      set_matrix_cell(i_row, i_col_to, get_matrix_cell(i_row, i_col_from, in), out);
    }
  }
  return out;
}

/***********************************************************************
 * Shuffle the positions of the motif
 ***********************************************************************/
void shuffle_motif
  (MOTIF_T* motif, mt_state* prng)
{
  int i, j, w;
  int *permute;

  w = get_motif_length(motif);
  permute = mm_malloc(sizeof(int) * w);
  // setup a permute list using an inside-out Fisher-Yates shuffle
  // http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_.22inside-out.22_algorithm
  for (i = 0;  i < w; i++) {
    j = rand_uint32(prng, i);
    if (j != i) {
      permute[i] = permute[j];
    }
    // in this case we are using a source which is a ordered list of numbers,
    // hence it is equal to i
    permute[j] = i; //permute[j] = source[i];
  }
  // permute the frequencies
  permute_matrix(motif->freqs, false, permute, w);
  permute_matrix(motif->scores, false, permute, w);
  // TODO recalculate other derived values?
  free(permute);
}

/***********************************************************************
 * Takes a matrix of letter probabilities and converts them into meme
 * score.
 *
 * The site count must be larger than zero and the pseudo count must be
 * positive.
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
   double pseudo_count) 
{
  int asize, length;
  double freq, score, total_count, counts, bg_freq;
  MATRIX_T *scores;
  int row, col;

  assert(alph != NULL);
  assert(freqs != NULL);
  assert(bg != NULL);
  assert(site_count > 0);
  assert(pseudo_count >= 0);

  length = get_num_rows(freqs);
  asize = alph_size_core(alph);

  scores = allocate_matrix(length, asize);
  total_count = site_count + pseudo_count;

  for (col = 0; col < asize; ++col) {
    bg_freq = get_array_item(col, bg);
    for (row = 0; row < length; ++row) {
      freq = get_matrix_cell(row, col, freqs);
      // apply a pseudo count
      freq = ((pseudo_count * bg_freq) + (freq * site_count)) / total_count;
      // if the background is correct this shouldn't happen
      if (freq <= 0) freq = 0.0000005;
      // the user might provide a background or motif file with 0 bg freqs
      if (bg_freq <= 0) bg_freq = 0.0000005;
      // convert to a score
      score = (log(freq / bg_freq) / log(2)) * 100;
      set_matrix_cell(row, col, score, scores);
    }
  }
  return scores;
}

/***********************************************************************
 * Takes a matrix of meme scores and converts them into letter 
 * probabilities.
 *
 * The site count must be larger than zero and the pseudo count must be
 * positive.
 *
 * The probablility can be got by:
 * p = (2 ^ (s / 100)) * bg
 *
 ***********************************************************************/
MATRIX_T* convert_scores_into_freqs
  (ALPH_T* alph,
   MATRIX_T *scores,
   ARRAY_T *bg,
   int site_count,
   double pseudo_count)
{
  int asize, length;
  double freq, score, total_count, counts, bg_freq;
  MATRIX_T *freqs;
  int row, col;

  assert(alph != NULL);
  assert(scores != NULL);
  assert(bg != NULL);
  assert(site_count > 0);
  assert(pseudo_count >= 0);

  length = get_num_rows(scores);
  asize = alph_size_core(alph);

  freqs = allocate_matrix(length, asize);
  total_count = site_count + pseudo_count;

  for (col = 0; col < asize; ++col) {
    bg_freq = get_array_item(col, bg);
    for (row = 0; row < length; ++row) {
      score = get_matrix_cell(row, col, scores);
      // convert to a probability
      freq = pow(2.0, score / 100.0) * bg_freq;
      // remove the pseudo count
      freq = ((freq * total_count) - (bg_freq * pseudo_count)) / site_count;
      if (freq < 0) freq = 0;
      else if (freq > 1) freq = 1;
      set_matrix_cell(row, col, freq, freqs);
    }
  }
  for (row = 0; row < length; ++row) {
    normalize_subarray(0, asize, 0.0, get_matrix_row(row, freqs));
  }

  return freqs;
}

/***********************************************************************
 * Turn a given motif into its own reverse complement.
 ***********************************************************************/
void reverse_complement_motif
  (MOTIF_T* a_motif)
{
  int i, temp_trim;
  ARRAY_T* left_freqs;
  ARRAY_T* right_freqs;

  assert(alph_has_complement(a_motif->alph));

  if (a_motif->freqs) {
    // Consider each row (position) in the motif.
    for (i = 0; i < (int)((a_motif->length + 1) / 2); i++) {
      left_freqs = get_matrix_row(i, a_motif->freqs);
      right_freqs = get_matrix_row(a_motif->length - (i + 1), a_motif->freqs);
      complement_swap_freqs(a_motif->alph, left_freqs, right_freqs);
    }
  }
  if (a_motif->scores) {
    // Consider each row (position) in the motif.
    for (i = 0; i < (int)((a_motif->length + 1) / 2); i++) {
      left_freqs = get_matrix_row(i, a_motif->scores);
      right_freqs = get_matrix_row(a_motif->length - (i + 1), a_motif->scores);
      complement_swap_freqs(a_motif->alph, left_freqs, right_freqs);
    }
  }
  //swap the trimming variables
  temp_trim = a_motif->trim_left;
  a_motif->trim_left = a_motif->trim_right;
  a_motif->trim_right = temp_trim;
  //swap the strand indicator
  //this assumes a ? is equalivant to +
  if (get_motif_strand(a_motif) == '-') {
    set_motif_strand('+', a_motif);
  } else {
    set_motif_strand('-', a_motif);
  }
}

/***********************************************************************
 * Apply a pseudocount to the motif pspm.
 ***********************************************************************/
void apply_pseudocount_to_motif
  (MOTIF_T* motif, ARRAY_T *background, double pseudocount)
{
  int pos, letter, len, asize, sites;
  double prob, count, total;
  ARRAY_T *temp;

  // no point in doing work when it makes no difference
  if (pseudocount == 0) return;
  assert(pseudocount > 0);
  // motif dimensions
  asize = alph_size_core(motif->alph);
  len = motif->length;
  // create a uniform background if none is given
  temp = NULL;
  if (background == NULL) {
    temp = get_uniform_frequencies(motif->alph, NULL);
    background = temp;
  }
  // calculate the counts
  sites = (motif->num_sites > 0 ? motif->num_sites : DEFAULT_SITE_COUNT);
  total = sites + pseudocount;
  for (pos = 0; pos < len; ++pos) {
    for (letter = 0; letter < asize; ++letter) {
      prob = get_matrix_cell(pos, letter, motif->freqs);
      count = (prob * sites) + (pseudocount * get_array_item(letter, background));
      prob = count / total;
      set_matrix_cell(pos, letter, prob, motif->freqs);
    }
    if (motif->flags & MOTIF_HAS_AMBIGS) {
      // recalculate ambiguous symbol values
      calc_ambigs(motif->alph, FALSE, get_matrix_row(pos, motif->freqs));
    }
  }
  if (temp) free_array(temp);
}

/***********************************************************************
 * Calculate the ambiguous letters from the concrete ones.
 ***********************************************************************/
void calc_motif_ambigs
  (MOTIF_T *motif)
{
  int i_row;
  resize_matrix(motif->length, alph_size_full(motif->alph), 0, motif->freqs);
  motif->flags |= MOTIF_HAS_AMBIGS;
  for (i_row = 0; i_row < motif->length; ++i_row) {
    calc_ambigs(motif->alph, FALSE, get_matrix_row(i_row, motif->freqs));
  }
}

/***********************************************************************
 * Normalize the motif's pspm
 ***********************************************************************/
void normalize_motif
  (MOTIF_T *motif, double tolerance)
{
  int i_row, asize;
  asize = alph_size_core(motif->alph);
  for (i_row = 0; i_row < motif->length; ++i_row) {
    normalize_subarray(0, asize, tolerance, get_matrix_row(i_row, motif->freqs));
  }
}

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
) {
  int i, len;

  len = a_motif->length;
  for (i = 0; i < len; ++i) {
    if (position_information_content(a_motif, i) >= threshold_bits) break;
  }
  a_motif->trim_left = i;
  if (i == len) {
    a_motif->trim_right = 0;
    return;
  }
  for (i = len-1; i >= 0; --i) {
    if (position_information_content(a_motif, i) >= threshold_bits) break;
  }
  a_motif->trim_right = len - i - 1;
}

/***********************************************************************
 * Compute the complexity of a motif as a number between 0 and 1.
 *
 * Motif complexity is the average K-L distance between the "motif
 * background distribution" and each column of the motif.  The motif
 * background is just the average distribution of all the columns.  The
 * K-L distance, which measures the difference between two
 * distributions, is the same as the information content:
 *
 *  \sum_i p_i log(p_i/f_i)
 *
 * This value increases with increasing complexity.
 ***********************************************************************/
double compute_motif_complexity
  (MOTIF_T* a_motif)
{
  double return_value;
  ARRAY_T* motif_background;  // Mean emission distribution.
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;

  num_cols = alph_size_core(a_motif->alph);
  num_rows = a_motif->length;

  // Compute the mean emission distribution.
  motif_background = get_matrix_col_sums(a_motif->freqs);
  scalar_mult(1.0 / (double)num_rows, motif_background);

  // Compute the K-L distance w.r.t. the background.
  return_value = 0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    ARRAY_T* this_emission = get_matrix_row(i_row, a_motif->freqs);
    for (i_col = 0; i_col < num_cols; i_col++) {
      ATYPE this_item = get_array_item(i_col, this_emission);
      ATYPE background_item = get_array_item(i_col, motif_background);

      // Use two logs to avoid handling divide-by-zero as a special case.
      return_value += this_item 
        * (my_log(this_item) - my_log(background_item));
    }
  }

  free_array(motif_background);
  return(return_value / (double)num_rows);
}

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
   MOTIF_T*  a_motif)
{
  // Make sure the given threshold is in the right range.
  if ((threshold < 0.0) || (threshold > 100.0)) {
    die(
      "Information threshold (%g) must be a percentage between 0 and 100.\n",
            threshold
    );
  }

  // Get the dimensions of the motif.
  int num_cols = alph_size_core(a_motif->alph);
  int num_rows = a_motif->length;

  // Compute and store the information content for each row
  // and the total information content for the motif.
  ATYPE total_information_content = 0.0;
  ARRAY_T* information_content = allocate_array(num_rows);
  int i_row;
  int i_col;
  for (i_row = 0; i_row < num_rows; i_row++) {
    ATYPE row_content = 0.0;
    ARRAY_T* this_emission = get_matrix_row(i_row, a_motif->freqs);
    for (i_col = 0; i_col < num_cols; i_col++) {
      ATYPE this_item = get_array_item(i_col, this_emission);
      ATYPE background_item = get_array_item(i_col, background);

      // Use two logs to avoid handling divide-by-zero as a special case.
      ATYPE partial_row_content = 
        this_item * (my_log(this_item) - my_log(background_item));

      row_content += partial_row_content;
      total_information_content += partial_row_content;

    }
    set_array_item(i_row, row_content, information_content);
  }

  // Search for the target position.
  int return_value = -1;
  ATYPE cumulative_content = 0.0;
  ATYPE percent = 0.0;
  if (from_start) {
    // Search from start for IC exceeding threshold.
    for (i_row = 0; i_row < num_rows; i_row++) {
      cumulative_content += get_array_item(i_row, information_content);
      percent = 100 *  cumulative_content / total_information_content;
      if (percent >= threshold) {
        return_value = i_row;
        break;
      }
    }
  }
  else {
    // Search from end for IC exceeding threshold.
    for (i_row = num_rows - 1; i_row >= 0; i_row--) {
      cumulative_content += get_array_item(i_row, information_content);
      percent = 100 *  cumulative_content / total_information_content;
      if (percent >= threshold) {
        return_value = i_row;
        break;
      }
    }
  }

  if (return_value == -1) {
    die(
      "Can't find a position that accounts for %g of information content.",
      threshold
    );
  }
  free_array(information_content);
  return(return_value);
}


/***********************************************************************
 * Returns the string that is the best possible match to the given motif.
 * Caller is responsible for freeing string.
 ***********************************************************************/
char *get_best_possible_match(MOTIF_T *motif) {
  int mpos, apos, asize; 
  char *match_string;
  int size;

  asize = alph_size_core(motif->alph);
  
  assert(motif != NULL);
  assert(motif->freqs != NULL);
  assert(motif->length == motif->freqs->num_rows);
  size = (motif->flags & MOTIF_HAS_AMBIGS ? 
      alph_size_full(motif->alph) : alph_size_core(motif->alph));
  assert(size == motif->freqs->num_cols); 

  match_string = mm_malloc(sizeof(char) * (motif->length + 1));

  // Find the higest scoring character at each position in the motif.
  for(mpos = 0; mpos < motif->length; ++mpos) {
    ARRAY_T *row = motif->freqs->rows[mpos];
    double max_v = row->items[0];
    int max_i = 0;
    for(apos = 1; apos < asize; ++apos) {
     if (row->items[apos] >= max_v) {
        max_i = apos;
        max_v = row->items[apos];
     }
    }
    match_string[mpos] = alph_char(motif->alph, max_i);
  }

  //  Add null termination
  match_string[motif->length] = '\0';

  return match_string;
}

/***********************************************************************
 * Duplicates the motif
 ***********************************************************************/
MOTIF_T* duplicate_motif
  (MOTIF_T *motif)
{
  MOTIF_T *motif_copy;
  motif_copy = mm_malloc(sizeof(MOTIF_T));
  copy_motif(motif, motif_copy);
  return motif_copy;
}

/***********************************************************************
 * Duplicates and reverse complements the motif
 ***********************************************************************/
MOTIF_T* dup_rc_motif
  (MOTIF_T *motif)
{
  MOTIF_T *rc_motif;
  rc_motif = mm_malloc(sizeof(MOTIF_T));
  copy_motif(motif, rc_motif);
  reverse_complement_motif(rc_motif);
  return rc_motif;
}

/***********************************************************************
 * Free dynamic memory used by one motif assuming the structure itself
 * does not need to be freed. 
 ***********************************************************************/
void free_motif
  (MOTIF_T *a_motif)
{
  /* Don't bother with empty motifs. */
  if (a_motif == NULL) 
    return;

  // Free dynamic members.
  myfree(a_motif->consensus);
  alph_release(a_motif->alph); 	// release the alphabet
  free_matrix(a_motif->freqs);
  free_matrix(a_motif->scores);
  myfree(a_motif->url);

  // Reset all member values.
  memset(a_motif, 0, sizeof(MOTIF_T));
}

/***********************************************************************
 * Free dynamic memory used by a given motif and free the structure.
 * To be useable by collections it takes a void * but expects
 * a MOTIF_T *.
 ***********************************************************************/
void destroy_motif
  (void * a_motif)
{
  free_motif((MOTIF_T*)a_motif);
  free((MOTIF_T*)a_motif);
}

/***********************************************************************
 * Convert a list of motifs into an array of motifs with a count.
 * This is intended to allow backwards compatibility with the older
 * version.
 ***********************************************************************/
void motif_list_to_array(ARRAYLST_T *motif_list, MOTIF_T **motif_array, int *num) {
  int count, i;
  MOTIF_T *motifs;
  count = arraylst_size(motif_list);
  motifs = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T) * count);
  for (i = 0; i < count; ++i) {
    copy_motif((MOTIF_T*)arraylst_get(i, motif_list), motifs+i);
  }
  *motif_array = motifs;
  *num = count;
}

/***********************************************************************
 * Convert a tree of motifs into an array of motifs with a count.
 * This is intended to allow backwards compatibility with the older
 * version.
 ***********************************************************************/
void motif_tree_to_array(RBTREE_T *motif_tree, MOTIF_T **motif_array, int *num) {
  int count, i;
  MOTIF_T *motifs;
  RBNODE_T *node;

  count = rbtree_size(motif_tree);
  motifs = mm_malloc(sizeof(MOTIF_T) * count);
  for (i = 0, node = rbtree_first(motif_tree); node != NULL; i++, node = rbtree_next(node)) {
    copy_motif((MOTIF_T*)rbtree_value(node), motifs+i);
  }
  *motif_array = motifs;
  *num = count;
}

/***********************************************************************
 * Get the motif at the selected index in the array. This is needed for
 * the older version which used arrays of motifs.
 ***********************************************************************/
MOTIF_T* motif_at(MOTIF_T *motif_array, int index) {
  return motif_array+index;
}

/***********************************************************************
 * Free an array of motifs
 ***********************************************************************/
void free_motif_array(MOTIF_T *motif_array, int num) {
  int i;
  for (i = 0; i < num; ++i) {
    free_motif(motif_array+i);
  }
  free(motif_array);
}

/***********************************************************************
 * free an array list and the contained motifs
 ***********************************************************************/
void free_motifs(
  ARRAYLST_T *motif_list
) {
  arraylst_destroy(destroy_motif, motif_list);
}

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 *
 * The reverse complement is always listed directly after the original.
 ***********************************************************************/
void add_reverse_complements
  (ARRAYLST_T* motifs)
{
  int i, count;

  count = arraylst_size(motifs);
  //double the array size
  arraylst_preallocate(count*2, motifs);
  arraylst_add_n(count, NULL, motifs);
  //move and reverse complement the original motifs
  for (i = count-1; i >= 0; --i) {
    //move the motif to its new position
    arraylst_swap(i, 2*i, motifs);
    //get the motif and the one that will become the reverse complement
    MOTIF_T *orig = arraylst_get(2*i, motifs);
    if (get_motif_strand(orig) == '?') set_motif_strand('+', orig);
    //copy and reverse complement the motif to the position after
    MOTIF_T *rc = dup_rc_motif(orig);
    arraylst_set(2*i + 1, rc, motifs);
  }
  assert(arraylst_size(motifs) == (2 * count));
}

/***********************************************************************
*
* Dump the frequencies of the motif to the output.
*
 ***********************************************************************/
void dump_motif_freqs(FILE *out, MOTIF_T* m){
  int i, j;
  ALPH_T *alph;
  MATRIX_T *freqs;
  alph = get_motif_alph(m);
  freqs = get_motif_freqs(m);
  fprintf(out, "%s", get_motif_st_id(m));
  for (i = 0; i < alph_size_core(alph); i++) {
    fprintf(out, "\t%c", alph_char(alph, i));
  }
  fprintf(out, "\n");
  for (i = 0; i < get_num_rows(freqs); i++) {
    fprintf(out, "%s", get_motif_st_id(m));
    for (j = 0; j < alph_size_core(alph); j++) {
      fprintf(out, "\t%.4f", get_matrix_cell(i, j, freqs));
    }
    fprintf(out, "\n");
  }
}

/***********************************************************************
 *
 * Check that all alphabets are the same or convertible in a set of
 * motif files
 *
 ***********************************************************************/
void read_motif_alphabets(ARRAYLST_T* motif_sources, bool xalph, ALPH_T** alph) {
  int i;
  char *db_source;
  ALPH_T *db_alph;
  for (i = 0; i < arraylst_size(motif_sources); i++) {
    char* db_source = (char*) arraylst_get(i, motif_sources);
    // open the motif file for reading
    MREAD_T *mread = mread_create(db_source, OPEN_MFILE);
    // get the alphabet from the motif file
    if (mread_has_motif(mread)) {
      db_alph = mread_get_alphabet(mread);
      if (*alph != NULL) {
        if (!alph_equal(*alph, db_alph)) {
          if (xalph) {
            switch(alph_core_subset(db_alph, *alph)) {
              case 0:
                die("The motifs in \"%s\" have the %s alphabet which is not "
                    "a subset of the %s alphabet.", db_source, 
                    alph_name(db_alph), alph_name(*alph));
                break;
              case -1:
                fprintf(stderr, "Warning: the alphabet expansion from %s to %s"
                    " requires changing complementing rules.", 
                    alph_name(db_alph), alph_name(*alph));
                break;
            }
            // since we're not applying pseudocounts the background doesn't matter
            // so we use the uniform background
            mread_set_conversion(mread, *alph, NULL);
          } else {
            die("The motifs in \"%s\" have the %s alphabet which is not "
                "the same as the expected %s alphabet.", db_source, 
                alph_name(db_alph), alph_name(*alph));
          }
        }
      } else {
        *alph = alph_hold(db_alph);
      }
    }
    // clean up motif reader
    mread_destroy(mread);
  }
  if (*alph == NULL) {
    die("No alphabet found in motif file(s).\n");
  }
} // read_motif_alphabets

/**************************************************************************
 * Compares two AP_T and orders largest probability to smallest probability.
 **************************************************************************/
int ap_cmp(const void *p1, const void *p2) {
  AP_T *a1, *a2;
  a1 = (AP_T*)p1;
  a2 = (AP_T*)p2;
  if (a1->prob == a2->prob) {
    return (a1->idx > a2->idx) - (a1->idx < a2->idx);
  }
  return (a1->prob < a2->prob ? 1 : -1);
}

/**************************************************************************
 * Compares two uint8_t and orders smallest to largest.
 **************************************************************************/
int idx_cmp(const void *p1, const void *p2) {
  return (*((uint8_t*)p1) > *((uint8_t*)p2)) - (*((uint8_t*)p1) < *((uint8_t*)p2));
}

/**************************************************************************
 * Attempts to make a reasonable consensus representation of a motif.
 **************************************************************************/
void motif2consensus(MOTIF_T* motif, STR_T* consensus, bool single_letter) {
  MATRIX_T *freqs;
  ARRAY_T *row;
  AP_T *labeled_row;
  uint8_t *idx_row;
  ALPH_T *alph;
  int i, j, k, ncomprise;
  bool match;
  alph = get_motif_alph(motif);
  labeled_row = mm_malloc(sizeof(AP_T) * alph_size_core(alph));
  idx_row = mm_malloc(sizeof(uint8_t) * alph_size_core(alph));
  freqs = get_motif_freqs(motif);
  for (i = 0; i < get_motif_length(motif); i++) {
    row = get_matrix_row(i, freqs);
    // create a sortable row where we can recover the labels
    for (j = 0; j < alph_size_core(alph); j++) {
      labeled_row[j].idx = j;
      labeled_row[j].prob = get_array_item(j, row);
    }
    // sort the row by the frequency, largest to smallest
    qsort(labeled_row, alph_size_core(alph), sizeof(AP_T), ap_cmp);
    // pick any symbols within 0.5 of the best
    idx_row[0] = labeled_row[0].idx;
    for (ncomprise = 1; ncomprise < alph_size_core(alph); ncomprise++) {
      if (labeled_row[ncomprise].prob < (labeled_row[0].prob * 0.5)) break;
      idx_row[ncomprise] = labeled_row[ncomprise].idx;
    }
    // if nothing else is anywhere near as good as the best one print it
    if (ncomprise == 1) {
      str_appendf(consensus, "%c", alph_char(alph, labeled_row[0].idx));
      continue;
    }
    // otherwise try to find an ambiguous symbol that exactly matches
    qsort(idx_row, ncomprise, sizeof(uint8_t), idx_cmp);
    match = false;
    for (j = alph_size_core(alph); j < alph_size_full(alph); j++) {
      if (alph_ncomprise(alph, j) > ncomprise) continue;
      else if (alph_ncomprise(alph, j) < ncomprise) break;
      for (k = 0; k < ncomprise; k++) {
        if (alph_comprise(alph, j, k) != idx_row[k]) break;
      }
      if (k == ncomprise) {
        str_appendf(consensus, "%c", alph_char(alph, j));
        match = true;
        break;
      }
    }
    // when we don't find a good match
    if (!match) {
      if (ncomprise * 2 > alph_size_core(alph)) {
        // when there are lots of symbols print the wildcard
        str_appendf(consensus, "%c", alph_wildcard(alph));
      } else if (single_letter) {
        // not too many letters, print the most common one
        str_appendf(consensus, "%c", alph_char(alph, labeled_row[0].idx));
      } else {
        // otherwise print up to 3 of the best matching symbols
        str_append(consensus, "[", 1);
        for (j = 0; j < ncomprise; j++) {
          if (j >= 3) break;
          str_appendf(consensus, "%c", alph_char(alph, labeled_row[j].idx));
        }
        str_append(consensus, "]", 1);
      }
    }
  }
}

