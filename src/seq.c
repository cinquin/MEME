/****************************************************************************
 * FILE: seq.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 06/24/2002
 * PROJECT: MHMM
 * DESCRIPTION: Biological sequences.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include "alphabet.h"
#include "utils.h"
#include "seq.h"

#define MAX_SEQ_NAME 100     // Longest sequence ID.
#define MAX_SEQ_COMMENT 128 // Longest comment.

// maximum count to keep track of when calculating a background
static const long BG_CALC_CHUNK = LONG_MAX;

// Instantiate the SEQ_T type.
struct seq {
  char  name[MAX_SEQ_NAME + 1];     // Sequence ID.
  char  desc[MAX_SEQ_COMMENT + 1];  // Sequence description.
  unsigned int length;              // Sequence length.
  unsigned int offset;              // Offset from the start of complete sequence.
  unsigned int starting_coord;      // Starting coord of the start of complete sequence.
  float weight;                 // Sequence weight.
  char* sequence;               // The actual sequence.
  int8_t* isequence;           // indexed sequence
  int*  intseq;                 // The sequence in integer format.
  unsigned int num_priors;      // Number of priors
  double* priors;               // The priors corresponding to the sequence.
  int*  gc;                     // Cumulative counts of the second complementary pair (G & C in DNA); note: 2Gb size limit.
  double total_gc;              // Total frequency of the second complementary pair (G & C in DNA) in sequence; Complementary alphabet with exactly 2 pairs only
  BOOLEAN_T is_complete;        // Is this the end of the sequence?
};


/****************************************************************************
 * Allocate one sequence object.
 ****************************************************************************/
SEQ_T* allocate_seq
  (char* name,
   char* description,
   unsigned int offset,
   char* sequence)
{
  SEQ_T* new_sequence;

  // Allocate the sequence object.
  new_sequence = (SEQ_T*)mm_malloc(sizeof(SEQ_T));

  // Fill in the fields.
  new_sequence->priors = NULL;
  new_sequence->length = 0;
  new_sequence->offset = offset;
  new_sequence->starting_coord = offset;
  new_sequence->sequence = NULL;
  new_sequence->isequence = NULL;
  new_sequence->weight = 1.0;
  new_sequence->intseq = NULL;
  new_sequence->gc = NULL;
  new_sequence->num_priors = 0;
  new_sequence->is_complete = FALSE;

  // Store the name, truncating if necessary.
  new_sequence->name[0] = 0;
  if (name != NULL) {
    strncpy(new_sequence->name, name, MAX_SEQ_NAME);
    new_sequence->name[MAX_SEQ_NAME] = '\0';
    if (strlen(new_sequence->name) != strlen(name)) {
      fprintf(stderr, "Warning: truncating sequence ID %s to %s.\n",
          name, new_sequence->name);
    }
  }

  // Store the description, truncating if necessary.
  new_sequence->desc[0] = 0;
  if (description != NULL) {
    new_sequence->desc[0] = '\0';
    if (description) {
      strncpy(new_sequence->desc, description, MAX_SEQ_COMMENT);
      new_sequence->desc[MAX_SEQ_COMMENT] = '\0';
    }
  }

  // Store the sequence.
  if (sequence != NULL) {
    copy_string(&(new_sequence->sequence), sequence);
    new_sequence->length = strlen(sequence);
  }

  return(new_sequence);
}

/****************************************************************************
 * Get and set various fields.
 ****************************************************************************/
char* get_seq_name
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->name);
}

char* get_seq_description
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->desc);
}

unsigned int get_seq_length
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->length);
}

unsigned int get_seq_offset
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->offset);
}

void set_seq_offset
  (unsigned int offset, 
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  a_sequence->offset = offset;
}

unsigned int get_seq_starting_coord
  (SEQ_T* a_sequence)
{
  return a_sequence->starting_coord;
}

void set_seq_starting_coord
  (unsigned int start,
  SEQ_T* a_sequence)
{
  a_sequence->starting_coord = start;
}

float get_seq_weight
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->weight);
}

void set_seq_weight
  (float  weight,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  a_sequence->weight = weight;
}

unsigned char get_seq_char
  (int index,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert(a_sequence->sequence != NULL);
  assert(index >= 0);
  assert(index <= a_sequence->length);
  return(a_sequence->sequence[index]);
}

void set_seq_char
  (int    index,
   char   a_char,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert(a_sequence->sequence != NULL);
  assert(index >= 0);
  assert(index <= a_sequence->length);
  a_sequence->sequence[index] = a_char;
}

int get_seq_int
  (int index,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert(a_sequence->sequence != NULL);
  assert(index >= 0);
  assert(index < a_sequence->length);
  return(a_sequence->intseq[index]);
}

bool has_seq_gc(SEQ_T* a_sequence) {
  return (a_sequence != NULL && a_sequence->gc != NULL);
}

int get_seq_gc
  (int index,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert(a_sequence->sequence != NULL);
  assert(index >= 0);
  assert(index < a_sequence->length);
  return(a_sequence->gc[index]);
}

/**************************************************************************
 * Return a pointer to the raw sequence data in a SEQ_T object.
 **************************************************************************/
char* get_raw_sequence
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->sequence);
}

char* get_raw_subsequence
  (int start, int stop, SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert((stop - start) >= 0);
  char *sequence = get_raw_sequence(a_sequence);
  char *subsequence = mm_malloc((stop - start + 2) * sizeof(char));
  strncpy(subsequence, sequence + start, stop - start + 1);
  subsequence[stop - start + 1] = 0;
  return(subsequence);
}

void set_raw_sequence(
  char *raw_sequence,
  BOOLEAN_T is_complete,
  SEQ_T* a_sequence
) {
  assert(a_sequence != NULL);
  a_sequence->sequence = raw_sequence;
  a_sequence->length = strlen(raw_sequence);
  a_sequence->is_complete = is_complete;

}

int8_t* get_isequence(SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  assert(a_sequence->isequence != NULL);
  return(a_sequence->isequence);
}

int* get_int_sequence
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->intseq);
}


int* get_gc_sequence
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->gc);
}

double get_total_gc_sequence
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->total_gc);
}

void set_total_gc_sequence
  (SEQ_T* a_sequence, double gc)
{
  assert(a_sequence != NULL);
  a_sequence->total_gc = gc;
}

/**************************************************************************
 * Return the number of priors currently stored in the SEQ_T object
 **************************************************************************/
unsigned int get_seq_num_priors
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->num_priors);
}

/**************************************************************************
 * Return a pointer to the priors associated with a SEQ_T object.
 * It may be a NULL pointer. If it is not NULL the number of priors
 * is the length member of the SEQ_T object.
 **************************************************************************/
double* get_seq_priors
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->priors);
}

/**************************************************************************
 * Set the priors data for a SEQ_T object.
 **************************************************************************/
void set_seq_priors
  (double *priors,
  SEQ_T* a_sequence) 
{
  assert(a_sequence != NULL);
  a_sequence->priors = priors;
}


/**************************************************************************
 * Copy a sequence object.  Memory must be freed by caller.
 **************************************************************************/
SEQ_T* copy_sequence
  (SEQ_T* source_sequence)
{
  // Allocate the sequence object.
  SEQ_T* new_sequence;

  new_sequence = (SEQ_T*)mm_malloc(sizeof(SEQ_T));

  // Store the name
  // moving from one buffer to another of identical size so no truncation
  strncpy(new_sequence->name, source_sequence->name, MAX_SEQ_NAME);
  new_sequence->name[MAX_SEQ_NAME] = '\0';

  // Store the description
  // moving from one buffer to another of identical size so no truncation
  strncpy(new_sequence->desc, source_sequence->desc, MAX_SEQ_COMMENT);
  new_sequence->desc[MAX_SEQ_COMMENT] = '\0';

  // Fill in the always set fields
  new_sequence->starting_coord = source_sequence->starting_coord;
  new_sequence->length = source_sequence->length;
  new_sequence->offset = source_sequence->offset;
  new_sequence->weight = source_sequence->weight;
  new_sequence->is_complete = source_sequence->is_complete;
  new_sequence->num_priors = source_sequence->num_priors;

  // copy over the sequence (note that either this or the isequence will be present)
  if (source_sequence->sequence != NULL) {
    new_sequence->sequence = mm_malloc(sizeof(char) * (source_sequence->length + 1));
    memcpy(new_sequence->sequence, source_sequence->sequence, (sizeof(char) * source_sequence->length));
    new_sequence->sequence[source_sequence->length] = '\0';
  } else {
    new_sequence->sequence = NULL;
  }
  // copy over the isequence (note that either this or the sequence will be present)
  if (source_sequence->isequence != NULL) {
    new_sequence->isequence = mm_malloc(sizeof(int8_t) * source_sequence->length);
    memcpy(new_sequence->isequence, source_sequence->isequence, (sizeof(int8_t) * source_sequence->length));
  } else {
    new_sequence->isequence = NULL;
  }
  // copy over the int sequence (this may be present, used by mhmm code, 
  // required to be int so that hashing of multiple characters can be done)
  if (source_sequence->intseq != NULL) {
    new_sequence->intseq = (int*)mm_malloc(sizeof(int) * source_sequence->length);
    memcpy(new_sequence->intseq, source_sequence->intseq, (sizeof(int) * source_sequence->length));
  } else {
    new_sequence->intseq = NULL;
  }
  // copy over the gc count (this may be present, used by mhmm code)
  if (source_sequence->gc != NULL) {
    new_sequence->gc = (int*)mm_malloc(sizeof(int) * source_sequence->length);
    memcpy(new_sequence->gc, source_sequence->gc, (sizeof(int) * source_sequence->length));
  } else {
    source_sequence->gc = NULL;
  }
  // copy over the priors
  if (source_sequence->priors && source_sequence->num_priors > 0) {
    new_sequence->priors = mm_malloc(sizeof(double) * source_sequence->num_priors);
    memcpy(new_sequence->priors, source_sequence->priors, (sizeof(double) * source_sequence->num_priors));
  } else {
    new_sequence->priors = NULL;
  }
  return(new_sequence);
}

/**************************************************************************
 * Convert a sequence to an indexed version. This is a one-way conversion
 * but you can choose to keep the source sequence if you need it
 **************************************************************************/
void index_sequence(SEQ_T* seq, ALPH_T* alph, int options) {
  char *in;
  int8_t *out;
  in = seq->sequence;
  if ((options & SEQ_KEEP) != 0 || sizeof(char) != sizeof(int8_t)) {
    out = mm_malloc(sizeof(int8_t) * seq->length);
    seq->isequence = out;
  } else {
    // use the same memory and overwrite the sequence with the indexed version
    out = (int8_t*)seq->sequence;
    seq->isequence = out;
    seq->sequence = NULL;
  }
  if ((options & SEQ_NOAMBIG) != 0) {
    for (; *in != '\0'; in++, out++) {
      *out = alph_indexc(alph, *in);
    }
  } else {
    for (; *in != '\0'; in++, out++) {
      *out = alph_index(alph, *in);
    }
  }
}


/**************************************************************************
 * Sometimes a sequence object contains only a portion of the actual
 * sequence.  This function tells whether or not the end of this
 * sequence object corresponds to the end of the actual sequence.
 **************************************************************************/
void set_complete
  (BOOLEAN_T is_complete,
   SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  a_sequence->is_complete = is_complete;
}

BOOLEAN_T is_complete
  (SEQ_T* a_sequence)
{
  assert(a_sequence != NULL);
  return(a_sequence->is_complete);
}


/***************************************************************************
 * Add or remove Xs from either side of the sequence.
 ***************************************************************************/
static void add_flanking_xs
  (SEQ_T* sequence, ALPH_T* alph)
{
  char*  new_seq = NULL;         // Pointer to copy of the sequence.

  new_seq = (char*)mm_calloc(sequence->length + 3, sizeof(char));
  strcpy(&(new_seq[1]), sequence->sequence);

  new_seq[0] = alph_wildcard(alph);
  new_seq[sequence->length + 1] = alph_wildcard(alph);
  new_seq[sequence->length + 2] = '\0';

  myfree(sequence->sequence);
  sequence->sequence = new_seq;
  sequence->length += 2;
}

void remove_flanking_xs
  (SEQ_T* sequence)
{
  char*  new_seq;         // Copy of the sequence.

  new_seq = (char*)mm_calloc(sequence->length - 1, sizeof(char));
  strncpy(new_seq, &(sequence->sequence[1]), sequence->length - 2);
  new_seq[sequence->length - 2] = '\0';

  myfree(sequence->sequence);
  sequence->sequence = new_seq;
  sequence->length -= 2;
}

/**************************************************************************
 * Prepare a sequence for recognition by
 *  - making sure it doesn't contain illegal characters,
 *  - adding flanking Xs to match START/END states, and
 *  - converting it to an integer format
 *  - computing cumulative GC counts for alphabets with 2 complementary pairs
 *  - if hard_mask is set, replace lower case characters with wildcard
 * In the integer form, each character in the sequence is replaced by
 * the index of that character in the alphabet. 
 **************************************************************************/
void prepare_sequence (
    SEQ_T* sequence, 
    ALPH_T* alph,
    BOOLEAN_T hard_mask
)
{
  int i_seq;        // Index in the sequence.
  int badchar;      // Number of characters converted.
  char wildcard;    // Wildcard character

  wildcard = alph_wildcard(alph);
  badchar = 0;

  for (i_seq = 0; i_seq < get_seq_length(sequence); i_seq++) {

    char c = sequence->sequence[i_seq];
    // Convert non-alphabetic characters to ambiguous.
    if (!alph_is_known(alph, c) || (hard_mask && islower(c)) ) {
      if (! (hard_mask && islower(c))) {
        fprintf(stderr, "%c -> %c\n", c, wildcard);
      }
      (sequence->sequence)[i_seq] = wildcard;
      badchar++;
    }
  }

  // Tell the user about the conversions.
  if (badchar) {
    fprintf(stderr, "Warning: converted %d non-alphabetic ", badchar);
    fprintf(stderr, "characters to %c in sequence %s.\n", wildcard, 
        get_seq_name(sequence));
  }

  // Add flanking X's.
  add_flanking_xs(sequence, alph);

  // Make the integer sequence.
  sequence->intseq = (int *)mm_malloc(sizeof(int) * get_seq_length(sequence));
  for (i_seq = 0; i_seq < get_seq_length(sequence); i_seq++) {
    (sequence->intseq)[i_seq]
      = alph_encode(alph, (sequence->sequence)[i_seq]);
  }

  //
  // Get cumulative GC counts.
  // Since this needs to be generic we don't use G and C but the second
  // complementary pair in the alphabet. The alphabet must only have 2
  // pairs so the problem space can be divided up on one dimension.
  //
  if (alph_size_pairs(alph) == 2) { 
    int len, x1, x2, y1, y2, a, *gc;
    char *seq;

    // first complementary pair
    x1 = 0;
    x2 = alph_complement(alph, x1);
    // second complementary pair (aka the GC pair for DNA)
    y1 = (x2 == 1 ? 2 : 1);
    y2 = alph_complement(alph, y1);

    len = get_seq_length(sequence);
    seq = sequence->sequence;
    sequence->gc = (int *)mm_malloc(sizeof(int) * get_seq_length(sequence));
    gc = sequence->gc;

    // set count at first position
    a = alph_indexc(alph, seq[0]);
    gc[0] = ((a == y1 || a == y2) ? 1 : 0);
    // set cumulative counts at rest of postitions
    for (i_seq = 1; i_seq < len; i_seq++) {
      a = alph_indexc(alph, seq[i_seq]);
      gc[i_seq] = gc[i_seq-1] + ((a == y1 || a == y2) ? 1 : 0);
    }
  }
}

/***************************************************************************
 * Remove the first N bases of a given sequence.
 ***************************************************************************/
void shift_sequence
  (int    offset,
   SEQ_T* sequence)
{
  char *new_sequence = NULL;
  double *new_priors = NULL;

  // Make a copy of the raw sequence for the overlap.
  assert(offset > 0);
  assert(offset <= sequence->length);

  memmove(
    sequence->sequence,
    sequence->sequence + offset,
    sequence->length - offset + 1
  );

  // Shift priors if needed.
  if (sequence->priors) {
    memmove(
      sequence->priors,
      sequence->priors + offset,
      (sequence->length - offset) * sizeof(double)
    );
  }

  sequence->offset += offset;

  // Free the integer version.
  myfree(sequence->intseq);
  sequence->intseq = NULL;

  // Free the GC counts.
  myfree(sequence->gc);
  sequence->gc = NULL;
}

/***************************************************************************
 * Get the maximum sequence length from a set of sequences.
 ***************************************************************************/
int get_max_seq_length
  (int     num_seqs,
   SEQ_T** sequences)
{
  int max_length;
  int this_length;
  int i_seq;

  max_length = 0;
  for (i_seq = 0; i_seq < num_seqs; i_seq++) {
    this_length = get_seq_length(sequences[i_seq]);
    if (this_length > max_length) {
      max_length = this_length;
    }
  }
  return(max_length);
}

/***************************************************************************
 * Get the maximum sequence ID length from a set of sequences.
 ***************************************************************************/
int get_max_seq_name
  (int     num_seqs,
   SEQ_T** sequences)
{
  int max_length;
  int this_length;
  int i_seq;

  max_length = 0;
  for (i_seq = 0; i_seq < num_seqs; i_seq++) {
    this_length = strlen(get_seq_name(sequences[i_seq]));
    if (this_length > max_length) {
      max_length = this_length;
    }
  }
  return(max_length);
}

/***************************************************************************
 * Set the sequence weights according to an external file.
 *
 * If the filename is "none," "internal," or NULL, then the weights are
 * set uniformly.
 ***************************************************************************/
void set_sequence_weights
  (char*    weight_filename, // Name of file containing sequence weights.
   int      num_seqs,        // Number of sequences.
   SEQ_T**  sequences)       // The sequences.
{
  ARRAY_T* weights;
  FILE *   weight_file;
  int      i_seq;

  /* Allocate the weights. */
  weights = allocate_array(num_seqs);

  /* Set uniform weights if no file was supplied. */
  if ((weight_filename == NULL) || (strcmp(weight_filename, "none") == 0) ||
      (strcmp(weight_filename, "internal") == 0)) {
    init_array(1.0, weights);
  }

  /* Read the weights from a file. */
  else {
    if (open_file(weight_filename, "r", FALSE, "weight", "weights",
		  &weight_file) == 0)
      exit(1);
    read_array(weight_file, weights);
    fclose(weight_file);

    /* Normalize the weights so they add to the total number of sequences. */
    normalize(0.0, weights);
    scalar_mult(num_seqs, weights);
  }

  /* Assign the weights to the sequences. */
  for (i_seq = 0; i_seq < num_seqs; i_seq++) {
    (sequences[i_seq])->weight = get_array_item(i_seq, weights);
  }

  /* Free the weights. */
  free_array(weights);
}


/****************************************************************************
 *  Return an array containing the frequencies in the sequence for each
 *  character of the alphabet. Characters not in the alphabet are not
 *  counted.
 ****************************************************************************/
ARRAY_T* get_sequence_freqs
  (SEQ_T* seq, ALPH_T* alph)
{
  int a_size, a_index, i;
  int total_bases;
  int *num_bases;
  int8_t *iseq;
  ARRAY_T* freqs;

  // Initialize counts for each character in the alphabet
  a_size = alph_size_core(alph);
  num_bases = mm_malloc(a_size * sizeof(int));
  for (i = 0; i < a_size; i++) {
    num_bases[i] = 0;
  }
  total_bases = 0;

  if (seq->isequence) {
    for (i = 0, iseq = seq->isequence; i < seq->length; i++, iseq++) {
      if (*iseq == -1 || *iseq >= a_size) continue;
      num_bases[*iseq]++;
      total_bases++;
    }
  } else {
    for (i = 0; i < seq->length; i++) {
      a_index = alph_index(alph, seq->sequence[i]);
      if (a_index == -1 || a_index >= a_size) continue;
      num_bases[a_index]++;
      total_bases++;
    }
  }

  freqs = allocate_array(a_size);
  for (i = 0; i < a_size; i++) {
    set_array_item(i, (double) num_bases[i] / (double) total_bases, freqs);
  }

  // Clean up the count of characters
  myfree(num_bases);

  return freqs;
}

/****************************************************************************
 * Free one sequence object.
 ****************************************************************************/
void free_seq
  (SEQ_T* a_sequence)
{
  if (a_sequence == NULL) {
    return;
  }
  myfree(a_sequence->sequence);
  myfree(a_sequence->isequence);
  myfree(a_sequence->intseq);
  myfree(a_sequence->gc);
  myfree(a_sequence->priors);
  myfree(a_sequence);
}

