#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "utils.h"
#include "alphabet.h"
#include "io.h"
#include "seq-reader-from-fasta.h"
#include "prior-reader-from-psp.h"
#include "seq.h"
#include "prealigned-io.h"

/****************************************************************************
 * Reads a single prealigned sequence int sequence line. Throws an error
 * if required width is not 0 and sequence length not equal to required width.
 ****************************************************************************/
void parse_prealigned_bodyline(char** sequence, char* line, int required_width) {
  if (!line) {
    *sequence = NULL;
  } else {
    if (required_width != 0 && strlen(line) != required_width) {
      die("Found a line with width != specified options width: %s. If specifying in phosphorylation sequences in prealigned or fasta, we assume that all sequences are centralized around a mod. Strlen: %d, but specified width is: %d", line, strlen(line), required_width);
    }
    *sequence = strdup(line);
  }
}

/****************************************************************************
 * Read all the sequences from a line-delimited file at once.
 Multiple files can be appended by calling this more than once.
 ****************************************************************************/
#define NUM_ALLOC 1000 /* Allocate how many sequences at once? */
void read_many_line_delimited_sequences
(ALPH_T*    alph,
 FILE*      file,
 int        max_seq, // Maximum sequence length.
 int*       num_seqs,	// Must be set to zero unless appending to sequences.
 int        required_width, // is there a required width of sequences? 0 if no.
 SEQ_T***   sequences)
{
  int i_seq;         /* Index of the current sequence. */
  int num_allocated; /* Number of pointers currently allocated. */
  
  /* Allocate initial memory. */
  num_allocated = *num_seqs + NUM_ALLOC;
  if (*num_seqs > 0) {
    *sequences = (SEQ_T**)mm_realloc(*sequences, sizeof(SEQ_T*) * num_allocated);
  } else {
    *sequences = (SEQ_T**)mm_malloc(sizeof(SEQ_T*) * num_allocated);
  }
  
  char* line = NULL;
  size_t len = 0;
  ssize_t read;
  
  /* Read the sequences one by one. */
  i_seq = *num_seqs;
  while ((read = getline(&line, &len, file)) != -1) {
    line = strtok(line, "\n");
    
    char* sequence = NULL;
    parse_prealigned_bodyline(&sequence, line, required_width);
    (*sequences)[i_seq] = allocate_seq("tempname", "tempdescription", 0, sequence);
    i_seq++;
    
    /* Allocate more space, if need be. */
    if (i_seq >= num_allocated) {
      num_allocated += NUM_ALLOC;
      *sequences = (SEQ_T**)mm_realloc(*sequences, sizeof(SEQ_T*) * num_allocated);
    }
  }
  
  /* Record the total number of sequences. */
  *num_seqs = i_seq;
  
  /* Complain if nothing was read. */
  if (i_seq == 0) {
    die("Failed to read a single sequence from the given FASTA file.\n");
  }
  
}
