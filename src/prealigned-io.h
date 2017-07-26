#include "seq.h"
#include "seq-reader-from-fasta.h"

#define MAX_SEQ 250000000   // Default maximum sequence length

/****************************************************************************
 * Reads a single prealigned sequence int sequence line. Throws an error
 * if required width is not 0 and sequence length not equal to required width.
 ****************************************************************************/
void parse_prealigned_bodyline(
  char** sequence,
  char* line,
  int required_width);

/****************************************************************************
 * Read all the sequences from a line-delimited file
 ****************************************************************************/
void read_many_line_delimited_sequences
(ALPH_T*    alph,
 FILE*      fasta_file,
 int        max_seq, // Maximum sequence length.
 int*       num_seqs,
 int        required_width,
 SEQ_T***   sequences);