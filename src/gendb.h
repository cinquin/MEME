#ifndef GENDB_H
#define GENDB_H

#include "mtwist.h"
#include "seq.h"

/*
 * Generate sequences given an alphabet, a background and the frequency
 * of ambiguous characters.
 */
SEQ_T** gendb(
  FILE *out, // Output stream; return output if null.
  mt_state *prng, // pseudo-random number generator state
  ALPH_T *alph, // alphabet
  ARRAY_T *bg, // background model
  double ambig, // frequency of ambiguous characters
  int nseqs, // Number of sequences to generate
  int min, // Shortest possible sequence
  int max // Longest possible sequence
);

#endif
