/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

/*      
        alphtype <alphabet>

        <alphabet>      MAST alphabet string

        Prints "DNA" or "PROTEIN" if the alphabet is one of the
        recognized ones.  Prints an error message to standard error
        otherwise.
*/
        
#include "alphabet.h"

static inline __attribute__((always_inline)) bool all_in_alph(ALPH_T *alph, const char *symbols) {
  const char *s;
  for (s = symbols; *s; s++) {
    if (!alph_is_known(alph, *s)) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  char *symbols;
  ALPH_T *dna_alph, *protein_alph;

  if (argc < 2) {
    printf("Usage:\n\talphtype <alphabet>\n");
    return 0;
  }

  symbols = argv[1];
  dna_alph = alph_dna();
  protein_alph = alph_protein();

  if (all_in_alph(dna_alph, symbols)) {
    printf("DNA\n");
  } else if (all_in_alph(protein_alph, symbols)) {
    printf("PROTEIN\n");
  } else {
    return 1;
  }
  alph_release(dna_alph);
  alph_protein(protein_alph);
  return 0;
}
  
