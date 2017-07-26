/*
 * $Id: mast.h 1425 2006-11-01 07:14:28Z tbailey $
 * 
 * $Log$
 * Revision 1.2  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2005/07/29 18:43:40  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef MAST_H 
#define MAST_H 

#include "array.h"
#include "meme.h"
#include "macros.h"
#include "general.h"
#include "banner.h"
#include "pssm-distr.h"
#include <sys/types.h>
#include <sys/stat.h>

/*
        Macro definitions
*/

/* compute the p-value of value of n samples of an extreme value whose 
   p-value is known and return it in r
*/ 
#define EV(pvalue, n, r) \
  { \
    double __ev_tmp__ = (pvalue) * (n); /* good approximation if small */ \
    (r) = __ev_tmp__ <= 0.001 ? \
    __ev_tmp__: 1.0 - pow((1.0 - (pvalue)), (double)(n)); \
  }

/*
        Data definitions 
*/

/* scores for each position in a sequence */
typedef struct score {
  double score;                         /* the score */
  bool ic;                           /* score of inverse complement strand */
} SCORE;

/*
        Procedures
*/

SCORE **score_it(
  ALPH_T *alph,         /* database alphabet */
  XLATE_T *xlate,       /* database is different alphabet to motifs */
  bool reverse_comp, /* same effect as reverse complementing the motifs */
  LO *los[],            /* array of pointers to log-odds matrices */
  int nmotifs,          /* number of motifs */
  char *sequence,       /* sequence */
  long length           /* length of the sequence */
);

SCORE **score_sequence(
  ALPH_T *alph,         /* database alphabet */
  XLATE_T *xlate,       /* database is different alphabet to motifs */
  STYPE stype,          /* handling of different DNA strands */
  bool neg_strand, /* strand to score, makes no difference to combine mode */
  char *sequence,       /* the sequence to score */
  long length,          /* length of the sequence */
  int nmotifs,          /* number of motifs */
  LO *los[]             /* array of pointers to log-odds matrices */
);

TILING score_tile_diagram(
  FILE *mast_out,       /* output */
  ALPH_T *alph,
  XLATE_T *xlate,
  char *sequence,       /* sequence to score and tile */
  long length,          /* length of sequence */
  LO *los[],            /* array of pointers to lo matrices */
  int nmotifs,          /* number motifs read */
  STYPE stype,          /* handling of different strands */
  bool neg_strand,   /* for dna sequences, score the negative strand */
  bool best_motifs,  /* show only best motifs in diagrams */
  bool print_p,      /* print p-value in block diagram */
  double **pv,          /* p-value tables for each motif */
  double m_thresh,      /* maximum motif p-value to print */
  double w_thresh,      /* maximum motif p-value-- weak hits */
  bool use_seq_p,    /* use sequence not position p-values */
  bool hit_list,     /* create hit list instead of diagram */
  char *name            /* name of sequence; only used if hit_list=true */
);

void print_diagram(
  char *dia,                            /* motif diagram string */
  char *hdr,                            /* prefix for each line of diagram */
  FILE *file                            /* destination file */
);

double qfast(
  int n,                        /* number of random variables in product */
  double k                      /* product of random variables */
);

void free_tiling(
  TILING tiling
); 

char *create_diagram(
  XLATE_T *xlate,
  STYPE stype,                          /* treatment of strands of DNA */
  bool best_motifs,                  /* diagrams have only best motif */
  bool print_p,                      /* print p-value in block */
  double thresh,                        /* strong hit threshold */
  int nmotifs,                          /* number of motifs */
  LO *los[],                            /* array of pointers to lo matrices */
  long length,                          /* length of sample */
  bool hit_list,                     /* create hit list instead of diagram */
  char *name,                           /* name of sequence; only used if hit_list=true */
  int skip,                             // skip first "skip" positions in sequence
  int offset,                           // offset of start of sequence
  TILING tiling                         /* tiling of sequence */
);

ARRAY_T *get_seq_comp(ALPH_T *alph, XLATE_T *xlate, char *sequence);

#endif

