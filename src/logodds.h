/*
 * $Id: logodds.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:41:11  nadya
 * Initial revision
 *
 */

#ifndef LOGODDS_H
#  define LOGODDS_H

#include "user.h"
#include "motif.h"
#include "array-list.h"

/*
  LOGODDS matrices
*/
 // number of log-odds matrices allowed 
#define MAXLO MAXG
/* 
  single-letter logodds matrix 
*/
#define logodds(a, b)   logodds[(a)][(b)]
#define logodds1(a, b)   logodds1[(a)][(b)]
#define logodds2(a, b)   logodds2[(a)][(b)]
//typedef double LOGODDSB; // base type of LOGODDS 
//typedef LOGODDSB **LOGODDS; // [pos][chash(letter)] 

/* 
  two-letter logodds matrix 
*/
//typedef double LOGODDS2B; // base type of LOGODDS2 
//typedef LOGODDS2B **LOGODDS2; // [pos][dhash(chash(l1),chash(l2),alen)] 

 // macro to convert from scaled score back to bit-score 
#define scaled_to_bit(x,w,scale,offset) (((x)/(scale)) + ((w)*(offset)))
#define bit_to_scaled(x,w,scale,offset) (NINT(((x) - ((w)*offset)) * (scale)))

/*
  macro to compute the three class bit log-odds score from the positive and
  negative motif scaled/offset scores:
        pos             positive motif scaled log-odds score:log_2(Pr(P)/Pr(B))'
        neg             negative motif scaled log-odds score: log_2(Pr(P)/Pr(N)'
        lo              pointer to logodds structure for motif
*/
#define score3class(pos, neg, lo) \
  ( -MEME_LOG_SUM( \
    -((scaled_to_bit(pos, lo->w, lo->scale, lo->offset) \
       * (Log2)) - (lo->ln_lambda1)), \
    -((scaled_to_bit(neg, lo->w, lo->scalen, lo->offsetn) \
       * (Log2)) - (lo->ln_lambda2))  \
  ) / (Log2) )

 // structure used by read_log_odds 
typedef struct lo {
  ALPH_T *alph;
  // this field exists because I can see a legitimate reason to want a logodds
  // table with and without entries for ambiguous symbols. I considered having
  // 2 booleans or an enumeration to note which length was in use but I decided
  // that was silly and people were responsible enough not to abuse this field
  // by putting in nonsensical values.
  int alen; // could equal alph_size_core, alph_size_wild or alph_size_full
  bool pair; // pair of motifs: m_ij/b_j followed by m_ij/n_ij 
  char *meme_name; // the ID from the meme file 
  char *meme_id2; // the ID2 from the meme file 
  int imotif; // loading order of the motif 
  int w; // width of motif 
  int ws; // width of motif in sequence 
  double thresh; // threshold for scores 
  double ev; // E-value of motif 
  double e_ll; // expected log likelihood of motif 
  double ic; // information content of motif 
  double sites; // number of occurrences of motif in dataset 
  bool pal; // motif is palindrome if true (requires complementable alphabet) 
  bool invcomp; // use reverse complement strand as well 
  double lambda; // lambda for motif 
  double L; // average length of sequences 
  char *best_seq; // best possible matching sequence 
  char *best_icseq; // inverse complement of best possible match sequence 
  double **logodds; // log-odds matrix 
  double **logodds2; // two-letter log-odds matrix 
  bool is_bad; // correlation with other motifs is bad 
  double scale; // scale factor (positive motif) for converting 2 bits
  double offset; // offset (positive motif) for converting 2 bits 
  double scalen; // scale factor (negative motif) for converting 2 bits
  double offsetn; // offset (negative motif) for converting 2 bits
  double scale3; // scale factor for converting 3-class score 2 bits 
  double offset3; // offset factor for converting 3-class score 2 bits 
  double ln_lambda1; // log( Pr(B)/Pr(P) ) 
  double ln_lambda2; // log( Pr(N)/Pr(P) ) 
} LO;

LO** convert_2_log_odds(
  ALPH_T *alph, // the motif alphabet
  XLATE_T *xlate, // translate the sequence to alph
  ARRAY_T *bgfreqs, // 0-order bg frequencies for alphabet (pointer) 
  ARRAYLST_T *motifs_in, // MOTIF_T motifs
  bool wildcard_only, // should the resulting logodds have the full alphabet
  int range, // scale entries in logodds matrices to [0..range] 
  int *nmotifs // number of motifs returned (ignored if NULL)
);

void free_log_odds(
  LO *lo // log-odds matrix
);

void destroy_log_odds(LO *lo);

void min_max(
  double **logodds, // log-odds matrix 
  int w, // width of motif 
  int a, // length of alphabet 
  double *minimum, // minimum score 
  double *maximum // minimum score 
);

MATRIX_T* motif_corr(
  int nmotifs, // number of motifs 
  LO *los[] // array of logodds structures 
);

void shuffle_cols(
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs // number of log-odds matrices in los 
);

void scale_lo(
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs, // number of log-odds matrices in los 
  int range // set entries in matrices to [0..range] 
);

void make_double_lo(
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs // number of log-odds matrices in los 
); 

#endif

