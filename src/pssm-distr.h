#ifndef pssm_distr_h
#define pssm_distr_h

#include "array.h"

double *calc_pssm_cdf(
  int w, // width of PSSM 
  int alen, // length of alphabet 
  int range, // largest value in PSSM 
  // scaled, integer PSSM: pssm[i][j] is score for j_th letter in i_th column
  // of motif entries in PSSM are in range [0..R
  double **pssm,
  ARRAY_T *prob // 0-order Markov background mode 
);

#endif
