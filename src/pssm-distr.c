#include "macros.h"
#include "pssm-distr.h"
#include "user.h"
#include <assert.h>

/**********************************************************************/
/*
        Calculate the theoretical distribution function for an
        (integer-valued) PSSM.
        Returns an array of pdf values:
          pvalue[x] = Pr(score == x)
        for 0 <= x <= range*w.

*/
/**********************************************************************/
static double *calc_pssm_pdf(
  int w, // width of PSSM 
  int alen, // length of alphabet 
  int range, // largest value in PSSM 
  // scaled, integer PSSM: pssm[i][j] is score for j_th letter in i_th column
  // of motif; entries in PSSM are in range [0..range]
  double **pssm,
  ARRAY_T *prob // at least a 0-order Markov background model with alen entries
)
{
  // Check validity of prob/psfm combination:
  assert(prob != NULL);
  assert(get_array_length(prob) >= alen);

  int i, j, k;
  int size = (w * range) + 1; // size of pdf array 
  double *pdf_old = NULL, *pdf_new = NULL;

  // set up the two arrays to hold probability density functions 
  Resize(pdf_old, size, double);
  Resize(pdf_new, size, double);

  // set probabilities of each new score to zero except for score 0 
  pdf_new[0] = 1;
  for (i = 1; i < size; i++) pdf_new[i] = 0;

  // recursively compute the pdf 
  for (i = 0; i < w; i++) { // loop over columns in motif 
    int max_score = i * range; // maximum possible cumulative score 
    SWAP(double *, pdf_new, pdf_old); // new column: swap old and new pdfs 
    // zero out the new pdf; new maximum score is old max + range 
    for (k = 0; k <= (max_score + range); k++) pdf_new[k] = 0;
    for (j = 0; j < alen; j++) { // loop over letters 
      int score = (int) pssm[i][j]; // get integer PSSM entry 
      for (k = 0; k <= max_score; k++) {
        if (pdf_old[k] != 0) {
          pdf_new[k + score] += pdf_old[k] * get_array_item(j, prob);
        }
        // sanity check
        if ((k + score) >= size) {
          fprintf(stderr,
            "calc_pssm_pdf error: i=%d j=%d k=%d max_score=%d score=%d size=%d\n",
            i, j, k, max_score, score, size);
          return NULL;
        }
      }
    }
  }

  // clean up 
  myfree(pdf_old);

  // return the pdf 
  return pdf_new;
}

/**********************************************************************/
/*
        calc_pssm_cdf

        Calculate the 1 minus the theoretical cumulative distribution
        function for an (integer-valued) PSSM given a 0-order Markov
        background model.  PSSM entries must be in range [0..range].

        Returns an array of 1-cdf values:
          pvalue[x] = Pr(score >= x)
        for 0 <= x <= range*w.
*/
/**********************************************************************/
double *calc_pssm_cdf(
  int w, // width of PSSM 
  int alen, // length of alphabet 
  int range, // largest value in PSSM 
  // scaled, integer PSSM: pssm[i][j] is score for j_th letter in i_th column 
  // of motif; entries in PSSM are in range [0..range]
  double **pssm,
  ARRAY_T *bg // 0-order Markov background model 
)
{
  int i, size;
  double *pdf; // probability distribution 

  pdf = calc_pssm_pdf(w, alen, range, pssm, bg);
  if (!pdf) return NULL;

  // compute 1-cdf from the pdf from the right to preserve right accuracy 
  size = (w * range) + 1; // size of cdf array 
  for (i = size - 2; i >= 0; i--) {
    pdf[i] += pdf[i + 1];
    pdf[i] = MIN(1.0, pdf[i]);
  }
  // return the cdf 
  return pdf;
} // calc_pssm_cdf 
