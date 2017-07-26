/***********************************************************************
 *                                                                      *
 *       MEME++                                                         *
 *       Copyright 1994-2006, The Regents of the University of California *
 *       Author: Timothy L. Bailey                                      *
 *                                                                      *
 ***********************************************************************/

/*
   Routines to map a sequence to a theta matrix.
*/

#include "meme.h"

/**********************************************************************/
/*
   init_map

   Set up an |A| x |A| sequence_letter_to_frequency_vector
   or sequence_letter_to_logodds_vector
   matrix where each column is the 
   mapping from a given letter to a frequency/logodds vector.
   An extra row and column is added for the wildcard character 'X', and
   set to the background frequencies (or 0 if, logodds).

   Two types of mapping are possible:
   Uni -   add n prior
   Pam -   Dayhoff PAM for proteins, transversion/transition
   PAM for DNA

   Returns the map.

*/
/**********************************************************************/
THETA init_map(
    // type of mapping:
    //  Uni - add n prior
    //  Pam - pam matrix
    MAP_TYPE type,                
    // degree of crispness, depends on type,
    //  Uni - pseudo count n to add
    //  Pam - pam distance
    double scale,                 
    ALPH_T *alph, // alphabet
    ARRAY_T *back, // background frequencies for wildcard
    bool lo // setup a logodds mapping if true 
    ) {
  int i, j, p, alen_core;
  THETA map; // the map 

  // dayhoff PAM 1 matrix; order of alphabet: ACDEFGHIKLMNPQRSTVWY 
  // dayhoff_ij = Pr(amino acid j --> amino acid i | time=1) 
  double dayhoff[20][20] = {
    { 9867, 3, 10, 17, 2, 21, 2, 6, 2, 4, 6, 9, 22, 8, 2, 35, 32, 18, 0, 2},
    { 1, 9973, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 5, 1, 2, 0, 3},
    { 6, 0, 9859, 53, 0, 6, 4, 1, 3, 0, 0, 42, 1, 6, 0, 5, 3, 1, 0, 0},
    { 10, 0, 56, 9865, 0, 4, 2, 3, 4, 1, 1, 7, 3, 35, 0, 4, 2, 2, 0, 1},
    { 1, 0, 0, 0, 9946, 1, 2, 8, 0, 6, 4, 1, 0, 0, 1, 2, 1, 0, 3, 28},
    { 21, 1, 11, 7, 1, 9935, 1, 0, 2, 1, 1, 12, 3, 3, 1, 21, 3, 5, 0, 0},
    { 1, 1, 3, 1, 2, 0, 9912, 0, 1, 1, 0, 18, 3, 20, 8, 1, 1, 1, 1, 4},
    { 2, 2, 1, 2, 7, 0, 0, 9872, 2, 9, 12, 3, 0, 1, 2, 1, 7, 33, 0, 1},
    { 2, 0, 6, 7, 0, 2, 2, 4, 9926, 1, 20, 25, 3, 12, 37, 8, 11, 1, 0, 1},
    { 3, 0, 0, 1, 13, 1, 4, 22, 2, 9947, 45, 3, 3, 6, 1, 1, 3, 15, 4, 2},
    { 1, 0, 0, 0, 1, 0, 0, 5, 4, 8, 9874, 0, 0, 2, 1, 1, 2, 4, 0, 0},
    { 4, 0, 36, 6, 1, 6, 21, 3, 13, 1, 0, 9822, 2, 4, 1, 20, 9, 1, 1, 4},
    { 13, 1, 1, 3, 1, 2, 5, 1, 2, 2, 1, 2, 9926, 8, 5, 12, 4, 2, 0, 0},
    { 3, 0, 5, 27, 0, 1, 23, 1, 6, 3, 4, 4, 6, 9876, 9, 2, 2, 1, 0, 0},
    { 1, 1, 0, 0, 1, 0, 10, 3, 19, 1, 4, 1, 4, 10, 9913, 6, 1, 1, 8, 0},
    { 28, 11, 7, 6, 3, 16, 2, 2, 7, 1, 4, 34, 17, 4, 11, 9840, 38, 2, 5, 2},
    { 22, 1, 4, 2, 1, 2, 1, 11, 8, 2, 6, 13, 5, 3, 2, 32, 9871, 9, 0, 2},
    { 13, 3, 1, 2, 1, 3, 3, 57, 1, 11, 17, 1, 3, 2, 2, 2, 10, 9901, 0, 2},
    { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 9976, 1},
    { 1, 3, 0, 1, 21, 0, 4, 1, 0, 1, 0, 3, 0, 0, 0, 1, 1, 1, 2, 9945}
  };

  // transition/transversion PAM 1 matrix; alphabet order "ACGT" 
  double trans[4][4] = {
    { 9900, 20, 60, 20}, 
    { 20, 9900, 20, 60},
    { 60, 20, 9900, 20},
    { 20, 60, 20, 9900}
  };

  // special PAM mutation frequencies for DNA
  double pam_dna_freq[] = { 
    0.25 /* A */,
    0.25 /* C */,
    0.25 /* G */,
    0.25 // T 
  };

  // special PAM mutation frequencies for proteins 
  double pam_prot_freq[] = { 
    0.096 /* A */,
    0.025 /* C */,
    0.053 /* D */,
    0.053 /* E */,
    0.045 /* F */,
    0.090 /* G */,
    0.034 /* H */,
    0.035 /* I */,
    0.085 /* K */,
    0.085 /* L */,
    0.012 /* M */,
    0.042 /* N */,
    0.041 /* P */,
    0.032 /* Q */,
    0.034 /* R */,
    0.057 /* S */,
    0.062 /* T */,
    0.078 /* V */,
    0.012 /* W */,
    0.030 /* Y */,
    0.0 // X 
  };

  alen_core = alph_size_core(alph);

  // create the array for the sequence to theta map 
  create_2array(map, double, alen_core+1, alen_core+1);

  switch (type) {
    case Uni: {
      double main_letter; // probability of main letter 
      double other; // probability of each other letter 
      main_letter = (1.0 + scale)/(1.0 + alen_core * scale);
      other = scale/(1.0 + alen_core * scale);
      if (VERBOSE) {printf("main= %g\n\n", main_letter);}
      // create a matrix of columns; each column gives mapping for a letter 
      for (i = 0; i < alen_core; i++) 
        for (j = 0; j < alen_core; j++)  
          map[i][j] = (i == j) ? main_letter : other;
    } break;
    case Pam: {
      double mul[20][20]; 
      if (!alph_is_builtin_dna(alph) && !alph_is_builtin_protein(alph)) {
        fprintf(stderr, "PAM mode only works with builtin alphabets.\n");
        exit(1);
      }
      // convert initial matrix to probabilities 
      for (i = 0; i < alen_core; i++)
        for (j = 0; j < alen_core; j++)
          map[i][j] = (alph_is_builtin_dna(alph) ? trans[i][j] : dayhoff[i][j]) / 10000;
      // take pam matrix to desired power 
      // copy: 
      for (i = 0; i < alen_core; i++)
        for (j = 0; j < alen_core; j++)
          mul[i][j] = map[i][j];
      // multiply: 
      while (--scale) {
        double result[20][20], sum;
        for (i = 0; i < alen_core; i++) {
          for (j = 0; j < alen_core; j++) {
            for (sum = p = 0; p < alen_core; p++) {
              sum += mul[i][p] * map[p][j];
            }
            result[i][j] = sum;
          }
        }
        for (j = 0; j < alen_core; j++) {
          for (i = 0; i < alen_core; i++) {
            // map[i][j] = result[i][j];
            RND(result[i][j], 8, map[i][j]);
          }
        }
      }
    }
  }

  // add last row and column for mapping from the wildcard
  for (i = 0; i < alen_core; i++) {
    map[alen_core][i] = map[i][alen_core] = get_array_item(i, back);
  }

  // convert to logodds matrix if requested 
  if (lo) {
    double *pfreq;
    double ufreq = 0, avg, x_avg;
    if (alph_is_builtin_dna(alph)) {
      pfreq = pam_dna_freq;
    } else if (alph_is_builtin_protein(alph)) {
      pfreq = pam_prot_freq;
    } else { 
      // custom alphabets get uniform "special PAM mutation frequencies"
      // since Tim couldn't remember how these fields were calculated
      // and we don't know if this would be better than the dataset background
      ufreq = 1.0 / alen_core;
      pfreq = NULL;
    }
    for (i = 0; i < alen_core; i++) {
      for (j = 0; j <= i; j++) {
        map[i][j] = map[j][i] = NINT(2 * LOG2(map[i][j] / (pfreq ? pfreq[i] : ufreq)));
      }
    }
    // do the last row and column for "X" matches: average of match to all chars 
    x_avg = 0; // average for X row 
    for (i = 0; i < alen_core; i++) {
      avg = 0; // average for row 
      for (j = 0; j < alen_core; j++) {
        avg += map[i][j];
      }
      x_avg += avg / alen_core;
      map[i][alen_core] = map[alen_core][i] = NINT(avg / alen_core);
    }
    map[alen_core][alen_core] = NINT(x_avg/alen_core);
#ifdef DEBUG
    for (i=0; i<=alen_core; i++) {
      for (j=0; j<=alen_core; j++) {
        printf("%3g ", map[i][j]);
      }
      printf("\n");
    }
#endif
  } // lo 

  return map;
} // init_map 


/**
 * convert_to_lmap
 *
 * Converts a matrix of sequence to theta probability mappings into a
 * matrix containing the int logs of those probabilities. Also sets up
 * a vector for mapping an "X" character to a vector of letter probabilities.
 * Those probabilities are uniform across all letters.
 */
void convert_to_lmap (
    THETA map,
    int lmap[MAXALPH][MAXALPH],
    ALPH_T *alph
    )
{  
  /* 
     Set up the matrix of frequency vectors for each letter in the alphabet.
     Column and row for the "match-anything" character X are set to 1.0/alength 
     so that such matches are neither favored nor disfavored, and where 
     they match is irrelevant:
     */
  int i,j;
  for (i = 0; i < alph_size_wild(alph); i++) {
    for (j = 0; j < alph_size_core(alph); j++) {
      lmap[i][j] = (i < alph_size_core(alph)) ? INT_LOG(map[j][i]) : INT_LOG(1.0 / alph_size_core(alph));
    }
    lmap[i][j] = INT_LOG(1.0 / alph_size_core(alph)); // X 
  }
}


/**********************************************************************/
/*
   init_theta

   Set theta to represent a consensus sequence by copying
   columns of the letter_to_theta_column_array map to theta.

   For columns in sequence containing 'X' or if the sequence
   is null, theta is set to the values in the extra (alength)
   column in the map.  Each column is length alength+1, because
   there is a row for 'X'.
   */
/**********************************************************************/
void init_theta(
    THETA theta, // theta 
    uint8_t *start, // integer encoded starting sequence 
    int w, // width of motif 
    THETA map, // frequency vectors for each letter 
    ALPH_T *alph
    )
{
  int i, j, c;

  // initialize the insite frequencies 
  for (i=0; i < w; i++) { // column in theta 
    c = (start != NULL ? start[i] : alph_wild(alph));
    for (j = 0; j < alph_size_wild(alph); j++) { // row in theta 
      theta(i, j) = map[c][j];
    }
  }
} // init_theta 


/**
 * convert_to_ltheta
 *
 * Convert the entries in the specified matrix motif model from doubles to INT
 * LOG values.
 */
void convert_to_ltheta (
    double matrix_ds[MAXSITE][MAXALPH], ///< The input matrix of doubles
    int matrix_il[MAXSITE][MAXALPH], ///< The output matrix of int logs
    int nrows,
    int ncols
    )
{
  int row_idx;
  for (row_idx = 0; row_idx < nrows; row_idx++) {
    int col_idx;
    for (col_idx = 0; col_idx < ncols; col_idx++) {
      matrix_il[row_idx][col_idx] = INT_LOG(matrix_ds[row_idx][col_idx]);
      fprintf(stderr, "%i ", matrix_il[row_idx][col_idx]);
    }
    fprintf(stderr, "\n");
  }
}

