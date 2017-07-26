/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994-2015, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

#include "macros.h"
#include "logodds.h"

static inline LO* convert_motif(ARRAY_T *bgfreqs, int alen, int nsyms, MOTIF_T *mo) {
  int i, j, k;
  LO *lo;

  assert(nsyms >= 1);

  lo = mm_malloc(sizeof(LO));
  memset(lo, 0, sizeof(LO));

  lo->alph = alph_hold(get_motif_alph(mo));
  lo->alen = alen;
  lo->imotif = get_motif_idx(mo);
  lo->meme_name = strdup(get_motif_id(mo));
  lo->meme_id2 = strdup(get_motif_id2(mo));
  lo->w = get_motif_length(mo);
  lo->ws = lo->w * nsyms;   // motif width in sequence 
  lo->ev = get_motif_evalue(mo);
  lo->is_bad = false;
  lo->pair = false; //according to Tim this isn't used anymore, I haven't investigated...

  //create the score matrix and best matching sequence
  lo->best_seq = mm_malloc(sizeof(char) * (lo->w + 1));
  lo->best_icseq = mm_malloc(sizeof(char) * (lo->w + 1));
  lo->logodds = mm_malloc(sizeof(double*) * lo->w);
  for (i = 0; i < lo->w; i++) {
    double best_score = LITTLE; // best score in row 
    int best_index = -1; // best index in row 
    lo->logodds[i] = mm_malloc(sizeof(double) * alen);
    // copy over the scores of the core alphabet
    for (j = 0; j < alph_size_core(lo->alph); j++) {
      lo->logodds[i][j] = get_motif_score(mo, i, j); 
      // I could allow for ambiguous characters in the best match
      // by adding this check to their loop as well
      if (lo->logodds[i][j] > best_score) {
        best_score = lo->logodds[i][j];
        best_index = j;
      }
    }
    // calculate weighted averages for the ambiguous letters
    // note that there should be at minimum a wildcard
    for (; j < alen; j++) {
      int index;
      double sum_weighted_score, sum_weights;
      sum_weighted_score = 0;
      sum_weights = 0;
      for (k = 0; k < alph_ncomprise(lo->alph, j); k++) {
        index = alph_comprise(lo->alph, j, k);
        sum_weighted_score += lo->logodds[i][index] * get_array_item(index, bgfreqs);
        sum_weights += get_array_item(index, bgfreqs);
      }
      lo->logodds[i][j] = (sum_weights > 0 ? sum_weighted_score / sum_weights : 0);

    }
    // build best matching seq 
    lo->best_seq[i] = alph_char(lo->alph, best_index);
  } // position in motif 

  lo->best_seq[i] = '\0'; // terminate best matching sequence 

  // create reverse complement of best sequence if complementable motif 
  if (alph_has_complement(lo->alph)) {                              // DNA motif: reverse comp. 
    strcpy(lo->best_icseq, lo->best_seq);
    invcomp_seq(lo->alph, lo->best_icseq, lo->w);
  } else { // if not complementable then just reverse 
    for (i=0; i<lo->w; i++) lo->best_icseq[i] = lo->best_seq[lo->w-1-i];
    lo->best_icseq[i] = '\0';
  }
  return lo;
  /*
  // these variables are not set
  double thresh; // threshold for scores 
  double e_ll; // expected log likelihood of motif 
  double ic; // information content of motif 
  double sites; // number of occurrences of motif in dataset 
  bool pal; // motif is palindrome if true (requires complementable alphabet) 
  bool invcomp; // use reverse complement strand as well 
  double lambda; // lambda for motif 
  double L; // average length of sequences 
  LOGODDS2 logodds2; // two-letter log-odds matrix 
  double *corr; // correlations with lower-numbered motifs 
  double scale; // scale factor (positive motif) for converting 2 bits
  double offset; // offset (positive motif) for converting 2 bits 
  double scalen; // scale factor (negative motif) for converting 2 bits
  double offsetn; // offset (negative motif) for converting 2 bits
  double scale3; // scale factor for converting 3-class score 2 bits 
  double offset3; // offset factor for converting 3-class score 2 bits 
  double ln_lambda1; // log( Pr(B)/Pr(P) ) 
  double ln_lambda2; // log( Pr(N)/Pr(P) ) 
  */
}

/**********************************************************************/
/*
        Convert the next log-odds matrices read using meme-io.
        Returns the number of matrices converted.
        Updates los and ff.  Sets up alphabet hash tables.

        Ambiguous letters are given logodds scores equal to the
        weighted average of the scores of the possible true letters using
        the background frequencies (ff) as weights.

        If range > 0,
          1) scale and round the log-odds matrices so that all entries are
             in range [0...range].
          2) Matrix entries will be rounded to log_10(range) significant digits.
          3) Bit score values  can be restored (with some loss of significance)
               score_bit = (score/scale) + (w*offset)

        See (install-path>/bin/make_logodds for format of logodds file.
*/
/**********************************************************************/
LO** convert_2_log_odds(
  ALPH_T *alph, // the motif alphabet
  XLATE_T *xlate, // the translator of the sequence alphabet to the motif alphabet
  ARRAY_T *bgfreqs, // 0-order bg frequencies for alphabet (pointer) 
  ARRAYLST_T *motifs_in, // MOTIF_T motifs
  bool wildcard_only, // should the resulting logodds have the full alphabet or just the core + wildcard
  int range, // scale entries in logodds matrices to [0..range] 
  int *nmotifs // number of motifs returned (ignored if NULL)
) {
  int n, alen, nsyms;
  LO **los; // list of log odds structures
  LO *lo; // log odds structure pointer 
  MOTIF_T *mo;
  
  // allocate space for log odds output
  los = mm_malloc(sizeof(LO*) * arraylst_size(motifs_in));

  // calculate alphabet length
  alen = (wildcard_only ? alph_size_wild(alph) : alph_size_full(alph));

  // calculate number of symbols per translated symbol
  nsyms = (xlate != NULL ? xlate_src_nsyms(xlate) : 1);

  // iterate over all MEME motifs
  for (n = 0; n < arraylst_size(motifs_in); n++) { // number of logodds matrix 
    mo = (MOTIF_T*)arraylst_get(n, motifs_in);
    los[n] = convert_motif(bgfreqs, alen, nsyms, mo);
  }

  // Convert motif matrices to integer so entries are in [0...range].
  if (range > 0) scale_lo(los, n, range);

  //Create a double-letter log-odds matrix for efficiency.
  make_double_lo(los, n);

  if (nmotifs != NULL) *nmotifs = n;
  return los;
}

/**********************************************************************/
/*
   Frees the memory allocated for a log odds object.
*/
/**********************************************************************/
void free_log_odds(LO *lo) {
  int nrows;
  // free the score matrices
  free_2array(lo->logodds, lo->w);
  if (lo->logodds2 != NULL) {
    nrows = (lo->w + 1) / 2;
    if (lo->pair) nrows *= 2;
    free_2array(lo->logodds2, nrows);
  }
  // free the best sequence snippits
  free(lo->best_seq);
  free(lo->best_icseq);
  // free the motif name
  free(lo->meme_name);
  free(lo->meme_id2);
  // dereference the alphabet
  alph_release(lo->alph);
}

/**********************************************************************/
/*
   Frees the memory allocated for a log odds object.
*/
/**********************************************************************/
void destroy_log_odds(LO *lo) {
  free_log_odds(lo);
  memset(lo, 0, sizeof(LO));
  free(lo);
}

/**********************************************************************/
/*
        read_log_odds

        Read in the next log-odds matrices from a file.
        Returns the number of matrices read.
        Updates los and ff.  Sets up alphabet hash tables.

  ALTERATION TO INTERFACE, on 11-08-06:
  If the new psfm is non-null, then the logodds are converted to
  letter frequencies by taking the exponent and multiplying by
  the background frequencies. If this is done, then the psfm
  MUST be of an appropriate size to contain all the frequency data.

        Ambiguous letters are given logodds scores equal to the
        weighted average of the scores of the possible true letters using
        the background frequencies (ff) as weights.

        If range > 0,
          1) scale and round the log-odds matrices so that all entries are
             in range [0...range].
          2) Matrix entries will be rounded to log_10(range) significant digits.
          3) Bit score values  can be restored (with some loss of significance)
               score_bit = (score/scale) + (w*offset)

        See (install-path>/bin/make_logodds for format of logodds file.
*/
/**********************************************************************/
/*
int read_log_odds(
  bool translate_dna, // DNA sequences and protein motifs 
  char *filename, // file name (output of make_logodds) 
  char *alphabet, // alphabet of log-odds matrices 
  char *blast_alphabet, // corresponding BLAST alphabet 
  int *p[MAXASCII], // alphabet permutation/substitution matrix 
  int range, // scale entries in logodds matrices to [0..range] 
  LO *los[MAXLO+1], // log-odds structures 
  double *ff, // null letter frequencies for alphabet (pointer) 
  double psfms[MAXLO][MAXSITE][MAXALPH]
 // If non-null, store the psfms for the pssms here 
)
{
  int i, j, k, len, tmp, pair;
  int nmotifs; // number of motifs 
  static FILE *fptr; // file pointer 
  LO *lo; // log odds structure pointer 
  int alen = strlen(blast_alphabet); // length of BLAST alphabet 

  // open the log-odds file
  fptr = fopen(filename, "r");
  if (fptr == NULL) {
    fprintf(stderr, "Cannot open file `%s'.\n", filename);
    exit(1);
  }

  // read in all the log-odds matrices
  for (nmotifs=0; ; nmotifs++) { // number of logodds matrix 
    lo = NULL; Resize(lo, 1, LO); // create a log odds structure 
    los[nmotifs] = lo; // put pointer to it in the array 

    // Too many logodds matrices?
    if (nmotifs > MAXLO) {
      fprintf(stderr,
        "Too many logodds matrices.  Recompile with larger MAXLO.\n");
      exit(1);
    }

    // read in header line
    if ( (i=fscanf(fptr, "%d %d %le %d",
      &(lo->w), &(lo->alen), &(lo->ev), &pair)) != 4 ) {
      if (nmotifs == 0) {
        fprintf(stderr, "Error reading log-odds matrix file `%s'.\n", filename);
        exit(1);
      } else {
        break;
      }
    } else {
      // there's no scanf specifier for bool so we have to read it in as int
      // and then assign it... it would probably have worked anyway
      // but I don't like warnings.
      lo->pair = (pair != 0);
    }

    // create the score matrix and best matching sequence
    lo->logodds = NULL; Resize(lo->logodds, lo->w, double *);
    lo->best_seq = NULL; Resize(lo->best_seq, lo->w+1, char);
    lo->best_icseq = NULL;
    for (i=0; i < lo->w; i++) {
      double tmp;
      double best_score = LITTLE; // best score in row 
      char best_letter = 'X'; // best letter in row 
      double *row = NULL; Resize(row, alen, double); // temporary row 
      lo->logodds[i] = NULL; Resize(lo->logodds[i], lo->alen, double); // row 
      for (j=0; j < lo->alen; j++) { // pos in old alph 
        if ( fscanf(fptr, "%lf", &tmp) != 1) { // read score 
          if (i!=0 || j!=0) {
            fprintf(stderr, "Error reading log-odds matrix file `%s'.\n",
              filename);
            exit(1);
          }
        }
        lo->logodds[i][j] = tmp;
      }

      // permute the columns of the row so they are in alphabetical order
      // and calculate the scores for any missing, ambiguous letters;
      // determine best matching letter for row
      for (j=0; j<alen; j++) { // position in new alphabet 
        double score = 0;
        double freq = 0;
        int old_p, new_p;
        // no need for weighted average if only one possible letter 
        if (p[j][1] < 0) {
          row[j] = lo->logodds[i][p[j][0]]; // substitute single score 
        } else {
          // calculate weighted average of matching letter scores 
          for (k=0; (old_p = p[j][k]) >= 0; k++) { // p list 
            new_p = hash(alphabet[old_p]); // new position of letter 
            score += lo->logodds[i][old_p] * ff[new_p]; // weighted sum 
            freq += ff[new_p]; // total frequency 
          } // p list 
          row[j] = score / freq; // weighted average 
        }
        if (p[j][1] < 0 && row[j] > best_score) { // update best sgl.let score 
          best_score = row[j]; // best score 
          best_letter = blast_alphabet[j]; // best letter 
        }
      } // position in new alphabet 
      myfree(lo->logodds[i]); // free old space 
      lo->logodds[i] = row; // replace row with permuted/substituted row 
      lo->best_seq[i] = best_letter; // build best matching seq 
    } // position in motif 

    // initialize lo flags 
    lo->best_seq[i] = '\0'; // terminate best matching sequence 
    lo->alen = alen; // set length of alphabet to new alphabet's 
    lo->dna = !strcmp(blast_alphabet, DNAB); // true if DNA motif 
    lo->is_bad = false;
    // measure the length of the string from nmotifs 
    len = 1, tmp = nmotifs+1;
    while (tmp >= 1) { 
      tmp /= 10;
      ++len;
    }
    // set the name of the motif to be the order in which it was loaded 
    lo->meme_name = mymalloc(len*sizeof(char));
    snprintf(lo->meme_name, len, "%d", nmotifs+1);
    lo->imotif = nmotifs+1;     // loading order of motif 

    if (lo->pair) lo->w /= 2; // double matrix-- divide width in half 
    lo->ws = lo->w * (translate_dna ? 3 : 1); // motif width in sequence 

    // create reverse complement of best sequence if DNA motif 
    Resize(lo->best_icseq, lo->w+1, char);
    if (lo->dna) { // DNA motif: reverse comp. 
      strcpy(lo->best_icseq, lo->best_seq);
      invcomp_dna(lo->best_icseq, lo->w);
    } else { // protein motif: reverse 
      for (i=0; i<lo->w; i++) lo->best_icseq[i] = lo->best_seq[lo->w-1-i];
      lo->best_icseq[i] = '\0';
    }
  } // motif 

  // 11-08-06: If psfms are required, then convert every pssm to its
  // corresponding psfm by using the background frequencies:
  if (psfms != NULL) {
    // Convert each log-odds matrix (ie pssm) to a psfm:
    int pssm_idx;
    for (pssm_idx = 0; pssm_idx < nmotifs; pssm_idx++) {
      convert_to_psfm(los[pssm_idx], ff, psfms[pssm_idx]);
    }
  }

  // Convert motif matrices to integer so entries are in [0...range].
  if (range>0) scale_lo(los, nmotifs, range);

  // Create a double-letter log-odds matrix for efficiency.
  make_double_lo(los, nmotifs);

  return nmotifs;
} // read_log_odds 
*/

/**********************************************************************/
/*
        min_max

        Find the lowest and highest scores possible with a
        given log-odds matrix.  Returns the single lowest entry in
        the log-odds matrix.
*/
/**********************************************************************/
void min_max(
  double **logodds, // log-odds matrix 
  int w, // width of motif 
  int a, // length of alphabet 
  double *minimum, // minimum score 
  double *maximum // minimum score 
)
{
  int i, j;
  double min_score = 0;
  double max_score = 0;

  for (i=0; i < w; i++) {
    double min = BIG;
    double max = LITTLE;
    for (j=0; j < a; j++) {
      min = MIN(min, logodds(i, j));
      max = MAX(max, logodds(i, j));
    }
    min_score += min;
    max_score += max;
  }
  *minimum = min_score;
  *maximum = max_score;
} // min_max 

/***********************************************************************/
/*
        motif_corr

        Compute correlations between pairs of motifs.
        The correlation between two motifs is the maximum sum of
        Pearson's correlation coefficients for aligned columns divided
        by the width of the shorter motif.  The maximum is found by
        trying all alignments of the two motifs.

        The correlations are saved in a matrix.
*/
/***********************************************************************/
MATRIX_T* motif_corr(
  int nmotifs, // number of motifs 
  LO *los[] // array of logodds structures 
)
{
  int i, j, k, l, m, n, o, alen_core;
  double **means; // means of motif columns 
  char *alphabet;
  ALPH_T *alph;
  MATRIX_T *corr;
  // initialise the correlation matrix
  corr = allocate_matrix(nmotifs, nmotifs);
  // only use non-ambiguous chars in correlation 
  alph = los[0]->alph;
  alen_core = alph_size_core(alph);

  // compute the motif column means 
  means = mm_malloc(sizeof(double*) * nmotifs);
  for (i = 0; i < nmotifs; i++) {
    int w = los[i]->w; // width of motif i 
    means[i] = mm_malloc(sizeof(double) * w);
    for (j = 0; j < w; j++) { // motif column 
      means[i][j] = 0;
      for (k = 0; k < alen_core; k++) {
        means[i][j] += los[i]->logodds[j][k];
      }
      means[i][j] /= alen_core;
    }
  }

  // compute the maximum sum of Pearson's correlation coefficient for motif pairs
  for (i = 0; i < nmotifs; i++) { // "from" motif 
    for (j = 0; j < i; j++) { // "to" motif 
      double rsum_max = LITTLE; // max. sum of r 
      LO *lo1, *lo2;
      int w1, w2;
      double *mu1, *mu2;
      double corr_val;
      for (o = 0; o < 2; o++) { // align i to j, then j to i 
        int m1, m2;
        if (o == 0) { // align motif i to motif j 
          m1 = i;
          m2 = j;
        } else { // align motif j to motif i 
          m1 = j;
          m2 = i;
        }
        lo1 = los[m1];
        lo2 = los[m2];
        w1 = lo1->w;
        w2 = lo2->w;
        mu1 = means[m1];
        mu2 = means[m2];
        for (k = 0; k < w2; k++) { // alignment 
          double rsum = 0; // sum of corr. coeffs 
          for (l = 0, m = k; l < w1 && m < w2; l++, m++) { // column pair 
            // correlation of column l in motif 1 and column m in motif 2 
            double sum1 = 0, sum2 = 0, sum3 = 0;
            double r, denom;
            for (n = 0; n < alen_core; n++) { // letter 
              double a = lo1->logodds[l][n] - mu1[l];
              double b = lo2->logodds[m][n] - mu2[m];
              sum1 += a * b;
              sum2 += a * a;
              sum3 += b * b;
            } // letter 
            denom = sqrt(sum2 * sum3);
            r = (denom ? (sum1 / denom) : 1);
            rsum += r;
          } // column pair 
          rsum_max = MAX(rsum_max, rsum);
        } // alignment 
      } // i to j, j to i 
      corr_val = rsum_max / MIN(los[i]->w, los[j]->w);
      set_matrix_cell(i, j, corr_val, corr);
      set_matrix_cell(j, i, corr_val, corr);
    } // to motif 
  } // from motif 

  // throw away means
  for (i = 0; i < nmotifs; i++) {
    free(means[i]);
    means[i] = NULL;
  }
  free(means);
  return corr;
} // motif_corr 

/**********************************************************************/
/*
        make_double_lo

        Create a double-letter logodds matrix for efficiency.
*/
/**********************************************************************/
void make_double_lo(
  LO **los, // array of pointers to log-odds matrices 
  int nmotifs // number of log-odds matrices in los 
)
{
  int i, j, k, imotif, isym;
  int w, alen, ncols, nrows;
  double **logodds, **logodds2;
  LO *lo;

  for (imotif = 0; imotif < nmotifs; imotif++) { // each motif 
    lo = los[imotif]; // motif 
    w = lo->w; // width of motif 
    alen = lo->alen; // length of motif alphabet (this is the entire alphabet that exists after hashing, so including wildcards)
    ncols = (alen+1)*(alen+1)+1; // columns (letters) in double matrix 
    nrows = (w+1)/2; // rows in hashed matrix 
    logodds = lo->logodds; // single-letter matrix 

    // Create the double-letter logodds matrix that gives scores for
    // each possible pair of letters.  Alphabet length is extended
    // by one for the "blank" character necessary when the sequence length
    // or motif length is odd.
    if (lo->pair) nrows *= 2; // double number of rows if 2 motifs 
    create_2array(logodds2, double, nrows, ncols);
    for (i = 0; i < w; i += 2) { // motif position 
      for (j = 0; j < alen; j++) { // letter in position+0 
        for (k = 0; k <= alen; k++) { // letter in position+1 
          // Note when k = alen it is the special "blank" character NOT THE WILDCARD!
          isym = (j * (alen + 1)) + k;
          logodds2[i/2][isym] = ((i == (w-1) || k == alen) ? logodds[i][j] : logodds[i][j] + logodds[i+1][k]);
          if (lo->pair) {
            logodds2[(nrows/2)+(i/2)][isym] = ((i == (w - 1) || k == alen) ? logodds[w+i][j] : logodds[w+i][j] + logodds[w+i+1][k]);
          }
        } // letter 1 
      } // letter 0 
    } // motif position 
    lo->logodds2 = logodds2;
  } // motif 

} // make_double_lo 

/**********************************************************************/
/*

        scale_lo

        Scale and round the log-odds matrices so that all entries are
        in range [0...range].
        Matrix entries will be rounded to log_10(range) significant digits.

        Bit score values  can be restored (with some loss of significance) by:
                score_bit = (score/scale) + (w*offset)
*/
/**********************************************************************/
void scale_lo(
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs, // number of log-odds matrices in los 
  int range // set entries in matrices to [0..range] 
)
{
  int i, j, imotif;

  /*
    Compute the scale and offset factors for each motif
    and apply them to the logodds matrices to scale them
    to [0..range]
  */
  for (imotif = 0; imotif < nmotifs; imotif++) {
    LO *lo = los[imotif]; // logodds structure 
    int a = lo->alen; // length of alphabet, either alph_size_core, alph_size_wild or alph_size_full
    int w = lo->w; // width of motif 
    int size = w*range+1; // largest possible score 
    double small = BIG; // smallest entry pos logodds matrix 
    double large = -BIG; // largest entry pos logodds matrix 
    double smalln = BIG; // smallest entry neg logodds matrix 
    double largen = -BIG; // largest entry neg logodds matrix 

    // find the smallest/largest entry in the logodds matrix 
    for (i=0; i<w; i++) {
      for (j=0; j<a; j++) {
        small = MIN(small, lo->logodds(i,j));
        large = MAX(large, lo->logodds(i,j));
        if (lo->pair) { // negative motif 
          smalln = MIN(smalln, lo->logodds(w+i,j));
          largen = MAX(largen, lo->logodds(w+i,j));
        }
      }
    }

    // skip this motif if it has no information 
    if (large==small || (lo->pair && largen==smalln)) {
      lo->scale = 0;
      continue;
    }

    // compute scale and offset so that matrix entries in [0..range] 
    lo->scale = range/(large-small); // positive motif 
    RND(lo->scale, RNDDIG, lo->scale); // round to RNDDIG places 
    lo->offset = small;
    if (lo->pair) { // negative motif 
      // set scale factor and offset for 3-class scores 
      lo->scalen = range/(largen-smalln);
      lo->offsetn = smalln;
      lo->ln_lambda1 = log(250000.0); // needed by score3class 
      lo->ln_lambda2 = 0; // needed by score3class 
      small = score3class(0, 0, lo);
      large = score3class(size-1, size-1, lo);
      lo->scale3 =  (size-1)/(large-small);
      RND(lo->scale3, RNDDIG, lo->scale3); // round to RNDDIG places 
      lo->offset3 = small;
    }

    // scale, offset and round logodds matrix entries to range [0..range] 
    for (i=0; i<w; i++) {
      for (j=0; j<a; j++) {
        lo->logodds(i,j) =
          bit_to_scaled(lo->logodds(i,j), 1, lo->scale, lo->offset);
        if (lo->pair) { // negative motif 
          lo->logodds(w+i,j) =
            bit_to_scaled(lo->logodds(w+i,j), 1, lo->scalen, lo->offsetn);
        }
      }
    }

  } // imotif 

} // scale_lo 

/**********************************************************************/
/*
        shuffle_cols

        Shuffle the columns of each motif.
*/
/**********************************************************************/
void shuffle_cols(
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs // number of log-odds matrices in los 
)
{
  int i, j, imotif, maxlen, w, a;
  double **tmp;
  int *permute;
  LO *lo;

  srand_mt(0); // FIXME I don't think it makes sense to seed the RNG inside this function.
  tmp = NULL;
  permute = NULL;
  maxlen = 0;
  a = (nmotifs > 0 ? los[0]->alen : 0); // length of alphabet 
  for (imotif = 0; imotif < nmotifs; imotif++) {
    lo = los[imotif];
    w = lo->w; // width of motif

    // check we have enough space in our arrays
    if (w > maxlen) {
      permute = mm_realloc(permute, sizeof(int) * maxlen);
      tmp = mm_realloc(tmp, sizeof(double*) * maxlen);
      for (i = maxlen; i < w; i++) {
        tmp[i] = mm_malloc(sizeof(double) * a);
      }
      maxlen = w;
    }

    // setup a permute list using an inside-out Fisher-Yates shuffle
    // http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_.22inside-out.22_algorithm
    for (i = 0;  i < w; i++) {
      j = (int)(drand_mt() * (i + 1));
      if (j != i) {
        permute[i] = permute[j];
      }
      // in this case we are using a source which is a ordered list of numbers,
      // hence it is equal to i
      permute[j] = i; //permute[j] = source[i];
    }

    // announce 
    printf("Permuting columns of motif %d: ", imotif+1);
    for (i = 0; i < w; i++) printf("%d ", permute[i]);
    printf("\n");

    // move logodds to tmp in permuted column order 
    for (i = 0; i < w; i++) {
      for (j = 0; j < a; j++) {
        tmp[i][j] = lo->logodds[permute[i]][j];
      }
    }

    // update logodds 
    for (i = 0; i < w; i++) {
      for (j = 0; j < a; j++) {
        lo->logodds[i][j] = tmp[i][j];
      }
    }
  }
  // free storage
  free(permute);
  for (i = 0; i < maxlen; i++) {
    free(tmp[i]);
  }
  free(tmp);
} // shuffle_cols 


/**
 * convert_to_psfm
 *
 * Convert a log odds matrix to a psfm by using the specified background
 * frequencies. The psfm parameter for storing the output must point to
 * a matrix that has had memory allocated prior to the call to this function.
 */
/*
void convert_to_psfm(
  LO *lo,            ///< Matrix of residue log-odds
  double *back,      ///< Background frequencies for determining freqs from los
  double psfm[MAXSITE][MAXALPH]  ///< Matrix in which to store the position
                     ///< specific freqs
)
{
  // Exponentiate each log-odds and then multiply by that letter's background
  // frequency in order to get the position-specific frequency...

  // LOGODDS data type possesses this structure: double[pos][chash(letter)]:
  LOGODDS lo_mat = lo->logodds;

  int col_idx;
  for (col_idx = 0; col_idx < lo->w; col_idx++) {
    int lett_idx;
    double col_sum = 0;
    double *curr_col = psfm[col_idx];
    for (lett_idx = 0; lett_idx < lo->alen; lett_idx++) {
      // The pssm entries are the logodds (base 2) * 100. Compensate for this:
      double curr_logodds_base2 = (lo_mat[col_idx][lett_idx])/100;
      double curr_logodds = curr_logodds_base2 * Log2;
      double curr_backfreq = back[lett_idx];

      // Calculate position-specific frequency from log-odds:
      double curr_psf = exp(curr_logodds)*curr_backfreq;
      col_sum += curr_psf;

      // Store value in appropriate element of output matrix:
      curr_col[lett_idx] = curr_psf;
    }

    // Check that column sums to approximately 1. If not, then report error.
    // Otherwise, normalise:
    if ((col_sum < 1 - ERR_EPSILON) || (col_sum > 1 + ERR_EPSILON)) {
      fprintf(stderr, "INPUT ERROR: probability column does not sum to 1. "
              "Pr = %f\nThis may be due to background frequencies that don't"
              " match this PSSM.\n",
              col_sum);
      exit(1);
    }
    else {
      for (lett_idx = 0; lett_idx < lo->alen; lett_idx++) {
        // Normalise current probability:
        curr_col[lett_idx] = curr_col[lett_idx]/col_sum;
      }
    }
  }
}
*/
