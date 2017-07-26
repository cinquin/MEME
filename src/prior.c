#include "meme.h"
#define MAXS    200
static BOOLEAN first_time = TRUE;
// for getting rid of gcc's annoying warn_unused_result
static inline void ignore_value (int i) { (void) i; }
static inline void ignore_ptr (void* p) { (void) p; } 
static inline void test_value (int test) { if (!test) die("IO function failed\n"); }

#define LogAddLog(x, y) ((x) > (y) ? LogAddLog1((x),(y)) : LogAddLog1((y),(x)))
#define LogAddLog1(x,y) ((x)-(y) > BITS ? (x) : (x) + log(1+exp((y)-(x))))

/*
 *  L       Number of distributions
 *  alph    Alphabet
 */
PriorLib *alloc_PriorLib(int l, ALPH_T *alph)
{
  PriorLib *temp;
  int i;

  temp = (PriorLib *)mymalloc(sizeof(PriorLib));
  temp->alph = alph_hold(alph);
  temp->L = l;

  temp->Mix = (double *)mymalloc(sizeof(double)*l);
  temp->B = (double *)mymalloc(sizeof(double)*l);
  temp->FullUpdate = (int *)mymalloc(sizeof(int)*l);
  temp->QUpdate = (int *)mymalloc(sizeof(int)*l);

  temp->StructID = (char **)mymalloc(sizeof(char *)*l);
  temp->Comment = (char **)mymalloc(sizeof(char *)*l);
  temp->Distr = (double **)mymalloc(sizeof(double *)*l);
  for (i = 0; i < l; i++) {
    temp->Distr[i] = (double *)mymalloc(sizeof(double) * (alph_size_wild(alph)));
    temp->StructID[i] = (char *)mymalloc(sizeof(char) * MAXS);
    temp->Comment[i] = (char *)mymalloc(sizeof(char) * MAXS);
  } // endfor 

  return(temp);
}

/*
 *    plib    PriorLib struct to deallocate
 */
void free_PriorLib(PriorLib *plib) 
{
  int i;
  if (plib == NULL) return;
  alph_release(plib->alph);
  free(plib->Mix);
  free(plib->B);
  free(plib->FullUpdate);
  free(plib->QUpdate);
  for (i = 0; i < plib->L; i++) {
    free(plib->Distr[i]);
    free(plib->StructID[i]);
    free(plib->Comment[i]);
  }
  // set pointers to null
  memset(plib->Distr, 0, sizeof(double *) * plib->L);
  memset(plib->StructID, 0, sizeof(char *) * plib->L);
  memset(plib->Comment, 0, sizeof(char *) * plib->L);
  free(plib->Distr);
  free(plib->StructID);
  free(plib->Comment);
  // set pointers to null
  memset(plib, 0, sizeof(PriorLib));
  free(plib);
}

/*
 *  plib_name       name of prior library file
 *  desired_beta    >  0, scale \beta_{i,j} so
 *                      \sum_{i=0}^L \lambda_i \sum_{j=1}^20 \beta_{i,j}
 *                      has this value
 *                  == 0, don't scale prior
 *                  <  0, just get alphabet
 */
PriorLib *read_PriorLib(char *plib_name, double desired_beta, ALPH_T *custom_alph)
{
  int i,j, line=0;
  int l;
  PriorLib *temp;
  char input[MAXS], foo[MAXS], alphabet[MAXALPH+1], checkstr[81], *token;
  double x;
  FILE *fp;
  ALPH_T *alph;

  // tlb 
  fp = fopen(plib_name, "r");
  if (!fp) {
    fprintf(stderr, "Can't find prior library %s\n", plib_name);
    exit(1);
  }

  token = "Alphabet="; line++;
  test_value(fscanf(fp,"%s %s\n", checkstr, alphabet) == 2);
  if (strcmp(checkstr, token)) {
    fprintf(stderr, "Line %d of prior library file \n %s \n"
        "should start with \"%s\" "
        "but it starts with \"%s\".\n", line, plib_name, token, checkstr);
    exit(1);
  }
  // determine alphabet
  if (custom_alph == NULL) {
    alph = alph_type(alphabet, 30);
    if (alph == NULL) {
      fprintf(stderr, "The partial alphabet specified in the prior library file"
          " does not match a built-in alphabet and no complete alphabet was specified.\n");
      exit(1);
    }
  } else {
    int alen_core;
    alen_core = strlen(alphabet);
    i = 0;
    if (alen_core == alph_size_core(custom_alph)) {
      for (i = 0; i < alen_core; i++) {
        if (!alph_is_concrete(custom_alph, alphabet[i])) break;
      }
    }
    if (i == 0 || i < alen_core) {
      fprintf(stderr, "The partial alphabet specified in the prior library file"
          " does not match the complete alphabet specified.\n");
      exit(1);
    }
    alph = alph_hold(custom_alph);
  }

  token = "NumDistr="; line++;
  test_value(fscanf(fp,"%s %d\n", checkstr, &l) == 2);
  if (strcmp(checkstr, token)) {
    fprintf(stderr, "Line %d of prior library file \n %s \n"
        "should start with \"%s\" "
        "but it starts with \"%s\"\n.", line, plib_name, token, checkstr);
    exit(1);
  }

  temp = alloc_PriorLib(l, alph);

  if (desired_beta < 0) {
    fclose(fp);
    return(temp);
  }

  for (i = 0; i < temp->L; i++) {
    // Get rid of number= 
    ignore_value(fscanf(fp,"%*s %*s\n"));
    // Mixture 
    ignore_value(fscanf(fp,"%*s"));
    test_value(fscanf(fp,"%lf\n", &x) == 1);
    temp->Mix[i] = x;
    // B (strength) 
    ignore_value(fscanf(fp,"%*s"));
    test_value(fscanf(fp,"%lf\n", &x) == 1);
    temp->B[i] = x;
    // Alpha 
    temp->Distr[i][0] = temp->B[i];
    ignore_value(fscanf(fp,"%*s"));
    for (j = 1; j < alph_size_wild(alph); j++) {
      test_value(fscanf(fp,"%lg", &x) == 1);
      temp->Distr[i][j] = x * temp->B[i];
    }
    // FullUpdate 
    ignore_value(fscanf(fp,"%*s"));
    test_value(fscanf(fp,"%d\n", &(temp->FullUpdate[i])) == 1);
    // QUpdate 
    ignore_value(fscanf(fp,"%*s"));
    test_value(fscanf(fp,"%d\n", &(temp->QUpdate[i])) == 1);
    // StructID 
    test_value(fgets(input, MAXS, fp) != NULL);
    test_value(sscanf(input,"%s",foo) == 1);
    input[strlen(input)-1] = '\0';
    strcpy( (temp->StructID[i]), (input + strlen(foo)) );
    // Comments 
    test_value(fgets(input, MAXS, fp) != NULL);
    test_value(sscanf(input,"%s",foo) == 1);
    strcpy( (temp->Comment[i]), (input + strlen(foo)) );
  }

  // tlb; scale beta to desired value 
  if (desired_beta > 0) {
    int i, j;
    double beta = 0;
    double scale;
    for (i=0; i<temp->L; i++) {
      beta += temp->Mix[i] * temp->B[i];
    }
    /*printf("beta = %10.6f\n", beta);*/
    scale = desired_beta/beta;
    for (i=0; i<temp->L; i++) {
      for (j=0; j < alph_size_wild(alph); j++) {
        temp->Distr[i][j] *= scale;
      }
    }
  }
  fclose(fp);

  alph_release(alph);
  return(temp);
}

/*
 *  freq  obs freq 
 *  Lib   priors
 *  reg   pseudo-counts
 */
void mixture_regularizer(double *freq, PriorLib *Lib, double *reg) {
  double f[MAXALPH+1], sum, tmp;
  int i, j;
  // double  logpajgy();

  // Put frequencies into array with f[0] = sum f_i 
  sum = 0.0;
  for (i = 0; i < alph_size_core(Lib->alph); i++) {
    sum += freq[i];
    f[i+1] = freq[i];
  }
  f[0] = sum;

  // Calculate probs 
  logpajgy(f, Lib, 0, 1);

  // Calculate new regularizer 
  for (i = 0; i < alph_size_core(Lib->alph); i++) {
    reg[i] = 0.0;
    for (j = 0; j < Lib->L; j++) {
      tmp = (exp(logpajgy(f, Lib, j, 0)))*
        ((Lib->Distr[j])[i+1]); // skip A0 
      reg[i] += tmp; 
    }
  }
}

/*
 * This function computes log(p(a^j|y)) used in the calculation of theta. 
 * It is defined to be
 * \log(\frac{q_j p(y given \alpha^j)}{\sum_k q_k p(y given \alpha^k})
 *
 * y        observed frequencies
 * Library  Library of priors
 * j        j'th prior to examine
 * Calc     if == 1 calculate probs
 */
double logpajgy(double *y, PriorLib *Library, int j, int Calc) {
  int i;
  double tmp;
  static double logprob[MAXS], logdenom; // Holders for probabilities 

  // Calculate log probs if not already done 
  if (Calc) {
    tmp = log(Library->Mix[0]) + logpygaj(y,Library->Distr[0],
        alph_size_core(Library->alph));
    logdenom = tmp;
    logprob[0] = tmp;

    // Do remaining terms 
    for (i = 1; i < Library->L; i++) {
      tmp = (log(Library->Mix[i]) + 
          logpygaj(y, Library->Distr[i], alph_size_core(Library->alph)));

      logdenom = LogAddLog(logdenom, tmp);
      logprob[i] = tmp;
    }
  }
  return(logprob[j] - logdenom);
}



#define MAXX1 100
#define DELTA1 .001
#define MAXX2 100000
#define DELTA2 1.0

// Change if MAXX or DELTA constants change */
#define SIZE1 100000 // MAXX1 / DELTA1
#define SIZE2 100000 // MAXX2 / DELTA2

static double lgam_array1[SIZE1 + 2];
static double lgam_array2[SIZE2 + 2];

static double lgam(double x) {
  if (x >= DELTA1 && x <= MAXX1) { 
    int i = (int) (x/DELTA1);
    return(lgam_array1[i] + (lgam_array1[i+1] - lgam_array1[i])/2);
  } else if (x > MAXX1 && x <= MAXX2) { 
    int i = (int) (x/DELTA2);
    return(lgam_array2[i] + (lgam_array2[i+1] - lgam_array2[i])/2);
  } else {
    return(lgamma(x));
  }
} // lgam 

/*
 * This function computes log(p(y|a^j)) used in the calculation of theta. 
 * It is defined to be
 * \log(\frac{\Gamma(n+1)\Gamma(\alpha_0)}{\Gamma(n+\alpha_0)}
 * \prod_{i=1}^{20}\frac{\Gamma(y_i+\alpha_i)}{\Gamma(y_i+1)\Gamma(\alpha_i)})
 *
 * y          observed frequencies
 * a          distribution parameters
 * AlphLength length of alphabet
*/
double logpygaj(double *y, double *a, int AlphLength) {
  int i;
  double temp;

  // set up array of values of lgamma to save time 
  if (first_time) {
    double x;
    for (i=1, x=0; i<=MAXX1/DELTA1 + 1; i++) { 
      x += DELTA1;
      lgam_array1[i] = lgamma(x);
      lgam_array1[i] = lgam_array1[i];
    }
    for (i=1, x=0; i<=MAXX2/DELTA2 + 1; i++) { 
      x += DELTA2;
      lgam_array2[i] = lgamma(x);
      lgam_array2[i] = lgam_array2[i];
    }
    first_time = FALSE;
  }

  temp=0.0;

  temp+= lgamma(y[0]+1.0);
  temp+= lgamma(a[0]);
  temp+= -lgamma(y[0]+a[0]);

  for (i=1; i<=AlphLength; i++) {
    temp+= lgamma(y[i]+a[i]);
    temp+= -lgamma(y[i]+1.0);
    temp+= -lgamma(a[i]);
  }

  return(temp);
}

