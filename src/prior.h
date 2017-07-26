 // This are the data structures used in the adaptive prior method 

#include "alphabet.h"
#include "user.h"
// #define Real double  
// #define RealPrec double

 // This structure stores parameters of the prior distributions 
typedef struct{
  ALPH_T *alph; // Alphabet
  int L; // Number of prior distributions 
  double *Mix; // Mixture coefficents for each prior 
  double *B; // strength of each prior 
  // Prior distributions. L Dirchlet's over AlphaChar positions:
  // Distribution[L][AlphaChar+1]
  double **Distr;        
  // !=0 re-estimate all alpha
  // ==0 re-estimate alpha 0 only
  int *FullUpdate;    
  // !=0 update mixture coefficents
  // ==0 do not update coefficents
  int *QUpdate;       
  char **StructID; // Structure Tag 
  char **Comment;
} PriorLib;

 // Allocate space for a prior library with L priors and Alpha characters 
PriorLib *alloc_PriorLib(int L, ALPH_T *alph);
void free_PriorLib (PriorLib *lib);

 // This reads prior information from a file into a PriorLib 
PriorLib *read_PriorLib(char *plib_name, double desired_beta, ALPH_T *custom_alph);

/* This calculates the regularizer given the observed
   frequencies, the prior library, and a weight of the priors */
extern void mixture_regularizer(
  double *freq, // obs freq 
  PriorLib *Lib, // priors 
  double *reg // pseudo-counts 
);

double logpajgy(double *y, PriorLib *Lib, int j, int Calc);

double logpygaj(double *y, double *a, int AlphLength);

