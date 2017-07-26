#ifndef MEME_H
#define MEME_H

#include "config.h"
#include "heap.h"
#include "user.h"
#include "macros.h"
#include "mtype.h"
#include "logs.h"
#include "prior.h"
#include "logodds.h"

// conventions:
//  ALL_CAPS  user defined type
//  Vname   enum type value
//  name()    macro

 // globals 
DEXTERN(int, PAGEWIDTH, 80); // page width for printing must be > MSN + 40 (see user.h) 
DEXTERN(bool, VERBOSE, false); // turn on verbose output mode 
DEXTERN(bool, TRACE, false); // print each start tried 
DEXTERN(bool, PRINT_FASTA, false); // print sites in BLOCKS format 
DEXTERN(bool, PRINTALL, false); // for debug 
DEXTERN(bool, PRINT_W, false); // print w_ij 
DEXTERN(bool, PRINT_Z, false); // print z_ij 
DEXTERN(bool, PRINT_LL, false); // print log likelihood 
DEXTERN(bool, PRINT_STARTS, false); // print starting subsequences 
DEXTERN(bool, NO_STATUS, false); // print run-time status to stderr 
DEXTERN(bool, TEXT_ONLY, false); // print documentation 
DEXTERN(bool, DOC, false); // print documentation 
EXTERN char *OFFSET_FILE; // current name of offset file 
// Type of timing:
//  0 : None
//  1 : Per loop
//  2 : Per start
DEXTERN(int, TIMER, 0);     

// macro to write a line of asterisks 
#define PSTARS(f) {int i; for (i=0;i<PAGEWIDTH;i++)fprintf(f, "*");fprintf(f, "\n");}

// type of negative motifs 
typedef enum {Pair, Blend} NEGTYPE;

// type of sequence to theta mapping 
typedef enum {Uni, Pam} MAP_TYPE;

// type of prior 
typedef enum {Mega, MegaP, Dmix, Dirichlet, Addone} PTYPE;

// type of handling of DNA strands in MAST 
typedef enum {Combine, Separate, Norc, Unstranded} STYPE;

// type of objective function 
typedef enum {Classic, NC, SmHG, CV, NZ, LL} OBJTYPE;
typedef enum {mHG, mRS} TESTTYPE;

// type of objective function for use during branching search 
typedef enum {LLR_POP, LLR, Likelihood, Class, GLAM} GLOB_OBJ_FN;

// type of point-wise branching moves to carry out 
typedef enum {NORMAL, X_ONLY, ALL, NO_POINT_B} POINT_BRANCHES;

// encapsulate how the [-m,...,-1,0,1,...,m] range of Z_i is stored
// The pointer to zi is offset to point to Z_i=0 
#define Zi(j) (zi[(j)])

// possible sites for a dataset, width combination 
#define ps(d, w) MAX((d)->n_samples,((d)->total_res-((d)->n_samples * ((w)-1))))

// tlb 6-18-99; added with wgt_total_res 
#define wps(d, w) ( MAX (wps1(d, w), 2) )
#define wps1(d, w) ( (d)->wgt_total_res - ((d)->n_samples * ( (w) - 1) ) )

// number of occurrences of a motif based on dataset, w, lambda 
#define nsites(d, w, l) ((l) * ps(d, w))

// scale min/max nsites based on dataset group
#define scale_min_nsites_to_group(d, g, min) \
  MAX(2, (int) (((double) (d)->n_group[g]/(d)->n_samples) * (min) + 0.5))
#define scale_max_nsites_to_group(d, g, max) \
  MIN(max, (int) (((double) (d)->n_group[g]/(d)->n_samples) * (max) + 0.5))

// number of independent columns in a model component 
#define ind_cols(w, pal) ((pal) ? ((w) + 1)/2 : (w))

// DNA palindrome enforcer: combine columns and complementary residues 
#define palindrome(theta1, theta2, w, alph)        \
{                   \
  int i, j;               \
  for (i=0; i<=(w+1)/2; i++) { /* create the palindrome */  \
    for (j=0; j< alph_size_core(alph); j++) {         \
      int ii = (w)-i-1;             \
      int jj = alph_complement((alph), (j));       \
      theta2[i][j] = (theta1[i][j] + theta1[ii][jj])/2;     \
      theta2[ii][jj] = MAX(0, theta2[i][j] - 1e-6);     \
    }                 \
  }                 \
}                 \

// Macros to control skipping of samples in different groups.
#define set_seq_groups_to_skip(dataset, skip0, skip1, skip2) \
  dataset->save_skip[0] = dataset->skip[0]; \
  dataset->save_skip[1] = dataset->skip[1]; \
  dataset->save_skip[2] = dataset->skip[2]; \
  dataset->skip[0] = skip0; \
  dataset->skip[1] = skip1; \
  dataset->skip[2] = skip2;

#define set_seq_groups_to_include(dataset, inc0, inc1, inc2) \
  set_seq_groups_to_skip(dataset, !(inc0), !(inc1), !(inc2))

#define restore_seq_groups_to_skip(dataset) \
  dataset->skip[0] = dataset->save_skip[0]; \
  dataset->skip[1] = dataset->save_skip[1]; \
  dataset->skip[2] = dataset->save_skip[2];

#define restore_seq_groups_to_include(dataset) \
  restore_seq_groups_to_skip(dataset)

#define skip_sample_if_required(dataset, sample) \
  if (dataset->skip[sample->group]) continue;

// dataset sample 
typedef struct sample {
  char *sample_name; // name of sample 
  char *descript; // description of sample 
  long length; // length of sequence 
  char *seq; // ascii sequence 
  uint8_t *res; // integer-coded sequence 
  uint8_t *resic; // integer-coded inverse complement 
  double sw; // sequence weight 
  double *weights; // Pr[pos not contained in previous site] 
  double *not_o; // P[no site in [x_{i,j}...x_{i,j+w}] 
  int *log_not_o; // log (not_o) 
  int **pY; // p(Y_j | theta_g) scratch spaces 
  char *pYic; // pY2 > pY1 (site on ic strand if != 0) 
  double *z; // tlb 6-21-99; E[z_ij] 
  double *z_buf; // memory allocated for z (as z can be offset)
  double *psp_original; // PSP for this sample; only used to set log_psp;
                        // derived for motif width: dataset->psp_w;
                        // NULL if using uniform priors in this sequence
  double *log_psp; // P_ij: log PSP for this sequence, 1 <= i <= length;
		   // P_i0: probability of not site in this sequence;
		   // normalized to width dataset->psp_current_w
  double max_log_psp; // always equal to largest value in log_psp for sequence
  double *log_psp_buf; // memory allocated for log_psp (as log_psp can be offset)
  double *counts; // counts of each character (X causes frac.) 
  double *logcumback; /*  log cumulative background probabilities
        logcumback[i] = 0 if i=1
                      = log Pr(s_{i-1} | H_0) otherwise
        */
  int nsites; // number of sites of current motif 
  int *sites; // list of sites of current motif 
  double minpv; // minimum p-value of sites of current motif 
  int group;	// sequence group: 0=main; 1=holdout1/noise1; 2=holdout2/noise2
  int orig_index;	// index of sample in input dataset (before shuffling)
} SAMPLE;

// Macro to compute the log background probability of 
// a subsequence from the log cumulative background probability array.
#define Log_back(lcb, i, w) (lcb[(i)+w] - lcb[i])

typedef double **THETA;
#define theta_ref(a, b, c)  (a)[b][c]
#define theta(b, c)   theta_ref(theta, b, c)
#define theta0(b, c)    theta_ref(theta0, b, c)
#define theta1(b, c)    theta_ref(theta1, b, c)
#define logtheta(b, c)    theta_ref(logtheta, b, c)
#define logtheta0(b, c)   theta_ref(logtheta0, b, c)
#define logtheta1(b, c)   theta_ref(logtheta1, b, c)
#define logtheta1_rc(b, c) theta_ref(logtheta1_rc, b, c)
#define obs(b, c)   theta_ref(obs, b, c)
#define obs1(b, c)    theta_ref(obs1, b, c)

// a site 
typedef struct p_prob *P_PROB;
typedef struct p_prob {
  int x; // sequence # 
  int y; // position # 
  bool ic; // on inverse complement strand 
  BOOLEAN negative;	// TRUE if from control dataset
  int rank;		// rank when sorted with control sites
  double ranksum;	// sum of ranks of primary sites
  double prob; // INT_LOG(probability of site) 
} p_prob;

// a model 
typedef struct Model {
  MOTYPE mtype; // type of model 
  int min_w; // minimum allowed width 
  int max_w; // maximum allowed width 
  bool all_widths; // consider all widths from min_w to max_w 
  double pw; // prior estimate of width 
  double min_nsites; // minimum allowed number of sites 
  double max_nsites; // maximum allowed number of sites 
  double psites; // prior estimate of number of sites 
  P_PROB maxima; // list of sites 
  bool pal; // motif is a DNA palindrome 
  bool invcomp; // use inverse complement DNA strand, too 
  int imotif; // number of motif 
  int w; // width of motif 
  THETA theta; // motif frequencies 
  THETA logtheta; // log of theta 
  THETA logtheta_rc; // log of reverse-complement of theta 
  THETA obs; // observed frequencies 
  double lambda; // lambda 
  double lambda_obs; // observed lambda 
  double nsites; // estimated number of sites 
  double nsites_obs; // observed number of sites 
  int nsites_dis; // number of sites after discretization 
  char cons[MAXSITE+1]; // consensus sequence of motif 
  char cons0[MAXSITE+1]; // character initial consensus 
  double rentropy[MAXSITE]; // relative entropy of each column of motif 
  double rel; // average relative entropy per col 
  double ic; // information content of motif 
  double ll; // log likelihood of all data under model 
  double mll_0; // motif log-likelihood under null model 
  double mll_1; // motif log-likelihood under motif model 
  double logpv; // log likelihood ratio of discrete motif 
  double logev; // log E-value of motif 
  double llr; // log likelihood ratio of motif 
  int iter; // number of EM iterations used 
  int ID; // processor id 
  int iseq; // start sequence 
  int ioff; // start sequence offset 
  double pc; // Performance coefficient
  OBJTYPE objfun; // objective function; needed for reduce_across_models()
  int alength; // length of the alphabet
} MODEL;

// user-input starting points 
typedef struct p_point {
  int c; // number of components 
  int w[MAXG]; // widths for motifs 
  double nsites[MAXG]; // nsites for motif 
  uint8_t *e_cons0[MAXG]; // integer encoded starting subsequence 
} P_POINT;

// starting point 
typedef struct s_point {
  double score; // log likelihood ratio of starting point 
  int iseq; // sequence number of starting point 
  int ioff; // sequence offset of starting point 
  int w0; // start width for motif 
  double nsites0; // start nsites0 for motif 
  double wgt_nsites; // effective (weighted) number of sites 
  uint8_t *e_cons0; // integer encoded starting subsequence 
  char *cons0; // character initial consensus 
  // This heap will contain the best seeds for this pair of
  // values, (w0, nsites0).
  HEAP *seed_heap;    
  // This indicates whether or not to store the result
  // at this s_point, when evaluating a seed under an
  // objective function.
  bool evaluate;   
  double sig; // significance of the LLR_POP "score" value 
} S_POINT;

// candidate final model 
typedef struct candidate {
  S_POINT *s_point; // starting point of model 
  int w; // final width of motif 
  bool pal; // palindrome flag 
  bool invcomp; // use inverse complement DNA strand, too 
  double lambda; // final lambda for motif 
  char cons[MAXSITE+1]; // final consensus of motif 
  double ic; // information content of motif 
  double rel; // relative entropy per col of each motif 
  double ll; // log-likelihood 
  double sig; // likelihood ratio test significance 
  double llr;           // LLR of model
} CANDIDATE;

 // summary of motif properties 
typedef struct motif_summary {
  int width; // width of motif 
  int num_sites; // num of sites of motif 
  int num_negative_sites; // num of sites of negative motif 
  double ic; // information content 
  double re; // relative entropy of motif 
  double llr; // log-likelihood ratio 
  double e_value_exp; // Exponent of E-value of motif 
  double e_value_mant; // Mantissa of E-value of motif 
  double bayes; // bayes threshold 
  double elapsed_time; // Time used to find motif 
  double **pssm; // log odds matrix 
  THETA psfm; // frequency matrix 
  char* regexp; // motif as regular expression 
  char* consensus; // motif as (possibly IUPAC) consensus
  P_PROB sites; // Pointer to array of motif sites 
} MOTIF_SUMMARY;

 // prior probabilities 
typedef struct Priors {
  PTYPE ptype; // type of prior to use 
  double prior_count[MAXALPH]; // ptype = Dirichlet: pseudo counts/letter 
  PriorLib *plib; // ptype = Dmix, Mega, MegaP: dirichlet mix 
  PriorLib *plib0; // ptype = MegaP: b=0 dirichlet mixture 
} PRIORS;

 // a known motif 
typedef struct motif {
  char name[MNAME]; // names of motif 
  int width; // (known) width of motif 
  int pos; // # positive samples this motif 
  double roc; // best roc for this motif 
  int shift; // best shift for this motif 
  int pass; // pass that found this motif 
  double recall; // best recall this motif 
  double precision; // best recall this motif 
  double min_thresh; // minimum threshold for motif 
  double max_thresh; // maximum threshold for motif 
  double pal; // motif is DNA palindrome 
  double like; // best likelihood this motif 
  double sig; // best significance this motif 
  double ic; // best info content this motif 
  double sites; // best nsites this motif 
  int w; // best width this motif 
  double thresh; // best threshold this motif 
  HASH_TABLE ht; // hash table of positives this motif 
} MOTIF;

// parameters controlling branching search 
typedef struct branch_params {
  int bfactor; // number of branches to perform 
  POINT_BRANCHES point_branch;  /* Controls what type of point branching
                                   (eg regular, X-only, none) to perform */
  bool w_branch; // Controls whether width-branching occurs 
} BRANCH_PARAMS;

// a dataset 
typedef struct Dataset {
 // set by read_seq_file 
  ALPH_T *alph; // alphabet 
  int total_res; // total size of dataset 
  double wgt_total_res; // weighted (sw*slen) total size of dataset
  int n_samples; // number samples 
  int n_group[3]; // number samples in groups 0, 1, 2
  SAMPLE **samples; // array of (pointers to) samples 
  SAMPLE **input_order; // samples in input order
  double *seq_weights; // array of sequence weights from input file
  int n_wgts; // number of sequence weights given
  long max_slength; // maximum length of sequences 
  long min_slength; // minimum length of sequences 
  int psp_w; // w defined by PSP file; 0 means no PSP file
  int log_psp_w; // log_psp is normalized for this width
  // 0 means need to normalize first time
  double *res_freq; // average letter frequencies 
  // *** MEME parameters ***
  // DNA palindrome flag:
  //  0 : no palindromes
  //  1 : force DNA palindromes
  bool pal;                  
  THETA map; // letter to frequency vector mapping 
  THETA lomap; // letter to logodds vector mapping 
  MOTIF motifs[NMOTIFS]; // known motifs in dataset 
  int nkmotifs; // number of known motifs in dataset 
  NEGTYPE negtype; // how to use negative examples 
  int back_order; // order of Markov background model 
  ARRAY_T *back; // Markov background model
  double log_total_prob; // total (log) cumulative background prob. 
  PRIORS *priors; // the prior probabilities model 
  P_POINT *p_point; // previously learned starting points 
  double wnsites; // weight on prior on nsites 
  bool ma_adj; // adjust width/pos. using mult. algn. method 
  double wg; // gap cost (initialization) 
  double ws; // space cost (extension) 
  bool endgaps; // penalize end gaps if true 
  double distance; // convergence radius 
  int nmotifs; // number of motifs to find 
  int maxiter; // maximum number of iterations for EM 
  double evt; // E-value threshold 
  char *mod; // name of model 
  char *mapname; // name of spmap 
  double map_scale; // scale of spmap 
  char *priorname; // name of type of prior 
  double beta; // beta for prior 
  int seed; // random seed 
  double seqfrac; // fraction of sequences to use 
  char *plib_name; // name of file with prior library 
  char *bfile; // name of background model file 
  char *datafile; // name of the dataset file 
  char *negfile; // name of negative examples file 
  char *command; // command line 
  OBJTYPE objfun; // objective function 
  THETA pairwise; // contains score matrix for pairwise scores 
  char *output_directory; // meme output directory 
  double max_time; // maximum allowed CPU time 
  int main_hs; // max size of heaps to be used throughout
  TESTTYPE test;		/* type of statistical test: Fisher (mHG) or rank-sum (mRS). */
  double ctfrac;		/* fraction of control sequences to use */
  int shuffle;			// preserve frequencies this size kmers when shuffling
  double hs_decrease;           /* the rate at which heap size decreases off
                                   central s_points */
  BRANCH_PARAMS *branch_params; // The branching params requested by the user 
  bool use_llr; // use true LLR fn during meme, cf POP 
  bool print_heaps; // print seed heap after each branch round 
  bool print_pllr; // print the LLR of the aligned planted sites 
  double pllr; // the LLR of the aligned planted sites 
  double param_V; // Parameter used by the "deme" objective fn 
  int kpos; // known position of motif in all sequences
  bool print_pred; // Print the sites predicted by MEME 
  char *correct_seed; // motif known to be in all sequences 
  double pseu[MAXALPH]; // letter pseudo-counts for use in GLAM fn 
  BOOLEAN skip[3];		// Which sequence groups to skip (0=main; 1,2: holdout/noise).
  BOOLEAN save_skip[3];		// Saves previous state of skip[].
  int last_seed_seq;		// Last sequence to search for seed words.
  int last_seed_pos;		// Last position to search for seed words.
  int max_words;		// Maximum number of seed words to test as EM starts.
  int imotif;                   /* # of the current motif being elucidated */
} DATASET;

// motif occurrence in sequence 
typedef struct {
  int seqno; // index in samples array 
  int pos; // character position of site in sequence 
  double zij; // value of z_ij 
  int invcomp; // on inverse complement strand 
} SITE;

// tiling of sequence with motifs 
typedef struct {
  // hits[i] = m, motif m occurs at position i in seq 
  // hits[i] = -m, motif m occurs at postion i on reverse strand of seq,
  // hits[i] = 0, means no hit
  int *hits;            
  double *pvalues; // pvalues[i] is p-value of match at position i 
  int  *svalues; // svalues[i] is scaled score of match at position i 
  double pv; // p-value of product of p-values of best hits 
  char *diagram;
} TILING;

 // subroutines
int meme_main(
  int argc,
  char *argv[],
  ARRAYLST_T* seq_array,
  int width,
  BOOLEAN_T eliminate_repeats
);

double exp(double x);
int em_subseq(
  THETA map, // freq x letter map 
  DATASET *dataset, // the dataset 
  MODEL *model, // the model 
  PRIORS *priors, // the priors 
  int w, // width of motif 
  int n_nsites0, // number of nsites0 values to try 
  double alpha, // sampling probability 
  P_POINT *p_point, // starting point for previous components 
  S_POINT s_points[] // array of starting points: 1 per nsites0 
);
void subseq7(
  MODEL *model, // the model
  DATASET *dataset, // the dataset 
  int w, // w to use
  int n_nsites0, // number of nsites0 values to try 
  S_POINT s_points[], // array of starting points: 1 per nsites0 
  // A hash table used for remembering which seeds have
  // been evaluated previously
  HASH_TABLE evaluated_seed_ht  
);
int pY_compare(
  const void *v1,
  const void *v2
);
void get_not_o(
  DATASET *dataset, // the dataset 
  int w // width of motif 
);
double get_log_sig(
  double llr, // log likelihood ratio 
  MOTYPE mtype, // type of model 
  int w, // width of motif 
  double wN, // weighted number of sites 
  int N, // number of sites 
  bool invcomp, // inv. compl. strand, too 
  bool pal, // motif is DNA palindrome 
  DATASET *dataset // the dataset 
);
void calc_entropy(
  MODEL *model, // the model 
  DATASET *dataset // the dataset 
);
double log_comb(
  int m, // length of sequence 
  int n // number of segments 
);
double get_log_nalign(
  MOTYPE mtype, // type of model 
  int w, // width of motif 
  int N, // number of occurrences 
  bool invcomp, // inv. compl. seq allowed 
  DATASET *dataset // the dataset 
);
void adjust_motif(
  MODEL *model, // the model 
  MODEL *scratch_model, // the scratch model 
  DATASET *dataset, // the dataset 
  PRIORS *priors, // the priors 
  double wnsites, // weight on prior on nsites 
  bool ma_adj, // adjust w using mult. algn. method  
  bool palindrome, // convert motif to palindrome 
  int c, // component of model to adjust 
  int min_w, // minimum width of motif allowed 
  int max_w, // maximum width of motif allowed 
  int maxiter, // maximum number iterations for EM 
  double distance, // stopping criterion 
  double wg, // gap cost (initialization) 
  double ws, // space cost (extension) 
  bool endgaps // penalize end gaps if true 
);
void init_theta(
  THETA theta, // theta 
  uint8_t *start, // integer encoded starting sequence 
  int w, // width of motif 
  THETA map, // frequency vectors for each letter 
  ALPH_T *alph // alphabet
);
S_POINT *get_starts(
  DATASET *dataset, // the dataset 
  MODEL *model, // the model 
  uint8_t *e_cons, // encoded consensus sequence 
  int *n_starts // number of starting points 
);
THETA init_map(
  // type of mapping:
  //  Uni - add n prior
  //  Pam - pam matrix
  MAP_TYPE type,    
  // degree of crispness, depends on type,
  //  Uni - add n prior (n)
  //  Pam - pam distance
  double scale,     
  ALPH_T *alph, // length of alphabet 
  ARRAY_T *back, // background frequencies 
  bool lo // create logodds matrix 
);
void convert_to_lmap (
  THETA map,
  int lmap[MAXALPH][MAXALPH],
  ALPH_T *alph
);
void convert_to_ltheta (
  double matrix_ds[MAXSITE][MAXALPH], // The input matrix of doubles
  int matrix_il[MAXSITE][MAXALPH], // The output matrix of int logs
  int nrows,
  int ncols
);
void copy_theta(
  THETA s, // source 
  THETA d, // destination 
  int w, // width of motif 
  int alength // length of alphabet 
);
void copy_model(
  MODEL *m1, // source 
  MODEL *m2, // destination 
  ALPH_T *alph // alphabet 
);
void free_model(
    MODEL *model // model 
);
void init_meme(
  int argc, // number of input arguments 
  char **argv, // input arguments 
  MODEL **model_p, // the model OUT 
  MODEL **best_model_p, // the best model OUT 
  MODEL **neg_model_p, // model of negative examples OUT 
  DATASET **dataset_p, // the dataset OUT 
  DATASET **neg_dataset_p, // dataset of negative examples OUT 
  char *text_filename, // name of the text output file IN 
  char **output_dirname, // name of the output directory OUT 
  FILE **text_output, // destination for text output OUT
  ARRAYLST_T* seq_array,
  int width,
  BOOLEAN_T eliminate_repeats
);

int exec_parallel(
  int argc,
  char *argv[]
);

MODEL *create_model(
  MOTYPE mtype, // type of model 
  bool invcomp, // use inv comp strand  
  int max_w, // maximum width of motif 
  ALPH_T *alph, // alphabet 
  OBJTYPE objfun // MEME objective function
);
double min_sites(
  double nu, // degrees of freedom 
  double alpha, // significance level 
  double max_h // maximum entropy 
);

int read_motifs (
  ALPH_T *alph,
  FILE *fdata, // opened dataset file 
  char *filename, // motif file 
  MOTIF motif[NMOTIFS], // motif info 
  bool save_dataset, // return dataset in memory 
  DATASET *dataset // integer-encoded dataset 
);

SITE *get_sites(
  DATASET *dataset, // the dataset 
  MODEL *model, // the model 
  int *n, // number of sites found 
  int *best_site // index of best site in array 
);

double **make_log_odds(
  THETA theta1, // motif theta 
  THETA theta0, // negative theta; use 0 if NULL 
  ARRAY_T *back, // background frequencies; use 0 if NULL 
  double q, // mixing parameter for theta0 and back 
  int w, // width of motif 
  int alength // length of alphabet 
);

int get_max(
  MOTYPE mtype, // the model type 
  DATASET *dataset, // the dataset 
  BOOLEAN negative,     // control dataset?
  int w, // length of sites 
  P_PROB maxima, // array of encoded site starts of local maxima 
  bool ic, // use inverse complement, too 
  bool sort // sort the maxima 
);

int align_top_subsequences(
  MOTYPE mtype, // type of model 
  int w, // width of motif 
  DATASET *dataset, // the dataset 
  int iseq, // sequence number of starting point 
  int ioff, // sequence offset of starting point 
  uint8_t *eseq, // integer encoded subsequence 
  char *name, // name of sequence 
  int n_nsites0, // number of nsites0 values to try 
  int n_maxima, // number of local maxima 
  P_PROB maxima, // sorted local maxima indices 
  double *col_scores, // column scores for last start point 
  S_POINT s_points[] // array of starting points 
);

double log_qfast(
  int n, // number of random variables in product 
  double logk // product of random variables 
);

void print_site_array(
  P_PROB sites, // An array of sites to be printed
  int nsites, // Length of the array
  FILE *outfile, // The stream for output
  int w, // The size of each of the sites
  DATASET *dataset // Contains the sequences which contain the sites
);

int get_first_siteloc (
  char *descript // The description of the sequence; contains site info
);

int get_n_strdiffs (
  char *str1, // A string in the aligned pair
  char *str2, // A string in the aligned pair
  int offset // Number of characters str1 is shifted to the right of
 // str2
);

void vector_subtract(
  int *col_a, // Positive column
  int *col_b, // Negative column
  int *diff_col, // Column to store the differences
  int col_len // Length of all columns
);

int *make_geometric_prog (
  int min_val, // Minimum integer value in the geometric progression
  int max_val, // Maximum integer value in the geometric progression
  double factor, // Factor specifying rate of increase in progression
  int *n_vals // The final number of values in the progression - OUT
);

int get_pred_sites (
  P_PROB psites, // An array to contain the predicted sites
  MOTYPE mtype, // Model type
  int w, // Length of seed being evaluated
  char *seed, // ASCII version of seed being evaluated
  int *lmotif[MAXSITE], // Storage space for the motif model
  int lmap[MAXALPH][MAXALPH], // The sequence to theta log map
  DATASET *dataset, // The dataset of sequences
  bool ic // Use inverse complement
);

void print_site_array(
  P_PROB sites, // An array of sites to be printed
  int nsites, // Length of the array
  FILE *outfile, // The stream for output
  int w, // The size of each of the sites
  DATASET *dataset // Contains the sequences which contain the sites
);

void renormalize (
  DATASET *dataset, // the dataset 
  int new_w, // new motif width 
  bool invcomp, // reverse complement? 
  MOTYPE mtype // OOPS, ZOOPS or TCM? 
);

#include "star.h"
#include "llr.h"
#include "em.h"
#include "read_seq_file.h"
#include "display.h"
#include "dpalign.h"
#include "histogram.h"
#include "message.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#endif
