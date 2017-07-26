/********************************************************************
 * FILE: ame.h
 * AUTHOR: Robert McLeay
 * CREATE DATE: 19/08/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * rsClover is a yet unpublished algorithm that seeks to assist in
 * determining whether a given transcription factor regulates a set
 * of genes that have been found by an experimental method that
 * provides ranks, such as microarray or ChIp-chip.
 *
 *
 * rsClover is a code name. The name will change prior to publication.
 ********************************************************************/
#ifndef __RSCLOVER_H__
#define __RSCLOVER_H__

#include <stdbool.h>
#include "html-monolith.h"
#include "motif.h"
#include "pssm.h"
#include "seq.h"

/*
 * ame-specific macros
 */

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a>b)?a:b

#define MAX_SEQ_LENGTH 1e6

/*
 * Define ame constants
 */

#define MEME_FORMAT 1   // input format is meme

#define UNIFORM_BG 0    //uniform.
#define MOTIF_BG 1      // background frequencies taken from motif file
#define FILE_BG 2       // background frequencies taken from specified file

#define QUICK_RS 0      //use my bodgy way of doing the ranksum
#define BETTER_RS 1     //use fabian's correct way of doing the ranksum

#define POS_FL 0  //use the fluorescence sort as positive indicator
#define POS_PWM 1 //use the PWM sort as positive indicator

#define RANKSUM_METHOD 0        //use the ranksum test
#define FISHER_METHOD  1        //use Fisher's exact test
#define MULTIHG_METHOD  2       //use Fisher's exact test modified to label with 0,1, or 2.
#define LONG_MULTIHG_METHOD  3  //use Fisher's exact test modified to label with 0,1,2, or 3.
#define LINREG_METHOD 4         //use linear regression test to minimise MSE.
#define SPEARMAN_METHOD 5       //use the spearman rank correlation co-efficient to calculate a score.

/*
 * Struct definitions
 */

typedef struct  {
  char* alph_filename; // file containing xalph
  char* bg_filename; //file to get base freqs from
  char* outputdir; // where to send outputs
  int bg_format; //whether it's fasta, meme, etc.
  char** motif_filenames; //filenames of the motif library
  int number_motif_files; //how many files comprise motif library
  char* sequence_filename; //(primary) input sequences in fasta format
  char* control_filename; //control input sequences in fasta format
  char* commandline; // command line with path stripped from program name
  float pseudocount; //add to the motif frequency counts
  int scoring; //AVG_ODDS or MAX_ODDS
  int verbose;
  int rs_method; //QUICK_RS or BETTER_RS
  int positive_list; //POS_FL or POS_PWM
  int pvalue_method; //RANKSUM_METHOD or FISHER_METHOD
  double fisher_pwm_threshold;
  double fisher_fasta_threshold;
  double pvalue_threshold; //threshold for the mhg test with pwms.
  double pvalue_report_threshold; //threshold for reporting a motif in output.
  bool length_correction;
  bool log_fscores;
  bool log_pwmscores;
  bool linreg_normalise;
  bool linreg_switchxy;
  bool clobber; // TRUE if we can replace existing output directory
  bool silent; // suppress version and command line output; for testing: discourage in normal use
  char* linreg_dump_dir;
  int fix_partition;
  // derived from command line:
  FILE *text_output;
  HTMLWR_T *html_output;
  JSONWR_T *json_output;
} rsc_arg_t;

typedef struct {
  char* filename;
  BOOLEAN_T revcomp;	// reverse complements included.
  MOTIF_T* motifs;	// If revcomp==TRUE, every other
			// index contains the RC copy of the 
			// preceding motif.
  PSSM_T* pssms;
  MATRIX_T** odds;
  int* db_idx;
  ARRAY_T** pv_lookup;
  int num_b4_rc;	// Number of motifs (not counting 
			// reverse complements).
  ARRAY_T* bg_freqs;
} rsc_motifs_t;

typedef struct {
  int db_idx;
  MOTIF_T *motif;
  int split;
  double pleft;
  double pright;
  double pboth;
  double u;
} rsc_result_t;

typedef struct {
  int f_rank;
  double f_score;
  int pwm_rank;
  double pwm_score;
} rsc_rank_t;


/*
 *  Function definitions
 */

void rsc_getopt(int argc, char *argv[]);
const char* rsc_get_usage();
ARRAY_T* rsc_load_background();
void rsc_load_motifs();
void rsc_scan_sequences();
void rsc_get_scores();
void rsc_usage();
void rsc_terminate(int status);
int rsc_compare_doubles (const void *a, const void *b);
int rsc_compare_scores (const void *a, const void *b);
int rsc_compare_mse (const void *a, const void *b);
int rsc_compare_ranks_f_rank (const void *a, const void *b);
int rsc_compare_ranks_pwm_score (const void *a, const void *b);
void rsc_dump_motif(MOTIF_T* m);
rsc_result_t* rsc_do_ranksum_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_fisher_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_multihg_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_linreg_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_spearman_test(rsc_rank_t** rankings, int motif_index);
unsigned long long choose(unsigned n, unsigned k);
long double logchoose(unsigned n, unsigned k);
double* rsc_init_fisher_factorials(int len);
double rsc_bonferroni_correction(double, double);

double rsc_score_sequence(MOTIF_T* motif, SEQ_T* seq, int scoring);

#else
#endif
