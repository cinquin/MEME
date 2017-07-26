/************************************************************************
*                                                                       *
*       MAST                                                            *
*       Author: Timothy L. Bailey                                       *
*                                                                       *
*       Copyright                                                       *
*       (1994 - 2001) The Regents of the University of California.      *
*       All Rights Reserved.                                            *
*                                                                       *
*       Permission to use, copy, modify, and distribute any part of     *
*       this software for educational, research and non-profit purposes,*
*       without fee, and without a written agreement is hereby granted, *
*       provided that the above copyright notice, this paragraph and    *
*       the following three paragraphs appear in all copies.            *
*                                                                       *
*       Those desiring to incorporate this software into commercial     *
*       products or use for commercial purposes should contact the      *
*       Technology Transfer Office, University of California, San Diego,*
*       9500 Gilman Drive, La Jolla, California, 92093-0910,            *
*       Ph: (619) 534 5815.                                             *
*                                                                       *
*       IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO     *
*       ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR         *
*       CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF   *
*       THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA  *
*       HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
*                                                                       *
*       THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*       UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE          *
*       MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*       THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND       *
*       EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*       INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        *
*       MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT    *
*       THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,           *
*       TRADEMARK OR OTHER RIGHTS.                                      *
************************************************************************/

/**********************************************************************/
/*
  mast <memefile> <database> [optional arguments...]
        See <install-path>/bin/mast for documentation.
*/
/**********************************************************************/

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
// for floating point exceptions
//#include <fenv.h>

#define DEFINE_GLOBALS
#include "mast.h"

#include "array.h"
#include "array-list.h"
#include "buffer.h"
#include "dir.h"
#include "io.h"
#include "motif-in.h"
#include "alphabet.h"
#include "diagram.h"
#include "projrel.h"
#include "read_sequence.h"
#include "string-builder.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

#define BAD_MOTIF_MARK 0

 // error flag 
bool errored = false;
 // printing 
VERBOSE_T verbosity = QUIET_VERBOSE;

const char *program_name = "mast"; // the program name 
static char *default_out_dir = "mast_out"; // default name of output 
static const char *XML_FILENAME = "mast.xml";
static const char *HTML_FILENAME = "mast.html";
static const char *TXT_FILENAME = "mast.txt";
static const char *MAST2TXT_FILENAME = "mast2txt";
    
typedef struct MAST_OPTIONS {
  time_t start_time; // the start time
  char *out_dir; // the output directory
  bool clobber; // should the output directory be clobbered
  char *bfile; // background file
  char *mo_source; // the file path or "-" indicating stdin. The file will be MEME motifs
  char *db_source; // the file path or "-" indicating stdin. The file will be either a fasta file or a list of fasta files
  bool db_list; // treat the db_source as a file containing a list of db files
  char *mo_name; // when not null, use to replace motif file name in results
  char *db_name; // when not null, use to replace sequence file name in results
  char *db_link; // when not null, use to create links for searching sequence names in results
  int first_motifs; // read only the first x motifs 
  double e_max; // maximum sequence p-value to output 
  double m_thresh; // maximum motif p-value to output 
  double w_thresh; // maximum motif p-value--weak hits; default: m_thresh * 10
  double min_ev; // minimum E-value of motifs to use 
  int min_seqs; // lower bound on total sequence count, allows for earlier filtering of results with a large evalue
  bool adj_hit_pvalue; // use sequence p-value for m_thresh 
  bool shuffle; // shuffle the motif columns
  bool sonly; // print product of spacing p-values 
  bool lump; // combine spacings into one p-value 
  bool status; // show progress 
  bool html; // generate html output 
  bool mast2txt; // run mast2txt on the output 
  STYPE stype; // handling of DNA strands 
  bool translate_dna; // when true expect the motifs to be protein and the sequences to be DNA
  bool best_motifs; // only print best motif in diagrams 
  bool use_seq_comp; // adjust E-values for actual compos. 
  bool rem_corr; // remove overly correlated motifs  
  bool hit_list; // print block diagrams  
  mt_state prng; // random number generator for motif column shuffling
  uint32_t seed; // random number generator seed
  RBTREE_T *ids; // user selected IDs
  RBTREE_T *idxs; // user selected indexes
  char *diagram;
} MAST_OPTIONS_T;

 // MAST DTD 
 // {{{
char* MAST_DTD = 
"<!DOCTYPE mast[\n"
"<!ELEMENT mast (model, alphabet, motifs, sequences, runtime)>\n"
"<!ATTLIST mast version CDATA #REQUIRED release CDATA #REQUIRED>\n"
"<!ELEMENT model (command_line, max_correlation, remove_correlated, strand_handling, translate_dna, max_seq_evalue,\n"
"    adj_hit_pvalue, max_hit_pvalue, max_weak_pvalue, host, when)>\n"
"<!ELEMENT command_line (#PCDATA)>\n"
"<!ELEMENT max_correlation (#PCDATA)>\n"
"<!ELEMENT remove_correlated EMPTY>\n"
"<!ATTLIST remove_correlated value (y|n) #REQUIRED>\n"
"<!ELEMENT strand_handling EMPTY>\n"
"<!ATTLIST strand_handling value (combine|separate|norc|protein) #REQUIRED>\n"
"<!ELEMENT translate_dna EMPTY>\n"
"<!ATTLIST translate_dna value (y|n) #REQUIRED>\n"
"<!ELEMENT max_seq_evalue (#PCDATA)>\n"
"<!ELEMENT adj_hit_pvalue EMPTY>\n"
"<!ATTLIST adj_hit_pvalue value (y|n) #REQUIRED>\n"
"<!ELEMENT max_hit_pvalue (#PCDATA)>\n"
"<!ELEMENT max_weak_pvalue (#PCDATA)>\n"
"<!ELEMENT host (#PCDATA)>\n"
"<!ELEMENT when (#PCDATA)>\n"
"<!ELEMENT alphabet (letter*)>\n"
"<!ATTLIST alphabet type (amino-acid|nucleotide) #REQUIRED bg_source (preset|file|sequence_composition) #REQUIRED bg_file CDATA #IMPLIED>\n"
"<!ELEMENT letter EMPTY>\n"
"<!ATTLIST letter symbol CDATA #REQUIRED ambig (y|n) \"n\" bg_value CDATA #IMPLIED>\n"
"<!ELEMENT motifs (motif*,correlation*,nos*)>\n"
"<!ATTLIST motifs source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED>\n"
"<!ELEMENT motif EMPTY>\n"
"<!-- num is simply the loading order of the motif, it's superfluous but makes things easier for XSLT -->\n"
"<!ATTLIST motif id ID #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED width CDATA #REQUIRED\n"
"   best_f CDATA #REQUIRED best_r CDATA #IMPLIED bad (y|n) \"n\">\n"
"<!-- for n > 1 motifs there should be (n * (n - 1)) / 2 correlations, obviously there are none for only 1 motif -->\n"
"<!ELEMENT correlation EMPTY>\n"
"<!ATTLIST correlation motif_a IDREF #REQUIRED motif_b IDREF #REQUIRED value CDATA #REQUIRED>\n"
"<!-- nos: Nominal Order and Spacing diagram, a rarely used feature where mast can adjust pvalues for an expected motif spacing -->\n"
"<!ELEMENT nos (expect*)>\n"
"<!-- length is in the same unit as the motifs, which is not always the same unit as the sequence -->\n"
"<!ATTLIST nos length CDATA #REQUIRED>\n"
"<!-- the expect tags are expected to be ordered by pos ascending -->\n"
"<!ELEMENT expect EMPTY>\n"
"<!ATTLIST expect pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED>\n"
"<!ELEMENT sequences (database*, sequence*)>\n"
"<!-- the database tags are expected to be ordered in file specification order -->\n"
"<!ELEMENT database EMPTY>\n"
"<!ATTLIST database id ID #REQUIRED num CDATA #REQUIRED source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED \n"
"    seq_count CDATA #REQUIRED residue_count CDATA #REQUIRED type (amino-acid|nucleotide) #REQUIRED link CDATA #IMPLIED>\n"
"<!-- the sequence tags are expected to be ordered by best combined p-value (of contained score tags) ascending -->\n"
"<!ELEMENT sequence (score*,seg*)>\n"
"<!ATTLIST sequence id ID #REQUIRED db IDREF #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED comment CDATA \"\" length CDATA #REQUIRED>\n"
"<!ELEMENT score EMPTY>\n"
"<!-- frame is the starting offset for translation of dna sequences which gives the lowest pvalues for the provided protein motifs -->\n"
"<!ATTLIST score strand (both|forward|reverse) #REQUIRED frame (a|b|c) #IMPLIED combined_pvalue CDATA #REQUIRED evalue CDATA #REQUIRED>\n"
"<!-- within each sequence the seg tags are expected to be ordered by start ascending -->\n"
"<!ELEMENT seg (data,hit*)>\n"
"<!ATTLIST seg start CDATA #REQUIRED>\n"
"<!ELEMENT data (#PCDATA)>\n"
"<!-- within each seg the hit tags are expected to be ordered by pos ascending and then forward strand first -->\n"
"<!ELEMENT hit EMPTY>\n"
"<!-- gap, while superfluous, makes creating motif diagrams for the text version much easier when using XSLT -->\n"
"<!ATTLIST hit pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED pvalue CDATA #REQUIRED strand (forward|reverse) \"forward\" \n"
"    match CDATA #REQUIRED translation CDATA #IMPLIED>\n"
"<!ELEMENT runtime EMPTY>\n"
"<!ATTLIST runtime cycles CDATA #REQUIRED seconds CDATA #REQUIRED>\n"
"]>\n";
 // }}}

/* these two macros are used to convert a preprocessor define into a string
 * you would logicaly think that it would work with just the one macro but
 * for some odd reason it doesn't */
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

#ifndef EXP
#define EXP 0
#else
#define EXP 1
#endif

 // default e-thresh 
#define EXPECT 10

#define MOTIF_THRESH_DEFAULT 1e-4

 // maximum pairwise motif correlation 
#define MAXCORR 0.6

 // number of different integral score values 
#define MAST_RANGE 100

 // the smallest size of a segment in the output 
#define SEG_CHUNK 75
#define SEG_FORMAT "%." STRINGIZE(SEG_CHUNK) "s\n"

/*
  Data Definitions (structs)
*/

typedef struct database_t DATABASE_T;
typedef struct strand_t STRAND_T;
typedef struct sseq_t SSEQ_T;

 // {{{
struct database_t {
  int index; // loading number of the database 
  char *source; // the source of the database "-" for stdin otherwise the file name 
  char *name; // the display name of the database 
  time_t last_mod; // the unix time the database was last modified 
  ALPH_T *alph; // the alphabet of the database
  int sequence_count; // the count of sequences in this database 
  long residue_count; // the count of residues in this database 
  char *link; // the link to search for further information on sequences 
  FILE *file; // the file to read from, may be stdin 
  FILE *save; // an temporary file to save sequences read from stdin 
};

 // sortable sequence score record 
struct strand_t {
  SSEQ_T *sequence; // details of the sequence 
  double Pvalue; // p-value of product of p-values 
  int strand; // -1 neg. strand, 0 both/protein, +1 pos. strand 
  double *best_scores; // the best score for each motif 
  int *best_location; // the location of the best scores for each motif 
};

struct sseq_t {
  int db_index; // identifies the database 
  int index; // loading number of the sequence 
  long  fp; // file pointer to beginning of sequence rec. 
  long length; // length of sequence 
  ARRAY_T *comp; // actual sequence composition or NULL if  
  int pv_alloc_len; // what is the length of the allocated pv matrix, or 0 if a reference which we don't need to free
  double **pv; // pvalue distribution for each motif (may be allocated or a reference)
  STRAND_T *pos_strand; // forward strand 
  STRAND_T *neg_strand; // reverse strand 
};
 // }}}

/************************************************************************/
/*
        error

        Prints an error message to stderr and sets a flag to exit on the next
        call to exit_on_error
*/
/************************************************************************/
static void error(char* format, ...) {
  va_list  argp;

  fprintf(stderr, "FATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  errored = true;
}

/************************************************************************/
/*
        exit_on_error

        Checks if the function error has been called and quits if it has
*/
/************************************************************************/
static void exit_on_error() {
  if (errored) {
#ifdef DEBUG
    abort();
#else
    exit(1);
#endif
  }
}

/**********************************************************************/
/*
        calc_p_value

        Calculate the p-value of the product of the individual motif 
        p-values < observed

        Returns the combined p-value.
*/
/**********************************************************************/
static double calc_p_value(
  int sample_index, // index of sample 
  long length, // length of the sequence 
  STYPE stype, // handling of DNA strands 
  int nmotifs, // number of motifs 
  LO *los[], // array of pointers to log-odds matrices 
  double *best_score, // array of best score for each motif 
  int *best_loc, // position of best match for each motif 
  double range, // number of different score values 
  double **pv, // p-value tables for each motif 
  bool sonly, // calculate p-value of observed spacings only 
  bool lump, // combine spacings into one p-value 
  int norder // number of motifs in diag 
)
{
  int i;
  double pvalue;
  int nspaces; // number of motif pairs spacings given for 
  double k_scores; // product of motif score p-values 
  double k_spacing; // product of spacing p-values 
  int smotifs; // number of motifs that got scored 

  // calculate the product of the motif score p-values 
  k_scores = 1.0;
  for (i=smotifs=0; i<nmotifs; i++) {
    int ws = los[i]->ws; // width of motif in sequence 
    int x = best_score[i]; // observed EV (extreme value) 
    long n = length - ws + 1; // number of samples in EV 
    double p; // p-value of EV 
    if (best_score[i] == LITTLE) continue; // motif wasn't scored 
    if (stype == Combine) n*=2; // combining both DNA strands 
    EV(pv[i][x], n, p); // compute sequence p-value of x 
    k_scores *= p; // product of p-values 
    smotifs++; // number of motifs scored 
  }

  // multiply by p-values of motif spacings if provided (norder > 1) 
  k_spacing = 1.0;
  nspaces = 0; // number of spaces p-valued 
  for (i=1; i<norder; i++) { // space to previous motif 
    if (space[i] >= 0) { // don't ignore this space 
      double p; // to hold p-value of obs. spacing 
      int err; // error in pos ith and i-1th motifs 
      int mi = order[i]; // current motif 
      int mim1 = order[i-1]; // previous motif  
      long npos; // number of positions for motif 
      int ws_avg; // average width of adjacent motifs 
      long mind; // minimum spacing that fits seq 
      long maxd; // maximum spacing that fits seq 

      // get absolute error: difference between observed and nominal spacing 
      err = abs(space[i] - (best_loc[mi] - best_loc[mim1]));

      // get the average motif width and maximum number of placements 
      ws_avg = (int) (los[mim1]->ws + los[mi]->ws)/2.0; // round down 
      npos = length - ws_avg + 1;

      mind = MAX(space[i] - err, ws_avg - length);
      maxd = MIN(space[i] + err, length - ws_avg);
      if (mind >= 0) {
        p = (maxd-mind+1) * (npos - (maxd+mind)/2.0);
      } else {
        p = -mind * (npos - (1-mind)/2.0) + (maxd+1) * (npos - maxd/2.0);
      }
      p /= npos * npos;
      if (p > 1.0 || p < 0.0) fprintf(stderr, 
        "\nerror in spacing p-value:seq_%d L=%8ld mind=%8ld E=%4d %-10g\n", 
         sample_index, length, mind, err, p);

      // skip if no possible spacing 
      if (maxd < mind) {
        continue;
      }

      k_spacing *= p; // product of p-values 
      nspaces++; // number of spacings used 

    } // don't ignore 
  } // space to previous motif 

 // finish the calculation 
  if (sonly) { // spacings only 
    pvalue = qfast(nspaces, k_spacing);
  } else if (lump && nspaces) { // lump spacings 
    pvalue = qfast(nspaces, k_spacing); // spacing pvalue 
    pvalue = qfast(smotifs+1, k_scores*pvalue); // spacing and motif p-value 
  } else { // spacings and motif scores 
    pvalue = qfast(smotifs+nspaces, k_scores*k_spacing);
  }

  return pvalue;
} // calc_p_value 

/**********************************************************************/
/*
        best_frame

        Find the predominant frame of the matches to the motifs for
        a sequence.  The best frame is the frame with the minimum
        product of p-values.

        Returns an character a, b or c
*/
/**********************************************************************/
static char best_frame(
  XLATE_T *xlate,
  int nmotifs, // number of motifs 
  long length, // length of sequence 
  TILING tiling // tiling of sequence 
) {
  double best_frame_p, prod_p;
  double *best_p;
  int i, m, best;
  long j;

  if (xlate_src_nsyms(xlate) > 26) {
    die("Something has gone wrong, the sequence translator reports that it "
        "converts more than 26 symbols, this number should be around 3...");
    return '?';
  }

  best_frame_p = 1; // p-value of best frame 
  best = 0; // best frame 

  // find the product of p-values for each frame 
  best_p = mm_malloc(sizeof(double) * nmotifs);
  for (i = 0; i < xlate_src_nsyms(xlate); i++) { // frame 
    // reset best p-values for each motif 
    for (j = 0; j < nmotifs; j++) best_p[j] = 1; // best p-value for motif j 
    prod_p = 1; // product of p-values 

    // find best p-value for each motif in this frame 
    for (j = i; j < length; j += xlate_src_nsyms(xlate)) { // position in sequence 
      if ((m = abs(tiling.hits[j]) - 1) >= 0) {
        best_p[m] = MIN(best_p[m], tiling.pvalues[j]);
      }
    }

    // form the product of p-values for this frame 
    for (j = 0; j < nmotifs; j++) prod_p *= best_p[j];

    // update best frame 
    if (prod_p < best_frame_p) {
      best_frame_p = prod_p;
      best = i;
    }
  } // frame 
  free(best_p);
  best_p = NULL;

  assert(best < 26);
  assert(best < xlate_src_nsyms(xlate));

  return (char)(((int)'a') + best);
} // best_frame 


/**********************************************************************/
/* 
 * best_scores
 * 
 * scores the sequence and pulls out the best scores for
 * each motif.
 * If there are existing scores in the best_score array it will consider
 * them in determining the best score.
 */
/**********************************************************************/
static void best_scores(
  ALPH_T *alph,
  XLATE_T *xlate, // database is DNA and motifs protein 
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs, // number of motifs 
  char *sequence, // sequence 
  long length, // length of the sequence 
  bool rc_seq, // should the scores be found for the reverse complment 
  double **best_score, // OUT: the best scores for each motif, caller must deallocate 
  int **best_loc // OUT: the location for each best score, caller must deallocate 
  ) 
{
  SCORE **scores;
  double *score; // best score per motif for positive strand 
  int *loc, i, j, last_j; // best location per motif for positive strand 

  score = *best_score;
  loc = *best_loc;
  if (score == NULL) {
    Resize(score, nmotifs, double);
    for (i = 0; i < nmotifs; i++) score[i] = LITTLE;
  }
  if (loc == NULL) {
    Resize(loc, nmotifs, int);
    for (i = 0; i < nmotifs; i++) loc[i] = 0;
  }

  scores = score_it(alph, xlate, rc_seq, los, nmotifs, sequence, length);
  

  for (i = 0; i < nmotifs; i++) { //motifs
    last_j = length - los[i]->ws;
    for (j = 0; j <= last_j; j++) {
      if (scores[i][j].score > score[i]) {
        score[i] = scores[i][j].score;
        loc[i] = j;
      }
    }
  }
  free_2array(scores, nmotifs);

  *best_score = score;
  *best_loc = loc;
}

/**********************************************************************/
/*
 * strand_destroy
 *
 * Destroys a strand score structure
 */
/**********************************************************************/
static void strand_destroy(STRAND_T *strand) {
    if (strand->best_scores) myfree(strand->best_scores);
    if (strand->best_location) myfree(strand->best_location);
    memset(strand, 0, sizeof(STRAND_T));
    myfree(strand);
}

/**********************************************************************/
/*
 * sseq_destroy
 *
 * Destroys a scored sequence
 * Note, also destroys connected
 * strand score structures.
 */
/**********************************************************************/
static void sseq_destroy(void *p) {
  SSEQ_T *sseq;
  sseq = (SSEQ_T*)p;
  // deallocate connected strand objects
  if (sseq->neg_strand) strand_destroy(sseq->neg_strand);
  if (sseq->pos_strand) strand_destroy(sseq->pos_strand);
  // free up optional allocations
  if (sseq->comp) free_array(sseq->comp);
  //free expected allocations
  if (sseq->pv_alloc_len) {
    free_2array(sseq->pv, sseq->pv_alloc_len);
  }
  //zero memory
  memset(sseq, 0, sizeof(SSEQ_T));
  //free struct
  myfree(sseq);
}

/**********************************************************************/
/*
 * add_strand
 *
 * Helper function to add a score for a strand (or both strands) to 
 * the list of scores
 */
/**********************************************************************/
static void add_strand(ARRAYLST_T *strand_scores, SSEQ_T *sseq, int strand, double pvalue, int nmotifs, double *best_score, int *best_loc) {
  STRAND_T *score = NULL;
  Resize(score, 1, STRAND_T);
  score->sequence = sseq;
  score->Pvalue = pvalue;
  score->strand = strand;
  if (sseq->comp) { // only needed if adjusting for sequence composition
    score->best_scores = mymalloc(nmotifs * sizeof(double));
    memcpy(score->best_scores, best_score, nmotifs * sizeof(double));
    score->best_location = mymalloc(nmotifs * sizeof(int));
    memcpy(score->best_location, best_loc, nmotifs * sizeof(int));
  } else {
    score->best_scores = NULL;
    score->best_location = NULL;
  }
  if (strand != -1) {
    sseq->pos_strand = score;
  } else {
    sseq->neg_strand = score;
  }
  arraylst_add(score, strand_scores);
}

/**********************************************************************/
/*
        get_scores

        Calculate the score and p-value for each sequence in the
        database. 

        Returns scores of sequences and sequence info in array lists
*/
/**********************************************************************/
static void get_scores(
  MAST_OPTIONS_T *opts,
  ALPH_T *alph,
  XLATE_T *xlate,
  DATABASE_T *database, // database 
  LO *los[], // array of pointers to log-odds matrices 
  int nmotifs, // number of log-odds matrices in los 
  double **pv, // p-value tables for each motif 
  int *nseqs, // total sequences in all databases 
  double *residues, // total number of residues in all seqs 
  ARRAYLST_T *sequences, // IN/OUT: sequences saved because they have interesting scores 
  ARRAYLST_T *strand_scores // IN/OUT: interesting strand scores 
)
{
  int i, j, imotif, sample_index;
  int est_seq; // estimated number of sequences
  long sample_len; // length of sample
  long file_pos; // needs to be long for files larger than 4Gb
  char *sample_name, *sample_seq, *sample_comment;
  
  double *best_score = NULL; // best score per motif for positive strand 
  int *best_loc = NULL; // best location per motif for positive strand 

  STRAND_T *seqs = NULL; // sortable array of seq/score records
  double pvalue; // expected number of sequences / total >= x 
  SCORE **scores; // scores for each motif vs seq. pos. 
  bool saved = false; // sequence already saved 
  SSEQ_T *sseq;

  sample_index = 0;
  // loop through each sequence
  file_pos = (database->save != NULL) ? ftell(database->save) : ftell(database->file);
  while (read_sequence(alph, database->file, &sample_name, &sample_comment, &sample_seq, &sample_len)) {
    // Skip over special MEME WEIGHTS extension
    if (sample_len == 0 && strcmp(sample_name, "WEIGHTS") == 0) {
      myfree(sample_seq);
      myfree(sample_comment);
      myfree(sample_name);
      file_pos = (database->save != NULL) ? ftell(database->save) : ftell(database->file);
      continue;
    }
    sseq = (SSEQ_T*)mymalloc(sizeof(SSEQ_T));
    sseq->db_index = database->index;
    sseq->index = sample_index;
    sseq->length = sample_len;
    sseq->fp = file_pos;
    sseq->comp = opts->use_seq_comp ? get_seq_comp(alph, xlate, sample_seq) : NULL;
    sseq->pv = pv;
    sseq->pv_alloc_len = 0;
    sseq->pos_strand = NULL;
    sseq->neg_strand = NULL;
    saved = false; // new sequence 

    // update size of database 
    if (opts->stype == Separate) { // treat each DNA sequence as two 
      (*nseqs) += 2; // number of sequences in database 
      *residues += 2*sample_len; // number of residues in database 
    } else if (opts->stype != Separate) { // treat each sequence as one 
      (*nseqs)++; // number of sequences in database 
      *residues += sample_len; // number of residues in database 
    }
    database->sequence_count += 1;
    database->residue_count += sample_len;

    //clear the scores prior to use
    if (best_score) for (i = 0; i < nmotifs; ++i) best_score[i] = LITTLE;
    if (best_loc) for (i = 0; i < nmotifs; ++i) best_loc[i] = 0;

    est_seq = MAX(opts->min_seqs, *nseqs);//get the current estimate sequence number

    //score the positive strand
    best_scores(alph, xlate, los, nmotifs, sample_seq, sseq->length, false, &best_score, &best_loc);
    if (opts->stype == Combine) {//and the negative strand in combined mode
      best_scores(alph, xlate, los, nmotifs, sample_seq, sseq->length, true, &best_score, &best_loc);
    }

    // calculate the combined p-value for the positive (or maybe both) strand(s)
    pvalue = calc_p_value(sseq->index, sseq->length, opts->stype, nmotifs, los, best_score, 
        best_loc, MAST_RANGE, pv, opts->sonly, opts->lump, norder);
    
    if (pvalue * est_seq < opts->e_max) {//save the pos strand
      saved = true;
      add_strand(strand_scores, sseq, (opts->stype == Combine || opts->stype == Unstranded ? 0 : 1), pvalue, nmotifs, best_score, best_loc);
    }
    

    if (opts->stype == Separate) {
      //clear the scores prior to use
      if (best_score) for (i = 0; i < nmotifs; ++i) best_score[i] = LITTLE;
      if (best_loc) for (i = 0; i < nmotifs; ++i) best_loc[i] = 0;
      //score negative strand
      best_scores(alph, xlate, los, nmotifs, sample_seq, sseq->length, true, &best_score, &best_loc);
      // calculate the combined p-value for the negative strand
      pvalue = calc_p_value(sseq->index, sseq->length, opts->stype, nmotifs, los, best_score, 
          best_loc, MAST_RANGE, pv, opts->sonly, opts->lump, norder);
      
      if (pvalue * est_seq < opts->e_max) {//save the neg strand
        saved = true;
        add_strand(strand_scores, sseq, -1, pvalue, nmotifs, best_score, best_loc);
      }
    }

    // save the sequence if a strand has an E-value that may be under threshold 
    if (saved) { 
      arraylst_add(sseq, sequences);
      // write the sequence to the save file if reading standard input
      if (database->file == stdin) {
        fprintf(database->save, ">%s %s\n%s\n", sample_name, sample_comment, sample_seq);
      }
    }

    // print progress report
    if (opts->status && *nseqs % 100 == 0) 
      fprintf(stderr, "\rsequences: %6d ", *nseqs);

    // free up space
    if (!saved) sseq_destroy(sseq);
    myfree(sample_seq);
    myfree(sample_comment);
    myfree(sample_name);

    file_pos = (database->save != NULL) ? ftell(database->save) : ftell(database->file);
    sample_index++;
  } // read_sequence 

  // free up space
  myfree(best_score);
  myfree(best_loc);

  /*
    recalculate the p-values based on actual sequence composition
  */
  for (i=0; opts->use_seq_comp && i < arraylst_size(sequences); i++) {
    double **pv;
    sseq = arraylst_get(i, sequences);

    //skip if E-value too big, or no scores for this sequence
    if ((sseq->pos_strand == NULL || (sseq->pos_strand->Pvalue * (*nseqs) > opts->e_max)) && 
        (sseq->neg_strand == NULL || (sseq->neg_strand->Pvalue * (*nseqs) > opts->e_max))) continue;       

    // calculate p-values of all integer score values in range [0...w*MAST_RANGE] 
    pv = (double**)mymalloc(nmotifs * sizeof(double*));
    for (j=0; j<nmotifs; j++) {
      pv[j] = calc_pssm_cdf(los[j]->w, los[j]->alen, MAST_RANGE, los[j]->logodds, sseq->comp);
    }

    // recalculate the combined p-value for this sequence using the actual
    // sequence composition as the background model
    if (sseq->pos_strand) {// for the positive strand
      STRAND_T *st = sseq->pos_strand;
      st->Pvalue = calc_p_value(sseq->index, sseq->length, opts->stype, nmotifs, los, st->best_scores, 
          st->best_location, MAST_RANGE, pv, opts->sonly, opts->lump, norder);
#ifdef DEBUG
      if (isnan(st->Pvalue)) {
        fprintf(stderr, "nan: i %d d-ind %d s-ind %d strand %d pvalue %f\n", i, sseq->db_index, 
            sseq->index, st->strand, st->Pvalue); 
        for (j=0; j<alen; j++) fprintf(stderr, "%c %f\n", unhash(j), get_array_item(j, sseq->comp));
        abort();
      }
#endif
    }
    if (sseq->neg_strand) { // for the negative strand
      STRAND_T *st = sseq->neg_strand;
      st->Pvalue = calc_p_value(sseq->index, sseq->length, opts->stype, nmotifs, los, st->best_scores, 
          st->best_location, MAST_RANGE, pv, opts->sonly, opts->lump, norder);
#ifdef DEBUG
      if (isnan(st->Pvalue)) {
        fprintf(stderr, "nan: i %d d-ind %d s-ind %d strand %d pvalue %f\n", i, sseq->db_index, 
            sseq->index, st->strand, st->Pvalue); 
        for (j=0; j<alen; j++) fprintf(stderr, "%c %f\n", unhash(j), get_array_item(j, sseq->comp));
        abort();
      }
#endif
    }
    /* 
      print progress report
    */
    if (opts->status && i && i % 100 == 0) 
      fprintf(stderr, "\rrecalc p-value sequences: %6d ", i);

    /* 
      free pv space if E-value too big, otherwise
      save composition-based pv distribution
    */
    if ((sseq->pos_strand == NULL || (sseq->pos_strand->Pvalue * (*nseqs) > opts->e_max)) &&
        (sseq->neg_strand == NULL || (sseq->neg_strand->Pvalue * (*nseqs) > opts->e_max))) {
      free_2array(pv, nmotifs);
    } else {
      sseq->pv = pv;
      sseq->pv_alloc_len = nmotifs;
    }
  } // recalculate p-values 

} // get_scores 


/**********************************************************************/
/*
    Process the motif & spacing diagram.

    NB. The motif spacing diagram is stored in the globals:
    'diagram', 'dptr', 'norder', 'order' and 'space'.
*/
/**********************************************************************/
static void process_diagram(int file_count, int used_count, LO **los) {
  RBTREE_T *los_by_idx;
  RBTREE_T *used_idxs;
  LO *lo;
  int i, j, m;
  // NB: diagram, dptr, norder, order AND space ARE GLOBALS DEFINED IN diagram.h

  // First check that we have a diagram to parse
  if (diagram == NULL) {
    norder = 0; // no diagram
    return;
  }

  // create a lookup for the logodds matrices
  los_by_idx = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);
  for (i = 0; i < used_count; i++) {
    rbtree_make(los_by_idx, &(los[i]->imotif), los[i]);
  }

  // keep track of the used indexes in the diagram so a motif can't be repeated
  used_idxs = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);

  // Parse the ordering and spacing diagram.
  dptr = 0; // diag read pointer 
  norder = 0; // lnumber of motifs in diag 
  for (i=0; i<MAXLO; i++) space[i] = 0; // set all spacers to zero  
  // parse diagram; sets globals norder, order, space */
  if (yyparse()) exit(1);                    
  // check that no unused motifs are in motif diagram 
  for (i = 0; i < norder; i++) {
    m = order[i]; // motif number
    if (m < 1 || m > file_count) {
      error("Unknown motif %d given in motif diagram.\n", m);
    } else if (rbtree_find(los_by_idx, &m) == NULL) {
      error("Filtered motif %d given in motif diagram.\n", m);
    }
    if (!rbtree_make(used_idxs, &m, NULL)) {
      error("Motif %d used more than once in motif diagram:\n  %s\n", m, diagram);
    }
  }
  // quit if there was a problem in the diagram 
  exit_on_error();

  // change format for spacers to be relative to start of prior motif 
  space[0] = -1; // ignore first spacer 
  for (i = 1; i < norder; i++) {
    m = order[i-1]; // previous motif in diagram  
    lo = (LO*)rbtree_get(los_by_idx, &m);
    space[i] += lo->w; // internal format is -1 
  }

  // convert spacing diagram to motif positions in los array
  for (i = 0; i < norder; i++) {
    for (j = 0; j < used_count; j++) {
      if (order[i] == los[j]->imotif) break; // found motif by it's loading number
    }
    order[i] = j; // position for motif in los array
  }

  // cleanup
  rbtree_destroy(los_by_idx);
  rbtree_destroy(used_idxs);
}


/**********************************************************************/
/*
    Read in the MEME formatted motifs from mo_source which may be standard
    input. Filter the motifs.
*/
/**********************************************************************/
static ARRAYLST_T* load_motifs(
  MAST_OPTIONS_T *opts,
  XLATE_T *xlate, // sequence translator
  ALPH_T **alph_p, // the motif alphabet (OUT)
  ARRAY_T **bg_p, // the motif background (OUT)
  int *nmotifs_p // total number of motifs in the file before filtering (OUT)
) {
  int i, row, col;
  int nmotifs;
  ARRAYLST_T *motifs; 
  MOTIF_T *motif;
  MATRIX_T *scores;
  MREAD_T *mread;
  ALPH_T *alph;
  MTYPE min, max, value;
  ARRAY_T *motif_bg;

  mread = mread_create(opts->mo_source, OPEN_MFILE);
  if (opts->bfile) mread_set_bg_source(mread, opts->bfile); 	// TLB added 2-Feb-2017
  alph = alph_hold(mread_get_alphabet(mread));
  if (xlate != NULL && !alph_equal(alph, xlate_dest_alph(xlate))) {
    die("Motif alphabet does not match the output alphabet of the sequence translator.\n");
  }
  motif_bg = mread_get_background(mread);
  nmotifs = 0;
  motifs = arraylst_create();
  while ((motif = mread_next_motif(mread)) != NULL) {
    nmotifs++;
    // filter the motifs
    if (opts->first_motifs != 0 && get_motif_idx(motif) > opts->first_motifs) goto skip_motif;
    if (rbtree_size(opts->ids) || rbtree_size(opts->idxs)) {
      i = get_motif_idx(motif);
      if (!rbtree_find(opts->ids, get_motif_id(motif)) && !rbtree_find(opts->idxs, &i)) {
        goto skip_motif;
      }
    } else if (get_motif_evalue(motif) >= opts->min_ev) {
      goto skip_motif;
    }
    // check that the motif can be scaled
    scores = get_motif_scores(motif);
    min = max = get_matrix_cell(0, 0, scores);
    for (row = 0, col = 1; row < get_num_rows(scores); row++, col = 0) {
      for (; col < get_num_cols(scores); col++) {
        value = get_matrix_cell(row, col, scores);
        if (value < min) min = value;
        if (value > max) max = value;
      }
    }
    if (min == max) goto skip_motif; // motif can't be scaled!
    // shuffle the motif if requested (EXPERIMENTAL OPTION)
    if (opts->shuffle) shuffle_motif(motif, &(opts->prng));
    // store the motif
    arraylst_add(motif, motifs);
    continue;
skip_motif:
    destroy_motif(motif);
  }
  mread_destroy(mread);

  if (nmotifs == 0) die("No motifs found in file \"%s\".\n", opts->mo_source);
  if (arraylst_size(motifs) == 0) {
    die("The motif selection criteria is too strict and all motifs were excluded.\n");
  }

  // return results
  if (alph_p != NULL) {
    *alph_p = alph;
  } else {
    // done with the alphabet
    alph_release(alph);
  }
  if (bg_p != NULL) {
    *bg_p = motif_bg;
  } else {
    free_array(motif_bg);
  }
  if (nmotifs_p != NULL) *nmotifs_p = nmotifs;
  return motifs;
}

/**********************************************************************/
/*
    Read in the background.
*/
/**********************************************************************/
static ARRAY_T* load_background(
  ALPH_T *alph,
  bool translate,
  STYPE stype,
  ARRAY_T *motif_bg,
  char **bg_source // IN/OUT background model source; may get reset if non-standard alphabet
) {
  int i, bgorder;
  double freq;
  ARRAY_T *bgfreqs;
  // initialize the background frequencies and alphabets 
  bgorder = -1;
  if (*bg_source == NULL || strcmp(*bg_source, "--nrdb--")==0) {
    // try to maintain the old behaviour
    // Note that these will not sum to 1 because there are small
    // frequencies assigned to ambiguous characters.
    // This is how MAST used to work so it remains for compatibility
    bgfreqs = get_mast_frequencies(alph, false, translate);
    //if (*bg_source) myfree(*bg_source);		// Get set to set to keyword; LEAK: this fails sometimes???
    if (bgfreqs) {	// Alphabet was DNA or protein
      *bg_source = strdup("--nrdb--");
    } else {		// Alphabet was not standard; use motif frequencies
      assert(motif_bg != NULL);
      bgfreqs = allocate_array(get_array_length(motif_bg));
      copy_array(motif_bg, bgfreqs);
      // Set the bg_source to "--motif--";
      *bg_source = strdup("--motif--");
    }
    bgorder = 0;
  } else if (strcmp(*bg_source, "motif-file")==0 || strcmp(*bg_source, "--motif--")==0) {
    assert(motif_bg != NULL);
    bgfreqs = allocate_array(get_array_length(motif_bg));
    copy_array(motif_bg, bgfreqs);
    bgorder = 0;
  } else if (strcmp(*bg_source, "--uniform--")==0) {
    bgfreqs = allocate_array(alph_size_core(alph));
    init_array(1.0/alph_size_core(alph), bgfreqs);
    bgorder = 0;
  } else {
    bgfreqs = load_markov_model(alph, &bgorder, *bg_source);
  }
  // round (just because it used to, I don't understand the benefit)
  for (i = 0; i < get_array_length(bgfreqs); i++) {
    freq = get_array_item(i, bgfreqs);
    RND(freq, 8, freq);
    set_array_item(i, freq, bgfreqs);
  }
  // average reverse complement probabilities together
  if (stype == Separate || stype == Combine) {
    average_rc_markov_model(alph, bgorder, bgfreqs);
  }
  // now for some reason I don't yet understand MAST wants values for ambig symbols 
  // but it's quite happy if they're zero
  resize_markov_model(alph_size_core(alph), alph_size_full(alph), bgfreqs, NULL);

  // return results
  return bgfreqs;
}

/************************************************************************/
/*
    Mark motifs that should be removed due to excessive correlation.

    Returns number of marked motifs.
*/
/************************************************************************/
static int mark_correlated_motifs(
  double maxcorr, // maximum allowed corellation 
  ARRAYLST_T *all_motifs, // original motif list
  LO **los, // array of logodds matrices 
  MATRIX_T *corr // matrix of correlations
)
{
  int i, j, ncorrelated, nmotifs;
  nmotifs = arraylst_size(all_motifs);
  ncorrelated = 0;
  for (i = 1; i < nmotifs; i++) { // from motif 
    for (j = 0; j < i; j++) { // to motif 
      if (!los[j]->is_bad && get_matrix_cell(i, j, corr) > maxcorr) {
        los[i]->is_bad = true;
        set_motif_mark((MOTIF_T*)arraylst_get(i, all_motifs), BAD_MOTIF_MARK);
        ncorrelated++;
        break; //no point in comparing this motif further
      }
    } // to motif 
  } // from motif 

  return ncorrelated;
}

/************************************************************************/
/*
    Create and print a list of non-overlapping hits for a single database.
*/
/************************************************************************/
static void print_db_hit_list(
  FILE *mast_out,
  MAST_OPTIONS_T *opts,
  ALPH_T *alph,
  XLATE_T *xlate,
  LO *los[], // array of pointers to lo matrices 
  int nmotifs, // number motifs read 
  double **pv, // p-value tables for each motif 
  DATABASE_T *db // source database 
) 
{
  long length; // length of sample 
  char *sample_name; // name of sample 
  char *sequence; // sequence of sample 
  char *id; // identifier text for sample 
  TILING tiling; // tiling and diagram of sequence 
  int nseqs = 0;

  if (opts->best_motifs) {
    fprintf(mast_out,"# Best single (non-overlapping) hit for each motif in all sequences from \"%s\".\n", db->source);
  } else {
    fprintf(mast_out,"# All non-overlapping hits in all sequences from \"%s\".\n", db->source);
  }
  fprintf(mast_out,"# sequence_name (strand+/-)motif id alt_id hit_start hit_end score hit_p-value\n");

  while (read_sequence(alph, db->file, &sample_name, &id, &sequence, &length)) {
    // Skip over special MEME WEIGHTS extension
    if (length == 0 && strcmp(sample_name, "WEIGHTS") == 0) {
      myfree(sequence);
      myfree(id);
      myfree(sample_name);
      continue;
    }
    nseqs++;
    if (opts->status && nseqs % 100 == 0) fprintf(stderr, "\rsequences: %6d ", nseqs);
    tiling = score_tile_diagram(mast_out, alph, xlate, sequence, length, los, nmotifs, opts->stype, false,
      opts->best_motifs, false, pv, opts->m_thresh, opts->w_thresh, false, true, sample_name);
    if (opts->stype == Separate) {
      tiling = score_tile_diagram(mast_out, alph, xlate, sequence, length, los, nmotifs, opts->stype, true,
        opts->best_motifs, false, pv, opts->m_thresh, opts->w_thresh, false, true, sample_name);
    }
    myfree(sample_name);
    myfree(sequence);
    myfree(id);
  } // read_sequence
  if (opts->status) fprintf(stderr, "\n");
} // print_db_hit_list

/************************************************************************/
/*
    Create and print a list of non-overlapping hits for all databases.
*/
/************************************************************************/
static void print_hit_list(
  int argc,
  char **argv,
  MAST_OPTIONS_T *opts,
  ALPH_T *alph,
  XLATE_T *xlate,
  LO *los[], // array of pointers to lo matrices 
  int nmotifs, // number motifs read 
  double **pv, // p-value tables for each motif 
  ARRAYLST_T *dbs
) {
  int i;
  // Create and print a hit list to stdout and exit
  for (i = 0; i < arraylst_size(dbs); ++i) {
    DATABASE_T *db = (DATABASE_T*)arraylst_get(i, dbs);
    print_db_hit_list(stdout, opts, alph, xlate, los, nmotifs, pv, db);
  }
  // output arguments 
  fprintf(stdout, "# mast");
  for (i = 1; i < argc; ++i) fprintf(stdout, " %s", argv[i]);
  fprintf(stdout, "\n");
}


/**********************************************************************/
/*
   to cater for multiple databases, mast can read a tab delimited file with
   one database record on each line.
  
   Returns array list of databases
 */
/**********************************************************************/
static ARRAYLST_T* init_databases(MAST_OPTIONS_T *opts, ALPH_T *alph) {
  int i;
  ARRAYLST_T *databases;
  DATABASE_T *db;
  time_t now;
  // check if multiple databases have been specified
  if (opts->db_list) {
    //read the tab delimited file
    FILE *fp;
    BUF_T *buffer;
    int col, delim;
    char *src, *name, *link, *token;
    databases = arraylst_create();
    if (strcmp(opts->db_source, "-") == 0) {
      fp = stdin;
    } else if (file_exists(opts->db_source)) {
      fp = fopen(opts->db_source, "r");
      if (!fp) die("Error reading file specified for database list \"%s\","
          " error was given as %s.\n", opts->db_source, strerror(errno));
    } else {
      fp = NULL;//stop warnings about "uninitilized variables"...
      die("File specified for database list \"%s\" doesn't exist!\n", opts->db_source);
    }
    buffer = buf_create(100);
    buf_flip(buffer);
    col = 0, src = NULL, name = NULL, link = NULL;
    //read the file token by token, only keep the first 3 tokens in a line
    while (!feof(fp) || buf_remaining(buffer)) {
      token = buf_fread_token(buffer, fp, buf_is_delim_letter, "\t\n\r", 
          false, NULL, 0, NULL);
      delim = buf_getc(buffer);
      switch (col++) {
        case 0: src = token; break;
        case 1: name = token; break;
        case 2: link = token; break;
        default: myfree(token);
      }
      if (delim == -1 || delim == '\n' || delim == '\r') {
        buf_fread_consume(buffer, fp, buf_is_delim_letter, "\n\r", false);
        if (!file_exists(src)) {
          error("File list contains a file \"%s\" "
              "that doesn't exist!\n", src);
        }
        db = mm_malloc(sizeof(DATABASE_T));
        memset(db, 0, sizeof(DATABASE_T));
        db->source = src;
        db->name = name;
        db->link = link;
        arraylst_add(db, databases);
        src = NULL, name = NULL, link = NULL;
        col = 0;
      }
    }
    if (src != NULL) myfree(src);
    if (name != NULL) myfree(name);
    if (link != NULL) myfree(link);
    arraylst_fit(databases);
    buf_destroy(buffer);
    fclose(fp);
    exit_on_error();
  } else {
    if (strcmp(opts->db_source, "-") != 0) {
      if (!file_exists(opts->db_source)) die("File specified for database \"%s\""
          " doesn't exist!\n", opts->db_source);
    }
    db = mm_malloc(sizeof(DATABASE_T));
    memset(db, 0, sizeof(DATABASE_T));
    copy_string(&(db->source), opts->db_source);
    copy_string(&(db->name), opts->db_name);
    copy_string(&(db->link), opts->db_link);
    databases = arraylst_create_sized(1);
    arraylst_add(db, databases);
  }
  now = time(NULL);//for the default last mod date of databases
  for (i = 0; i < arraylst_size(databases); i++) {
    db = arraylst_get(i, databases);
    db->index = i;
    db->alph = alph_hold(alph);
    db->last_mod = now;
    db->sequence_count = 0;
    db->residue_count = 0;
    if (strcmp(db->source, "-") == 0) {
      db->file = stdin;
      db->save = tmpfile();
      if (db->save == NULL) die("Error creating temporary file for "
          "storing stdin, error given as: %s\n", strerror(errno));
    } else {
      db->file = fopen(db->source, "r");
      db->save = NULL;
      if (db->file == NULL) die("Error opening sequence file \"%s\", "
          "error given as: %s\n", db->source, strerror(errno));
#ifdef UNIX
      {//update the last mod date to the correct value
        struct stat stbuf;
        stat(db->source, &stbuf);
        db->last_mod = stbuf.st_mtime;
      }
#endif
    }
  }
  return databases;
}

static void destroy_database(void *p) {
  DATABASE_T *db;
  db = (DATABASE_T*)p;
  alph_release(db->alph);
  myfree(db->source);
  myfree(db->name);
  myfree(db->link);
  if (db->save) fclose(db->save);
  if (db->file != stdin) fclose(db->file);
  memset(db, 0, sizeof(DATABASE_T));
  myfree(db);
}

/**********************************************************************/
/* 
 * score_asc
 *
 * Compares two strands to allow sorting in the array list.
 * Sorts primarily by combined p-value ascending and then breaks
 * ties by sorting on database, file position in database and then
 * on strand. This should mean that no two strands are considered equal.
 */
/**********************************************************************/
static int score_asc(const void *v1, const void *v2) {
  STRAND_T *hit1, *hit2;
  double diff;
  int db_diff;
  long fp_diff;
  hit1 = *((STRAND_T**)v1);
  hit2 = *((STRAND_T**)v2);
  //order by score
  diff = hit1->Pvalue - hit2->Pvalue;
  if (diff < 0) {
    return -1;
  } else if (diff > 0) {
    return 1;
  }
  //order by database
  db_diff = hit1->sequence->db_index - hit2->sequence->db_index;
  if (db_diff < 0) {
    return -1;
  } else if (db_diff > 0) {
    return 1;
  }
  //order by file position
  fp_diff = hit1->sequence->fp - hit2->sequence->fp;
  if (fp_diff < 0) {
    return -1;
  } else if (fp_diff > 0) {
    return 1;
  }
  //must be the same sequence
  //order by strand
  if (hit1->strand > hit2->strand) {
    return -1;
  } else if (hit1->strand < hit2->strand) {
    return 1;
  }
  //I don't think this is possible
  return 0;
}

/**********************************************************************/
/*
 * check_chunk
 *
 * checks a part of the hits array of size chunk starting at start
 * for motif hits. Checks if the hits found would overlap the end
 * of the chunk and makes a prediction on how many more chunks would
 * need to be checked.
 *
 * Returns 0 if there is nothing in this chunk, returns 1 if there is
 * something in this chunk but it does not overlap, returns n if there
 * is something in this chunk and at minimum it would overlap a further
 * n-1 chunks. Assumes that hits can not overlap the end of the sequence.
 */
/**********************************************************************/
static int check_chunk(int length, int *hits, LO **motifs, int start, int chunk) {
  int i, end, found;
  found = 0;
  end = MIN(length, start + chunk);
  for (i = start; i < end; ++i) {
    if (hits[i]) {
      found = 1;
      i += motifs[abs(hits[i]) - 1]->ws;
      if (i > end) {
        found = ((int)((i - end) / chunk)) + 2;
      }
      --i;//because of the loop inc
    }
  }
  return found;
}

/**********************************************************************/
/*
 * create_hit_match_patterns
 *
 * Creates the patterns showing which characters mast believes are good
 * matches and the translation from DNA to Protein to compare with the 
 * motif (only if translating dna).
 */
/**********************************************************************/
static void create_hit_match_patterns(ALPH_T *alph, XLATE_T *xlate, char *sequence,
    int pos, bool neg_strand, LO *motif,
    char *matches, char *translation, int buflen) {
  char *s;
  int i, inc, w, ws;
  double scale, offset;
  double **logodds;
  if (xlate != NULL) alph = xlate_dest_alph(xlate);
  // double check that we can do what we want
  if (buflen <= motif->w) die("Buffer is too short to generate motif match patterns.\n");
  // the width and width in sequence of the motif
  w = motif->w;
  ws = motif->ws; //should be w * inc
  // the motif log odds matrix
  logodds = motif->logodds;
  // the details required for converting log odds scores
  scale = motif->scale;
  offset = motif->offset;
  // increment in sequence for each entry in motif
  inc = (xlate != NULL ? xlate_src_nsyms(xlate) : 1);
  // for each entry in the motif, calculate if the sequence is a good match
  // and if the underlying sequence is a different alphabet then
  // translate the sequence into the motif alphabet
  for (i = 0, s = sequence+pos; i < w; i++, s += inc) {
    int cx;
    double score;
    if (xlate != NULL) {
      cx = xlate_index(xlate, neg_strand, s);
    } else {
      if (neg_strand) {
        cx = alph_complement(alph, alph_index(alph, *s));
      } else {
        cx = alph_index(alph, *s);
      }
    }
    score = neg_strand ? logodds[w-i-1][cx] : logodds[i][cx]; 
    matches[i] = (scaled_to_bit(score, 1, scale, offset) > 0) ? '+' : ' ';
    if (xlate != NULL) translation[i] = alph_char(alph, cx);
  }
  matches[w] = '\0';
  if (xlate != NULL) translation[w] = '\0';
}

/**********************************************************************/
/*
 * output_mast_xml_hit
 *
 * Outputs a hit for the current position if one exists.
 */
/**********************************************************************/
static inline void output_mast_xml_hit(
    FILE *out, ALPH_T *alph, XLATE_T *xlate, char *sequence, int pos,
    TILING tiling, LO **motifs, RBTREE_T *idx2pos,
    char *matches_buffer, char *translation_buffer, 
    int buffer_len, STYPE stype) {
  int hit, gap;
  LO *motif;
  bool neg_strand;
  double pvalue;
  hit = tiling.hits[pos];
  //check that a hit exists at that location
  if (!hit) return; //no motif, return
  motif = motifs[abs(hit) - 1];
  neg_strand = (hit < 0);
  pvalue = tiling.pvalues[pos];
  create_hit_match_patterns(alph, xlate, sequence, pos, neg_strand, motif,
      matches_buffer, translation_buffer, buffer_len);
  fprintf(out, "\t\t\t\t<hit pos=\"%d\" idx=\"%d\"", pos + 1, 
      *((int*)rbtree_get(idx2pos, &(motif->imotif))));
  if (stype != Unstranded && stype != Norc) {
    fprintf(out, " rc=\"%s\"", (neg_strand ? "y" : "n"));
  }
  fprintf(out, " pvalue=\"%.1e\" match=\"%s\"", pvalue, matches_buffer);
  if (xlate != NULL) fprintf(out, " translation=\"%s\"", translation_buffer);
  fprintf(out, "/>\n");
}

/**********************************************************************/
/*
 * Outputs a motif in XML
 */
/**********************************************************************/
static void print_xml_motif(FILE *file, int db, MOTIF_T *motif, char *indent, char *tab) {
  int i, j, len;
  char *id, *alt, *url;
  STR_T *b;
  ALPH_T *alph;
  MATRIX_T *scores;
  double A, C, G, T, log_ev;

  b = str_create(10);
  id = get_motif_id(motif);
  alph = get_motif_alph(motif);
  alt = get_motif_id2(motif);
  len = get_motif_length(motif);
  url = get_motif_url(motif);
  scores = get_motif_scores(motif);
  log_ev = get_motif_log_evalue(motif);

  fprintf(file, "%s<motif db=\"%d\" ", indent, db);
  fprintf(file, "id=\"%s\" ", xmlify(id, b, true));
  if (alt && alt[0] != '\0') fprintf(file, "alt=\"%s\" ", xmlify(alt, b, true));
  fprintf(file, "length=\"%d\"", len);
  if (get_motif_nsites(motif) > 0)
    fprintf(file, " nsites=\"%g\"", get_motif_nsites(motif));
  if (log_ev > -HUGE_VAL && log_ev < HUGE_VAL)
    fprintf(file, " evalue=\"%s\"", str_evalue(b, log_ev, 1));
  if (test_motif_mark(motif, BAD_MOTIF_MARK)) {
    fprintf(file, " bad=\"y\"");
  }
  if (url && url[0] != '\0') {
    fprintf(file, "\n%s%s%surl=\"%s\"", indent, tab, tab, xmlify(url, b, true));
  }
  fprintf(file, ">\n");
  for (i = 0; i < len; ++i) {
    fprintf(file, "%s%s<pos", indent, tab);
    for (j = 0; j < alph_size_core(alph); j++) {
      fprintf(file, " %s=\"%g\"", alph_xml_id(alph, j, b),
          get_matrix_cell(i, j, scores));
    }
    fprintf(file, "/>\n");
  }
  fprintf(file, "%s</motif>\n", indent);
  str_destroy(b, false);
}

/**********************************************************************/
/*
 * output_mast_xml
 *
 * Outputs mast xml results.
 *
 * returns the file name, caller is responsible for freeing file name.
 */
/**********************************************************************/
void output_mast_xml(
    MAST_OPTIONS_T *opts, ALPH_T *alph, XLATE_T *xlate, ARRAY_T *bg, ARRAY_T *motif_bg,
    int nmotifs, int nbad, LO **good_motifs, ARRAYLST_T *all_motifs, MATRIX_T *corr,
    time_t motifs_last_mod,
    int nos_norder, int *nos_order, int *nos_space,
    int nseqs, ARRAYLST_T *databases, ARRAYLST_T *sequence_scores, 
    int argc, char **argv) {
  int i, j, seg_start, seg_width, lpos, longest_motif_len;
  char *xml_path, *strand_handling_str, *bg_source_str, *matches, *translation;
  FILE *out;
  STR_T *b;
  double cycles;
  RBTREE_T *idx2pos;
  // open file for output
  if (opts->out_dir != NULL) {
    // configure output
    if (create_output_directory(opts->out_dir, opts->clobber, opts->status)) {
      // Failed to create output directory.
      exit(1);
    }
    xml_path = make_path_to_file(opts->out_dir, XML_FILENAME);
    out = fopen(xml_path, "w");
  } else {
    xml_path = NULL;
    out = stdout;
  }
  // map motif loading index to the output ordering
  idx2pos = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, rbtree_intcpy, free);
  for (i = 0; i < arraylst_size(all_motifs); i++) {
    int idx = get_motif_idx((MOTIF_T*)arraylst_get(i, all_motifs));
    rbtree_put(idx2pos, &idx, &i);
  }
  // create xmlify buffer
  b = str_create(10);

  // MAST
  fprintf(out, "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n");
  fprintf(out, "<mast version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n");
  
  // COMMAND LINE
  fprintf(out, "\t<command_line><arg>mast</arg>");
  for (i = 1; i < argc; ++i) fprintf(out, "<arg>%s</arg>", xmlify(argv[i], b, false));
  fprintf(out, "</command_line>\n");

  // MISCELLANEOUS SETTINGS
  switch (opts->stype) {
    case Combine:     strand_handling_str = "combine";    break;
    case Separate:    strand_handling_str = "separate";   break;
    case Norc:        strand_handling_str = "norc";       break;
    case Unstranded:  strand_handling_str = "unstranded"; break;
    default:          
      strand_handling_str = NULL;
      die("Impossible value for strand_handling.\n");
  }
  fprintf(out, "\t<settings"
      " strand_handling=\"%s\""
      " max_correlation=\"%.2f\""
      " remove_correlated=\"%c\""
      " max_seq_evalue=\"%.2g\""
      " adjust_hit=\"%c\""
      " max_hit_pvalue=\"%g\""
      " max_weak_pvalue=\"%g\""
      "/>\n",
      strand_handling_str,
      MAXCORR,
      (opts->rem_corr ? 'y' : 'n'),
      opts->e_max,
      (opts->adj_hit_pvalue ? 'y' : 'n'),
      opts->m_thresh,
      opts->w_thresh);

  // ALPHABETS
  // main alphabet
  alph_print_xml(alph, "alphabet", "\t", "\t", out);
  if (xlate != NULL) {
    // sequence alphabet (if different)
    alph_print_xml(xlate_src_alph(xlate), "sequence_alphabet", "\t", "\t", out);
    // translations symbols in to symbols out
    fprintf(out, "\t<translate num_sequence=\"%d\" num_motif=\"%d\"/>\n", xlate_src_nsyms(xlate), 1);
  }
  // BACKGROUND
  //get the alphabet background source
  bg_source_str = (opts->use_seq_comp ? "--sequence--" : opts->bfile);
  fprintf(out, "\t<background source=\"%s\"", bg_source_str);
  for (i = 0; i < alph_size_core(alph); i++) {
    fprintf(out, " %s=\"%g\"", alph_xml_id(alph, i, b), get_array_item(i, bg));
  }
  fprintf(out, "/>\n");
  //if (opts->bfile != NULL) {
    //fprintf(out, " file=\"%s\"/>\n", xmlify(opts->bfile, b, true));
  //} else {
  //  fprintf(out, "/>\n");
  //}

  // MOTIF DBs
  fprintf(out, "\t<motif_dbs>\n");
  fprintf(out, "\t\t<motif_db source=\"%s\" ", xmlify(opts->mo_source, b, true));
  if (opts->mo_name) fprintf(out, "name=\"%s\" ",  xmlify(opts->mo_name, b, true));
  fprintf(out, "last_mod_date=\"%s\">\n", strtok(ctime(&motifs_last_mod), "\n"));
  fprintf(out, "\t\t\t<background");
  for (i = 0; i < alph_size_core(alph); i++) {
    fprintf(out, " %s=\"%g\"", alph_xml_id(alph, i, b), get_array_item(i, motif_bg));
  }
  fprintf(out, "/>\n");
  fprintf(out, "\t\t</motif_db>\n");
  fprintf(out, "\t</motif_dbs>\n");

  // MOTIFS
  fprintf(out, "\t<motifs>\n");
  for (i = 0; i < arraylst_size(all_motifs); i++) {
    print_xml_motif(out, 0, (MOTIF_T*)arraylst_get(i, all_motifs), "\t\t", "\t");
  }
  fprintf(out, "\t</motifs>\n");

  // MOTIF CORRELATIONS
  fprintf(out, "\t<correlations>\n");
  for (i = 1; i < arraylst_size(all_motifs); ++i) {
    for (j = 0; j < i; ++j) {
      fprintf(out, "\t\t<correlation idx_a=\"%d\" idx_b=\"%d\""
          " value=\"%.2f\"/>\n", j, i, get_matrix_cell(i, j, corr));
    }
  }
  fprintf(out, "\t</correlations>\n");

  // NOMINAL ORDER AND SPACING OF MOTIFS
  if (nos_norder) {
    int gap, pos, prev_len;
    // It appears from reading process_diagram that the spacer before the first 
    // motif is ignored and the spacer after the last motif isn't loaded.
    // Additionally it appears that the length of the gaps between the motifs 
    // is assumed to be in the same alphabet as the motifs.
    fprintf(out, "\t<nos>\n");
    for (i = 0, prev_len = 0, gap = 0; i < nos_norder; ++i) {
      // note that nos_order is the position in the good_motifs array,
      // we need to convert it to the position in all_motifs
      pos = *((int*)rbtree_get(idx2pos, &(good_motifs[nos_order[i]]->imotif)));
      if (i == 0) {
        fprintf(out, "\t\t<expect idx=\"%d\"/>\n", pos);
      } else {
        // length is not a bug, this diagram is not scaled to the units of the sequence
        gap = nos_space[i] - prev_len;
        fprintf(out, "\t\t<expect gap=\"%d\" idx=\"%d\"/>\n", gap, pos);
      }
      prev_len = get_motif_length((MOTIF_T*)arraylst_get(pos, all_motifs));
    }
    fprintf(out, "\t</nos>\n");
  }

  // SEQUENCE DBs
  fprintf(out, "\t<sequence_dbs>\n");
  // output each database
  for (i = 0; i < arraylst_size(databases); ++i) {
    DATABASE_T *db;
    db = (DATABASE_T*)arraylst_get(i, databases);
    fprintf(out, "\t\t<sequence_db source=\"%s\"", xmlify(db->source, b, true));
    if (db->name) fprintf(out, " name=\"%s\"", xmlify(db->name, b, true));
    //note that ctime puts a newline on the end of the time string, so I use strtok to null it
    fprintf(out, " last_mod_date=\"%s\" seq_count=\"%d\" residue_count=\"%ld\"",
        strtok(ctime(&(db->last_mod)), "\n"), db->sequence_count, db->residue_count);
    if (db->link) fprintf(out, " link=\"%s\"", xmlify(db->link, b, true));
    fprintf(out, "/>\n");
  }
  fprintf(out, "\t</sequence_dbs>\n");

  // SEQUENCES
  fprintf(out, "\t<sequences>\n");
  //calculate longest motif used in scoring so I know what size buffers to allocate
  longest_motif_len = 0;
  for (i = 0; i < nmotifs; ++i) {
    LO *motif = good_motifs[i];
    if (motif->w > longest_motif_len) longest_motif_len = motif->w;
  }
  matches = mymalloc(sizeof(char) * (longest_motif_len + 1));
  translation = mymalloc(sizeof(char) * (longest_motif_len + 1));
  //output each high scoring sequence
  for (i = 0; i < arraylst_size(sequence_scores); i++) {
    FILE *seq_file;
    char *sample_name, *sample_comment, *sample_sequence;
    long sample_length;
    TILING score_tiling, other_tiling;
    SSEQ_T *sequence;
    STRAND_T *score, *other;
    DATABASE_T *db;
    bool neg_strand;

    //get the next score and its sequence.
    sequence = ((score = (STRAND_T*)arraylst_get(i, sequence_scores))->sequence);

    //check that this sequence hasn't already been printed
    if (opts->stype == Separate) {
      neg_strand = (score->strand == -1);
      // if there was the possibility of two strands then get both
      // and check that the sequence hasn't been printed already.
      if (!neg_strand) other = sequence->neg_strand;
      else other = sequence->pos_strand;
      if (other) {
        // check if this sequence has already been printed
        if (other->Pvalue < score->Pvalue || (other->Pvalue == score->Pvalue && score->strand == -1)) 
          continue;
        // check that the evalue is acceptable
        if (other->Pvalue * nseqs > opts->e_max)
          other = NULL;
      }
    } else {
      neg_strand = false;
      other = NULL;
    }

    //get the database associated with the sequence
    db = (DATABASE_T*)arraylst_get(sequence->db_index, databases);
    //get the file that the sequence should be loaded from
    seq_file = (db->save ? db->save : db->file);

    // reload the sequence from file
    fseek(seq_file, sequence->fp, SEEK_SET);
    read_sequence(alph, seq_file, &sample_name, &sample_comment, &sample_sequence, &sample_length);
    // we should have already skipped the MEME WEIGHTS extension
    if (sample_length == 0 && strcmp(sample_name, "WEIGHTS") == 0) {
      error("Found MEME WEIGHTS section that we should have already skipped!");
      exit_on_error();
    }

    fprintf(out, "\t\t<sequence db=\"%d\" name=\"%s\" ",
        sequence->db_index, xmlify(sample_name, b, true));
    fprintf(out, "comment=\"%s\" length=\"%ld\">\n",
        xmlify(sample_comment, b, true), sample_length);

    // score, tile and diagram the sequence with each of the motifs 
    score_tiling = score_tile_diagram(NULL, alph, xlate, sample_sequence,
        sample_length, good_motifs, nmotifs, opts->stype, neg_strand, 
        false, false, sequence->pv, opts->m_thresh, opts->w_thresh, 
        opts->adj_hit_pvalue, false, sample_name);

    fprintf(out, "\t\t\t<score strand=\"%s\" combined_pvalue=\"%.2e\" evalue=\"%.2g\"", 
        (score->strand == 0 ? "both" : 
         (score->strand == 1 ? "forward" : "reverse")), 
        score->Pvalue, score->Pvalue * nseqs);
    if (xlate != NULL) fprintf(out, " frame=\"%c\"", best_frame(xlate, nmotifs, sample_length, score_tiling));
    fprintf(out, "/>\n");
    if (opts->stype == Separate) {
      // score, tile and diagram the sequence with each of the motifs 
      other_tiling = score_tile_diagram(NULL, alph, xlate, sample_sequence,
          sample_length, good_motifs, nmotifs, opts->stype, !neg_strand,
          false, false, sequence->pv, opts->m_thresh, opts->w_thresh,
          opts->adj_hit_pvalue, false, sample_name);
    } else { // avoid compilier complaints
      memset(&other_tiling, 0, sizeof(TILING));
    }

    if (other) {
      fprintf(out, "\t\t\t<score strand=\"%s\" combined_pvalue=\"%.2e\" evalue=\"%.2g\"", 
          (other->strand == 1 ? "forward" : "reverse"), other->Pvalue, other->Pvalue * nseqs);
      if (xlate != NULL) fprintf(out, " frame=\"%c\"", best_frame(xlate, nmotifs, sample_length, other_tiling));
      fprintf(out, "/>\n");
    } 
    //now output sequence segments in multiplies of SEG_CHUNK skipping segments that don't have hits
    seg_start = 0;
    while (seg_start < sample_length) {
      {
        int chunk, chunks, primary, secondary, start;
        for (chunk = 0, chunks = 0, start = seg_start; chunk == 0 || chunk < chunks; ++chunk, start += SEG_CHUNK) {
          primary = check_chunk(sample_length, score_tiling.hits, good_motifs, start, SEG_CHUNK);
          if (opts->stype == Separate) {
            secondary = check_chunk(sample_length, other_tiling.hits, good_motifs, start, SEG_CHUNK);
            chunks = MAX(chunks, chunk + MAX(primary, secondary));
          } else {
            chunks = MAX(chunks, chunk + primary);
          }
        }
        if (chunks) {
          seg_width = MIN(sample_length - seg_start, SEG_CHUNK * chunks);
        } else {
          seg_width = -SEG_CHUNK;
        }
      }
      if (seg_width > 0) {
        int offset, pos, seg_end;
        fprintf(out, "\t\t\t<seg start=\"%d\">\n", seg_start + 1);
        fprintf(out, "\t\t\t\t<data>\n");
        for (offset = 0; offset < seg_width; offset += SEG_CHUNK) {
          fprintf(out, SEG_FORMAT, sample_sequence+(seg_start+offset));
        }
        fprintf(out, "\t\t\t\t</data>\n");
        pos = seg_start;
        seg_end = seg_start + seg_width;
        for (pos = seg_start; pos < seg_end; ++pos) {
          output_mast_xml_hit(out, alph, xlate, sample_sequence, pos,
              score_tiling, good_motifs, idx2pos, matches, translation, 
              longest_motif_len + 1, opts->stype);
          if (opts->stype == Separate)
            output_mast_xml_hit(out, alph, xlate, sample_sequence, pos, 
                other_tiling, good_motifs, idx2pos, matches, translation, 
                longest_motif_len + 1, opts->stype);
        }
        fprintf(out, "\t\t\t</seg>\n");
      }
      seg_start += abs(seg_width);
    }

    //clean up TILING objects, which have allocations in them.
    myfree(score_tiling.hits);
    myfree(score_tiling.pvalues);
    myfree(score_tiling.svalues);
    myfree(score_tiling.diagram);
    if (opts->stype == Separate) {
      myfree(other_tiling.hits);
      myfree(other_tiling.pvalues);
      myfree(other_tiling.svalues);
      myfree(other_tiling.diagram);
    }
    //clean up
    myfree(sample_name);
    myfree(sample_comment);
    myfree(sample_sequence);
    fprintf(out, "\t\t</sequence>\n");
  }
  myfree(matches);
  myfree(translation);
  fprintf(out, "\t</sequences>\n");
  //the run time has to be printed last, or it wouldn't be the run time
  //as it has to be in the xml then it can't include the time to create the html
  cycles = myclock();
  fprintf(out, "\t<runtime host=\"%s\" when=\"%s\" cycles=\"%.0f\" seconds=\"%.3f\"/>\n",
      hostname(), strtok(ctime(&(opts->start_time)),"\n"), cycles, cycles / CLOCKS_PER_SEC);
  fprintf(out, "</mast>\n");
  rbtree_destroy(idx2pos);
  str_destroy(b, false);
  // finish xml output
  if (xml_path != NULL) {
    fclose(out);
    free(xml_path);
  }
}

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage = 
    "\n"
    "Usage: mast [options] <motif file> <sequence file>\n"
    "\n"
    "   Options:\n"
    "     -bfile <bf>    read background frequencies from <bf>\n"
    "     -dblist        the file specified as database contains a list of databases\n"
    "     -o <dir>       directory to output mast results\n"
    "     -oc <dir>      directory to output mast results with overwriting allowed\n"
    "     -hit_list      print only a list of non-overlapping hits to stdout\n"
    "     -remcorr       remove highly correlated motifs from query\n"
    "     -m <id>        use only motif(s) named <id> (overrides -mev);\n"
    "                      can be repeated\n"
    "     -mi <m>        use only motif(s) numbered <m> (overrides -mev);\n"
    "                      can be repeated\n"
    "     -c <count>     only use the first <count> motifs or all motifs\n"
    "                      when <count> is zero\n"
    "     -mev <mev>     use only motifs with E-values less than <mev>\n"
    "     -diag <diag>   nominal order and spacing of motifs\n"
    "     -norc          do not score reverse complement DNA strand\n"
    "     -sep           score reverse complement DNA strand as a separate sequence\n"
    "     -dna           translate DNA sequences to protein\n"
    "     -comp          adjust p-values and E-values for sequence composition\n"
    "     -ev <ev>       print results for sequences with E-value < <ev>;\n"
    "                      default: %g\n"
    "     -mt <mt>       show motif matches with p-value < mt; default: %g\n"
    "     -w             show weak matches (mt < p-value < mt*10) in angle brackets\n"
    "     -best          include only the best motif in diagrams;\n"
    "                      hit_list mode only\n"
    "     -seqp          use SEQUENCE p-values for motif thresholds\n"
    "                      default: use POSITION p-values\n"
    "     -mv <mf>       in results use <mf> as motif file name\n"
    "     -df <df>       in results use <df> as database name; ignored when\n"
    "                      option -dblist is specified\n"
    "     -dl <dl>       in results use <dl> as link to search sequence names;\n"
    "                      ignored when -dblist specified\n"
    "     -minseqs <ms>  lower bound on number of sequences in db\n"
    "     -nostatus      do not print progress report\n"
    "     -notext        do not generate text output\n"
    "     -nohtml        do not generate html output\n"
    "     -version       print the version and exit\n"
    ;
  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, usage, EXPECT, MOTIF_THRESH_DEFAULT);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts { OPT_BFILE, OPT_DBLIST, OPT_O, OPT_OC, OPT_HITLIST, OPT_REMCORE,
  OPT_M, OPT_MI, OPT_C, OPT_MEV, OPT_DIAG, OPT_NORC, OPT_SEP, OPT_DNA, OPT_COMP,
  OPT_EV, OPT_MT, OPT_W, OPT_BEST, OPT_SEQP, OPT_MF, OPT_DF, OPT_DL,
  OPT_MINSEQS, OPT_NOSTATUS, OPT_NOTEXT, OPT_NOHTML, OPT_VERSION, OPT_SHUFFLE,
  OPT_SONLY, OPT_LUMP };

/**************************************************************************
 * Checks the consistency of arguments and loads them into the options 
 * structure for easy access.
 **************************************************************************/
static void process_arguments(int argc, char **argv, MAST_OPTIONS_T *options) {
  // Note options are specified as follows:
  // <name> <has argument> <(not used)> <int to return>
  struct option mast_options[] = {
    {"bfile", required_argument, NULL, OPT_BFILE},
    {"dblist", no_argument, NULL, OPT_DBLIST},
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"hit_list", no_argument, NULL, OPT_HITLIST},
    {"remcorr", no_argument, NULL, OPT_REMCORE},
    {"m", required_argument, NULL, OPT_M},
    {"mi", required_argument, NULL, OPT_MI},
    {"c", required_argument, NULL, OPT_C},
    {"mev", required_argument, NULL, OPT_MEV},
    {"diag", required_argument, NULL, OPT_DIAG},
    {"norc", no_argument, NULL, OPT_NORC},
    {"sep", no_argument, NULL, OPT_SEP},
    {"dna", no_argument, NULL, OPT_DNA},
    {"comp", no_argument, NULL, OPT_COMP},
    {"ev", required_argument, NULL, OPT_EV},
    {"mt", required_argument, NULL, OPT_MT},
    {"w", no_argument, NULL, OPT_W},
    {"best", no_argument, NULL, OPT_BEST},
    {"seqp", no_argument, NULL, OPT_SEQP},
    {"mf", required_argument, NULL, OPT_MF},
    {"df", required_argument, NULL, OPT_DF},
    {"dl", required_argument, NULL, OPT_DL},
    {"minseqs", required_argument, NULL, OPT_MINSEQS},
    {"nostatus", no_argument, NULL, OPT_NOSTATUS},
    {"notext", no_argument, NULL, OPT_NOTEXT},
    {"nohtml", no_argument, NULL, OPT_NOHTML},
    {"version", no_argument, NULL, OPT_VERSION},
    // experimental options
    {"shuffle", no_argument, NULL, OPT_SHUFFLE}, // shuffle columns of motifs
    {"sonly", no_argument, NULL, OPT_SONLY}, // use only spacing p-values in product
    {"lump", no_argument, NULL, OPT_LUMP}, // combine spacings into one p-value
    {NULL, 0, NULL, 0} //boundary indicator
  };
  int idx;
  bool weak;
  weak = false;

  // set option defaults
  memset(options, 0, sizeof(MAST_OPTIONS_T));
  options->start_time = time(NULL);
  options->out_dir = default_out_dir;
  options->clobber = true;

  options->bfile = NULL;
  options->mo_source = "";
  options->db_source = "";
  options->db_list = false;
  options->diagram = NULL;
  options->mo_name = NULL;
  options->db_name = NULL;
  options->db_link = NULL;
  options->first_motifs = 0;
  options->e_max = EXPECT;
  options->m_thresh = MOTIF_THRESH_DEFAULT;
  options->min_ev = BIG;
  options->min_seqs = 0;
  options->adj_hit_pvalue = false;
  options->shuffle = false;
  options->sonly = false;
  options->lump = false;
  options->status = true;
  options->html = true;
  options->mast2txt = true;
  options->stype = Combine;
  options->translate_dna = false;
  options->best_motifs = false;
  options->use_seq_comp = false;
  options->rem_corr = false;
  options->hit_list = false;
  options->seed = mts_seed(&(options->prng)); // currently no way to set seed
  options->ids = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
  options->idxs = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);

  // process arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", mast_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_BFILE:
        options->bfile = optarg;
        break;
      case OPT_DBLIST:
        options->db_list = true;
        break;
      case OPT_O:
      case OPT_OC:
        options->clobber = (opt == OPT_OC);
        options->out_dir = optarg;
        break;
      case OPT_HITLIST:
        options->hit_list = true;
        break;
      case OPT_REMCORE:
        options->rem_corr = true;
        break;
      case OPT_M:
        rbtree_make(options->ids, optarg, NULL);
        break;
      case OPT_MI:
        idx = atoi(optarg);
        rbtree_make(options->idxs, &idx, NULL);
        break;
      case OPT_C:
        options->first_motifs = atoi(optarg);
        break;
      case OPT_MEV:
        options->min_ev = atof(optarg);
        break;
      case OPT_DIAG:
        options->diagram = optarg;
        break;
      case OPT_NORC:
        options->stype = Norc;
        break;
      case OPT_SEP:
        options->stype = Separate;
        break;
      case OPT_DNA:
        options->translate_dna = true;
        break;
      case OPT_COMP:
        options->use_seq_comp = true;
        break;
      case OPT_EV:
        options->e_max = atof(optarg);
        break;
      case OPT_MT:
        options->m_thresh = atof(optarg);
        break;
      case OPT_W:
       weak = true;
        break;
      case OPT_BEST:
        options->best_motifs = true;
        break;
      case OPT_SEQP:
        options->adj_hit_pvalue = true;
        break;
      case OPT_MF:
        options->mo_name = optarg;
        break;
      case OPT_DF:
        options->db_name = optarg;
        break;
      case OPT_DL:
        options->db_link = optarg;
        break;
      case OPT_MINSEQS:
        options->min_seqs = atoi(optarg);
        break;
      case OPT_NOSTATUS:
        options->status = false;
        break;
      case OPT_NOTEXT:
        options->mast2txt = false;
        break;
      case OPT_NOHTML:
        options->html = false;
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      // EXPERIMENTAL OPTIONS (not in usage message)
      case OPT_SHUFFLE: // shuffle columns of motifs
        options->shuffle = true;
        break;
      case OPT_SONLY: // use only spacing p-values in product
        options->sonly = true;
        break;
      case OPT_LUMP: // combine spacings into one p-value
        options->lump = true;
        break;
      case '?':
        usage(NULL);
        break;
      default:
        die("Unhandled option! Should never get here!");
    }
  }
  if (optind >= argc) usage("Motif file not specified");
  options->mo_source = argv[optind];
  optind++;
  if (optind >= argc) usage("Sequence file not specified");
  options->db_source = argv[optind];
  optind++;
  if (optind < argc) usage("Unused values specified");
  options->w_thresh = options->m_thresh * (weak ? 10 : 1);
  if (options->hit_list) options->e_max = BIG;
  if (strcmp(options->db_source, options->mo_source) == 0) {
    usage("The sequence database and motif file have the same source.\n");
  }
}

/**********************************************************************/
/*
  main routine
*/
/**********************************************************************/
int main(int argc, char *argv[])
{
  //feenableexcept(FE_DIVBYZERO | FE_INVALID);// FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW
  MAST_OPTIONS_T options;
  ARRAYLST_T *databases; // the list of seq databases 
  time_t mo_last_mod; // last modification date of motifs 
  ALPH_T *alph; // alphabet (as loaded from motifs)
  ALPH_T *seq_alph; // if translating then the source alph, otherwise just an alias for alph
  XLATE_T *xlate; // translator
  int i, j;
  int nseqs; // total sequences in database 
  double residues = 0; // total number of residues in seqs 
  ARRAYLST_T * all_motifs; // all motifs, not just the good ones 
  int file_count; // number of motifs in the file
  int used_count; // number of motifs after applying filters
  LO **los; // array of logodds matrices
  int nbad; // number of bad (correlated) motifs 
  double **pv = NULL; // p-value tables for each motif 
  ARRAY_T *bgfreqs, *motif_bg; // background letter probability distribution 

  (void) myclock(); // record CPU time 

  process_arguments(argc, argv, &options);

  if (strcmp(options.mo_source, "-") != 0) {
    //update the last mod date to the correct value
    struct stat stbuf;
    stat(options.mo_source, &stbuf);
    mo_last_mod = stbuf.st_mtime;
  } else {
    mo_last_mod = time(NULL);
  }

  xlate = NULL;
  if (options.translate_dna) {
    // create the DNA to protein translator
    xlate = xlate_dna2protein();
  }

  // load the motifs
  all_motifs = load_motifs(&options, xlate, &alph, &motif_bg, &file_count);
  used_count = arraylst_size(all_motifs);

  // load the background the MAST way; TLB added motif_bg 2-Feb-2017
  bgfreqs = load_background(alph, xlate != NULL, options.stype, motif_bg, &options.bfile);

  // convert the motifs to LO
  los = convert_2_log_odds(alph, xlate, bgfreqs, all_motifs, false, MAST_RANGE, NULL);

  // compute the pairwise motif correlations (similarities) 
  MATRIX_T *los_corr;
  los_corr = motif_corr(used_count, los);

  // get list of motifs that should be removed 
  nbad = mark_correlated_motifs(MAXCORR, all_motifs, los, los_corr);

  // remove highly correlated motifs from query 
  if (options.rem_corr && nbad) {
    for (i = 0, j = 0; i < used_count; i++) {
      if (!los[i]->is_bad) {
        los[j++] = los[i];
      }
    }
    used_count -= nbad;
  }

  // set global variable diagram in diagram.h
  diagram = options.diagram;
  // parse the ordering and spacing diagram (from globals in diagram.h)
  // this must happen after we've remove the unused motifs
  process_diagram(file_count, used_count, los);

  // get the sequence alphabet
  seq_alph = xlate != NULL ? xlate_src_alph(xlate) : alph;

  // Check that we can handle strands as requested
  if (!alph_has_complement(seq_alph)) { // protein database 
    if (options.stype == Separate || options.stype == Norc) {
      die("You may not specify -sep or -norc with a unstranded database.\n");
    }
    options.stype = Unstranded;
  } // database type 

  // calculate p-values of all integer score values in range [0...w*MAST_RANGE]
  pv = mm_malloc(sizeof(double*) * used_count);
  for (i = 0; i < used_count; i++) {
    pv[i] = calc_pssm_cdf(los[i]->w, los[i]->alen,
      MAST_RANGE, los[i]->logodds, bgfreqs);
    if (!pv[i]) {
      error("There is something wrong with motif %s\n", los[i]->meme_name);
    }
  }
  exit_on_error();

  databases = init_databases(&options, seq_alph);

  if (options.hit_list) {
    print_hit_list(argc, argv, &options, alph, xlate, los, used_count, pv, databases);
    return 0; // all done! 
  }

  // get the scores and p-values for all sequences in the database
  ARRAYLST_T *sequences = arraylst_create();
  ARRAYLST_T *hits = arraylst_create();
  for (i = 0, nseqs = 0; i < arraylst_size(databases); ++i) {
    DATABASE_T *db = (DATABASE_T*)arraylst_get(i, databases);
    get_scores(&options, alph, xlate, db, los, used_count, 
      pv, &nseqs, &residues, sequences, hits);
  }

  // bail if no sequences were read succesfully
  if (nseqs == 0) {
    die("Quitting due to errors or empty database.\n"); 
  }

  //sort by score
  arraylst_qsort(score_asc, hits);
  //remove anything with an evalue worse than the threshold
  for (i = arraylst_size(hits)-1; i >= 0; --i) {
    STRAND_T *st = arraylst_get(i, hits);
    double evalue = st->Pvalue * nseqs;
    if (evalue <= options.e_max) break;
  }
  if (++i < arraylst_size(hits)) {//TODO FIXME check if the hits are freed somewhere
    arraylst_remove_range(i, arraylst_size(hits) - i, NULL, hits);
  }

  //output
  output_mast_xml(
      &options, alph, xlate, bgfreqs, motif_bg,
      used_count, nbad, los, all_motifs, los_corr,
      mo_last_mod, 
      norder, order, space,
      nseqs, databases, hits, 
      argc, argv);

  // finish status report
  if (options.status) fprintf(stderr, "\n");

  //cleanup
  arraylst_destroy(sseq_destroy, sequences);
  arraylst_destroy(destroy_database, databases);
  for (i = 0; i < used_count; i++) {
    destroy_log_odds(los[i]);
  }
  free(los); // subset of all_los so no need to free contents
  free_matrix(los_corr);
  free_motifs(all_motifs);
  arraylst_destroy(NULL, hits);
  free_array(bgfreqs);
  free_array(motif_bg);
  if (xlate) xlate_destroy(xlate);
  free_2array(pv, used_count);
  alph_release(alph);
  rbtree_destroy(options.ids);
  rbtree_destroy(options.idxs);

  //output alternate formats
  bool output_error = false;
  if (options.out_dir != NULL && options.html) {
    // generate HTML
    char *prog;
    prog = get_meme_bin_file("mast_xml_to_html");
    if (prog != NULL) {
      STR_T *cmd;
      int ret;
      cmd = str_create(0);
      str_append2(cmd, prog);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, options.out_dir, XML_FILENAME);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, options.out_dir, HTML_FILENAME);

      ret = system(str_internal(cmd));

      if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
        output_error = true;
        fprintf(stderr, "Warning: mast_xml_to_html exited abnormally and may "
            "have failed to create HTML output.\n");
      }

      str_destroy(cmd, FALSE);
      free(prog);
    } else {
      output_error = true;
      fprintf(stderr, "Warning: could not find mast_xml_to_html. "
          "The HTML output could not be created.\n");
    }
  }
  if (options.out_dir != NULL && options.mast2txt) {
    // generate HTML
    char *prog;
    prog = get_meme_bin_file("mast_xml_to_txt");
    if (prog != NULL) {
      STR_T *cmd;
      int ret;
      cmd = str_create(0);
      str_append2(cmd, prog);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, options.out_dir, XML_FILENAME);
      str_append(cmd, " ", 1);
      str_append_path(cmd, 2, options.out_dir, TXT_FILENAME);

      ret = system(str_internal(cmd));

      if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
        output_error = true;
        fprintf(stderr, "Warning: mast_xml_to_txt exited abnormally and may "
            "have failed to create text output.\n");
      }

      str_destroy(cmd, FALSE);
      free(prog);
    } else {
      output_error = true;
      fprintf(stderr, "Warning: could not find mast_xml_to_txt. "
          "The text output could not be created.\n");
    }
  }

  return 0;
} // main 

