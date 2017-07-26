#include "array-list.h"
#include "momo.h"

// Structure for tracking momo command line algorithm.
typedef enum modl_op {
  INIT_STEP, // initial state
  MERGE,
  ADD,
  APPEND,
  REPLACE,
  REMOVE
} MODL_OP_T;

// If op == merge: remove mtfidx1, mtfidx2, merge and append results to end
// If op == add: use position and residue
// If op == append: add the results of merging mtfIdx1 and position/residue to end. Remove mtfidx2 (unless mtfidx is the last element of the motif. In this case, we do not remove a motif.
// If op == replace: do same as append, but remove mtfIdx1 as well.
// If op == remove: remove mtfidx1
typedef struct modl_step {
  MODL_OP_T op;
  int mtfIdx1;
  int mtfIdx2;
  int position;
  int residue;
  double score;
} MODL_STEP_T;

// This is only used to rank the motifs in the final motif set
// We rank the motifs in the increase in MDL when they are removed
typedef struct final_regex_motif_scores {
  REGEX_MOTIF_T* motif;
  double score;
} FINAL_REGEX_MOTIF_SCORES_T;

void create_modl_motifs(SUMMARY_T* summary,
                        MOMO_OPTIONS_T* options,
                        MOD_INFO_T * mod_info,
                        SEQ_T** all_sequences,
                        int num_sequences);

void cleanup_modl_step(MODL_STEP_T* modl_step);

/***********************************************************************
 Functions for output
 ***********************************************************************/
double remove_motif(ARRAYLST_T* regexmotifs,
                    MATRIX_T* bg_freqs,
                    int max_motifs,
                    int number,
                    MOMO_OPTIONS_T* options,
                    SUMMARY_T* summary,
                    MOD_INFO_T* mod_info,
                    BOOLEAN update);

void do_step(MODL_STEP_T* step,
             ARRAYLST_T* regexmotifs,
             MATRIX_T* bg_freqs,
             int max_motifs,
             MOMO_OPTIONS_T* options,
             SUMMARY_T* summary,
             MOD_INFO_T* mod_info);

char* regexmotif_to_string(REGEX_MOTIF_T* regexmotif,
                           MOD_INFO_T* mod_info,
                           SUMMARY_T* summary,
                           MOMO_OPTIONS_T* options);
