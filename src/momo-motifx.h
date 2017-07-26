#include "array-list.h"
#include "momo.h"

/**
 * Creates and stores a motifs using the motif-x
 * algorithm for a given mod until no more are left.
 */
void create_motifx_motifs(SUMMARY_T* summary,
                          MOMO_OPTIONS_T* options,
                          MOD_INFO_T * mod_info,
                          SEQ_T** all_sequences,
                          int num_sequences);
