#include "array-list.h"
#include "momo.h"

// Structure for tracking momo command line algorithm.
typedef enum motifx_status {
  ACTIVE,
  INACTIVE,
  DELETED
} MOTIFX_STATUS_T;

typedef struct regexmotif {
  BOOLEAN_T* conserved;
  MATRIX_T* residues;
//  BOOLEAN_T* representation;
//  ARRAYLST_T* residues;
}  REGEX_MOTIF_T;

MATRIX_T* get_count_matrix(MATRIX_T* count,
                           ARRAYLST_T* sequences,
                           MOTIFX_STATUS_T** status,
                           MOMO_OPTIONS_T* options,
                           SUMMARY_T* summary);

void convert_bg_freqs_to_binomial(MATRIX_T* phospho_count,
                                  MATRIX_T* bg_freqs,
                                  int num_phospho_seqs);

void print_matrix_to_terminal(MATRIX_T* matrix, MOMO_OPTIONS_T* options, SUMMARY_T* summary);

void create_motifs(MOMO_OPTIONS_T* options,
                   SUMMARY_T* summary,
                   SEQ_T** all_sequences,
                   int num_sequences);