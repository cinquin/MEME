/********************************************************************
 * MOMO Portal
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "display.h"
#include "display_globals.h"
#include "dir.h"
#include "momo.h"
#include "momo-output.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "binomial.h"
#include "meme.h"
#include "read_seq_file.h"
#include "momo-algorithm.h"

/* local variables */
#define DATA_HASH_SIZE 100003

/**
 * Creates a motif for a given mod using a simple frequency matrix.
 */
void create_simple_motif(SUMMARY_T* summary,
                         MOMO_OPTIONS_T* options,
                         MOD_INFO_T * mod_info) {
  int i;
  int j;
  
  const char* alph_letters = summary->alph_letters;
  
  // Create the frequency matrix
  MATRIX_T* freqs = NULL;
  freqs = get_count_matrix(freqs, mod_info->seq_list, NULL, options, summary);
  normalize_rows(0.0, freqs);
  
  // Create the motif
  MOTIF_INFO_T* motifinfo = mm_malloc(sizeof(MOTIF_INFO_T));
  motifinfo->motif = allocate_motif(mod_info->mod_name, "", summary->alph, freqs, NULL);
  motifinfo->seqs = arraylst_create();
  for (i = 0; i < arraylst_size(mod_info->seq_list); ++i) {
    SEQ_T* seqobject = options->eliminate_repeats ? hash_get_entry_value(arraylst_get(i, mod_info->seq_list)) : arraylst_get(i, mod_info->seq_list);
    arraylst_add(get_raw_sequence(seqobject), motifinfo->seqs);
  }
  motifinfo->fg_size = arraylst_size(mod_info->seq_list);
  arraylst_add(motifinfo, mod_info->motifinfos);
  
  // clean up
  free_matrix(freqs);
}
