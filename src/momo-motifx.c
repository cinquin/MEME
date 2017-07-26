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
#include "momo-motifx.h"

/* local variables */
#define DATA_HASH_SIZE 100003

/**
 * Given a bg frequency matrix and a phospho count matrix, we will convert the bg
 * freq matrix into a binomial matrix. No reason to keep a bg freq matrix, so we will recycle.
 */
void convert_bg_freqs_to_binomial(MATRIX_T* phospho_count, MATRIX_T* bg_freqs, int num_phospho_seqs) {
  int i;
  int j;
  for (i = 0; i < get_num_rows(bg_freqs); ++i) {
    for (j = 0; j < get_num_cols(bg_freqs); ++j) {
      double p = get_matrix_cell_defcheck(i, j, bg_freqs);
      double cxj = get_matrix_cell_defcheck(i, j, phospho_count);
      int n = num_phospho_seqs;
      double log_value = log_betai(cxj, n-cxj+1, p);
      set_matrix_cell_defcheck(i, j, log_value, bg_freqs);
    }
  }
}


/**
 * Given a matrix, find the row idx and column idx of the minimum value
 * that passes the score and count threshold. If none is found, then
 * row_idx, col_idx, and lowest will not be modified.
 */
void find_most_significant_within_matrix(MATRIX_T* binomial,
                                         MATRIX_T* phospho_count,
                                         int* row_idx,
                                         int* col_idx,
                                         double* lowest,
                                         MOMO_OPTIONS_T* options) {
  int i;
  int j;
  for (i = 0; i < get_num_rows(binomial); ++i) {
    for (j = 0; j < get_num_cols(binomial); ++j) {
      double binomial_value = get_matrix_cell_defcheck(i, j, binomial);
      double count_value = get_matrix_cell_defcheck(i, j, phospho_count);
      // 1. value cannot be nan/inf
      // 2. log(value) < log(options->score_threshold)
      // 3. value should not be the modified location
      if (i != options->width/2 && !isnan(binomial_value) && !isinf(binomial_value) && binomial_value < log(options->score_threshold) && count_value >= options->count_threshold) {
        if (binomial_value < *lowest || isnan(*lowest)) {
          *row_idx = i;
          *col_idx = j;
          *lowest = get_matrix_cell_defcheck(i, j, binomial);
        }
      }
      
    }
  }
}

/**
 * Remove sequences do not match a pattern from phospho and bg lists and update their respective count matrix
 */
void remove_sequences_and_update_matrix(char letter,
                                        int pos,
                                        ARRAYLST_T* seqs,
                                        MOTIFX_STATUS_T** status_array,
                                        int* num_active,
                                        MATRIX_T* count,
                                        SUMMARY_T* summary,
                                        MOMO_OPTIONS_T* options) {
  
  int i;
  const char* alph_letters = summary->alph_letters;
  
  // Look through phospho_seqs and remove sequences. Update phospho_seqs
  for (i = 0; i < arraylst_size(seqs); ++i) {
    char* curr_seq = get_raw_sequence((SEQ_T*) (options->eliminate_repeats ? hash_get_entry_value((HASH_TABLE_ENTRY*) arraylst_get(i, seqs)) : arraylst_get(i, seqs)));
    // For anything active that does not match the pattern, turn it inactive.
    MOTIFX_STATUS_T status = (*status_array)[i];
    if (status == ACTIVE && curr_seq[pos] != letter) {
      *num_active = *num_active - 1;
      (*status_array)[i] = INACTIVE;
    }
  }
  count = get_count_matrix(count, seqs, status_array, options, summary);
}

/**
 * Delete the active sequences and turn the inactive sequences to inactive.
 */
void delete_sequences(MOTIFX_STATUS_T** status_list, int size) {
  int i;
  for (i = 0; i < size; ++i) {
    MOTIFX_STATUS_T status = (*status_list)[i];
    if (status == ACTIVE) {
      (*status_list)[i] = DELETED;
    } else if (status == INACTIVE) {
      (*status_list)[i] = ACTIVE;
    }
  }
}

/**
 * Recursive function. While there is a significant position/residue pair,
 * will update pattern and return a count matrix for this pattern.
 */
MATRIX_T* add_to_pattern(char* pattern,
                         ARRAYLST_T* phospho_seqs,
                         ARRAYLST_T* bg_seqs,
                         MOTIFX_STATUS_T** phospho_status,
                         MOTIFX_STATUS_T** bg_status,
                         int* num_active,
                         int* num_bg_active,
                         MATRIX_T* phospho_count,
                         MATRIX_T* bg_count,
                         double* motif_score,
                         SUMMARY_T* summary,
                         MOMO_OPTIONS_T* options) {
  int i;
  const char* alph_letters = summary->alph_letters;
  
  // Set binomial to bg count, normalize to frequencies, and then convert it into a binomial matrix
  MATRIX_T* binomial = duplicate_matrix(bg_count);	// at this point it is a bg count matrix
  normalize_rows(0.0, binomial);                    // at this point it is a bg freq matrix
  convert_bg_freqs_to_binomial(phospho_count, binomial, *num_active);
  
  // Find minimum, and update if passes thresholds.
  int row_idx = -1;
  int col_idx = -1;
  double minimum = NAN;
  find_most_significant_within_matrix(binomial, phospho_count, &row_idx, &col_idx, &minimum, options);
  
  // clean up binomial after usage.
  free_matrix(binomial);
  
  // If one of the position/residue pairs pass, update pattern and statuses, then try to add another character to the pattern.
  if (row_idx >= 0 && col_idx >= 0) {
    pattern[row_idx] = alph_letters[col_idx];
    *motif_score = *motif_score - (minimum / log(10));
    
    // remove sequences an update matrices for phospho and bg
    remove_sequences_and_update_matrix(alph_letters[col_idx],
                                       row_idx,
                                       phospho_seqs,
                                       phospho_status,
                                       num_active,
                                       phospho_count,
                                       summary,
                                       options);
    remove_sequences_and_update_matrix(alph_letters[col_idx],
                                       row_idx,
                                       bg_seqs,
                                       bg_status,
                                       num_bg_active,
                                       bg_count,
                                       summary,
                                       options);
    
    return add_to_pattern(pattern,
                          phospho_seqs,
                          bg_seqs,
                          phospho_status,
                          bg_status,
                          num_active,
                          num_bg_active,
                          phospho_count,
                          bg_count,
                          motif_score,
                          summary,
                          options);
  }
  return phospho_count;
}

/**
 * Recursive function. Creates and stores a motif using the motif-x
 * algorithm until no more are left.
 */
void create_motifx_motif(ARRAYLST_T* phospho_seqs,
                         ARRAYLST_T* bg_seqs,
                         MOTIFX_STATUS_T** phospho_status,
                         MOTIFX_STATUS_T** bg_status,
                         MATRIX_T* phospho_count,
                         MATRIX_T* bg_count,
                         int* num_active,
                         int* num_bg_active,
                         char* modname,
                         MOD_INFO_T* mod_info,
                         MOMO_OPTIONS_T* options,
                         SUMMARY_T* summary) {
  int i;
  int j;
  
  const char* alph_letters = summary->alph_letters;
  
  // Initialize pattern, sequence count, bg sequence count, and overall score for this motif.
  char* pattern = mm_malloc(options->width + 1);
  for (i = 0; i < options->width; ++i) {
    pattern[i] = 'X';
  }
  pattern[options->width] = '\0';
  int* num_active_copy = mm_malloc(sizeof(int));
  *num_active_copy = *num_active;
  int* num_bg_active_copy = mm_malloc(sizeof(int));
  *num_bg_active_copy = *num_bg_active;
  double* motif_score = mm_malloc(sizeof(double));
  *motif_score = 0;
  
  // Set the pattern, num active copy, num bg active copy, motif score, and get a count of the sequences
  MATRIX_T* result_count_matrix = add_to_pattern(pattern, phospho_seqs, bg_seqs, phospho_status, bg_status, num_active_copy, num_bg_active_copy, phospho_count, bg_count, motif_score, summary, options);
  
  // If any of the characters are not X, then we have found a pattern
  BOOLEAN_T found_pattern = FALSE;
  for (i = 0; i < options->width; ++i) {
    if (pattern[i] != 'X') {
      found_pattern = TRUE;
    }
  }
  
  // If there is a pattern, store the pattern and call create_motifx_motif again.
  if (found_pattern) {
    // fill out the rest of the pattern (e.g. if you have pattern ..ASAAA, and realize the actual pattern is A.ASAAA
    for (i = 0; i < options->width; i++) {
      for (j = 0; j < strlen(alph_letters); j++) {
        if ((int) get_matrix_cell_defcheck(i, j, result_count_matrix) == *num_active_copy) {
          pattern[i] = alph_letters[j];
        }
      }
    }
    
    // create the pattern name
    char* pattern_name = mm_malloc(strlen(pattern) + strlen(modname) + 3);
    pattern_name[0] = '\0';
    strncat(pattern_name, pattern, strlen(pattern)/2);
    strncat(pattern_name, "_", 1);
    strncat(pattern_name, modname, strlen(modname));
    strncat(pattern_name, "_", 1);
    strncat(pattern_name, pattern + strlen(pattern)/2 + 1, strlen(pattern)/2);
    
    // convert this count matrix into frequencies
    normalize_rows(0.0, result_count_matrix);
    
    // Store this motif
    MOTIF_INFO_T* motifinfo = mm_malloc(sizeof(MOTIF_INFO_T));
    MOTIF_T* motif = allocate_motif(pattern_name, "", summary->alph, result_count_matrix, NULL);
    set_motif_nsites(motif, *num_active_copy);
    motifinfo->motif = motif;
    motifinfo->seqs = arraylst_create();
    motifinfo->score = *motif_score;
    motifinfo->fg_match = *num_active_copy;
    motifinfo->fg_size = *num_active;
    motifinfo->bg_match = *num_bg_active_copy;
    motifinfo->bg_size = *num_bg_active;
    for (i = 0; i < arraylst_size(phospho_seqs); ++i) {
      MOTIFX_STATUS_T status = (*phospho_status)[i];
      if (status == ACTIVE) {
        SEQ_T* active_sequence = (options->eliminate_repeats) ? hash_get_entry_value(arraylst_get(i, phospho_seqs)) : arraylst_get(i, phospho_seqs);
        arraylst_add(get_raw_sequence(active_sequence), motifinfo->seqs);
      }
    }
    arraylst_add(motifinfo, mod_info->motifinfos);
    
    // delete the sequences from this motif. turn inactive into active.
    delete_sequences(phospho_status, arraylst_size(phospho_seqs));
    delete_sequences(bg_status, arraylst_size(bg_seqs));
    
    // update the count of number of actives
    *num_active = *num_active - *num_active_copy;
    *num_bg_active = *num_bg_active - *num_bg_active_copy;
    
    // recalculate phospho count and bg count.
    phospho_count = get_count_matrix(phospho_count, phospho_seqs, phospho_status, options, summary);
    bg_count = get_count_matrix(bg_count, bg_seqs, bg_status, options, summary);
    
    // free up space
    myfree(pattern);
    myfree(num_active_copy);
    myfree(num_bg_active_copy);
    myfree(motif_score);
    myfree(pattern_name);
    
    // try to create another motif.
    create_motifx_motif(phospho_seqs,
                        bg_seqs,
                        phospho_status,
                        bg_status,
                        phospho_count,
                        bg_count,
                        num_active,
                        num_bg_active,
                        modname,
                        mod_info,
                        options,
                        summary);
  }
  // free up space
  myfree(pattern);
  myfree(num_active_copy);
  myfree(num_bg_active_copy);
  myfree(motif_score);
}

/**
 * Initializes and returns a MOTIFX_STATUS_T array of the given size
 */
MOTIFX_STATUS_T* init_status_array(int size) {
  int i;
  MOTIFX_STATUS_T* result = (MOTIFX_STATUS_T*) mm_malloc(sizeof(MOTIFX_STATUS_T) * size);
  for (i = 0; i < size; ++i) {
    result[i] = ACTIVE;
  }
  return result;
}

/**
 * Creates and stores a motifs using the motif-x
 * algorithm for a given mod until no more are left.
 */
void create_motifx_motifs(SUMMARY_T* summary,
                          MOMO_OPTIONS_T* options,
                          MOD_INFO_T * mod_info,
                          SEQ_T** all_sequences,
                          int num_sequences) {
  
  int i;
  ARRAYLST_T* phospho_seqs = mod_info->seq_list;
  ARRAYLST_T* bg_seqs = mod_info->bg_seq_list;
  
  // Initialize status
  MOTIFX_STATUS_T* phospho_status = init_status_array(arraylst_size(phospho_seqs));
  MOTIFX_STATUS_T* bg_status = init_status_array(arraylst_size(bg_seqs));
  
  // Initialize sequence counter. This contains the number of ACTIVE sequences in status
  int* num_active = mm_malloc(sizeof(int));
  int* num_bg_active = mm_malloc(sizeof(int));
  *num_active = arraylst_size(phospho_seqs);
  *num_bg_active = arraylst_size(bg_seqs);
  
  // Initialize count matrices
  MATRIX_T* phospho_count = NULL;
  MATRIX_T* bg_count = NULL;
  phospho_count = get_count_matrix(phospho_count, phospho_seqs, NULL, options, summary);
  bg_count = get_count_matrix(bg_count, bg_seqs, NULL, options, summary);
  
  // Recursively create motifs
  create_motifx_motif(phospho_seqs,
                      bg_seqs,
                      &phospho_status,
                      &bg_status,
                      phospho_count,
                      bg_count,
                      num_active,
                      num_bg_active,
                      mod_info->mod_name,
                      mod_info,
                      options,
                      summary);
  
  // Clean up
  myfree(num_active);
  myfree(num_bg_active);
  myfree(phospho_status);
  myfree(bg_status);
  free_matrix(phospho_count);
  free_matrix(bg_count);
}
