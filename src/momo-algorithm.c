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
#include "momo-simple.h"
#include "momo-motifx.h"
#include "momo-modl.h"

/* local variables */
#define DATA_HASH_SIZE 100003

/**
 * Using an array list of sequences (or hash table entries), creates and returns a count matrix.
 */
MATRIX_T* get_count_matrix(MATRIX_T* count,
                           ARRAYLST_T* sequences,
                           MOTIFX_STATUS_T** status,
                           MOMO_OPTIONS_T* options,
                           SUMMARY_T* summary) {
  int i;
  int j;
  
  BOOLEAN_T sequences_are_hash_entries = options->eliminate_repeats;
  const char* alph_letters = summary->alph_letters;
  
  if (!count) {
    count = allocate_matrix(options->width, strlen(alph_letters));
  }
  init_matrix(0.0, count);
  
  for (i = 0; i < arraylst_size(sequences); ++i) {
    if (!status || (*status)[i] == ACTIVE) {
      SEQ_T* seqobject = sequences_are_hash_entries ? hash_get_entry_value(arraylst_get(i, sequences)) : arraylst_get(i, sequences);
      char* raw_sequence = get_raw_sequence(seqobject);
      
      for (j = 0; j < options->width; j++) {
        if (strchr(alph_letters, raw_sequence[j])) {
          int aa_idx = strchr(alph_letters, raw_sequence[j]) - alph_letters;
          set_matrix_cell_defcheck(j, aa_idx, get_matrix_cell_defcheck(j, aa_idx, count) + 1.0, count);
        }
      }
    }
  }
  return count;
}

/**
 * Prints a matrix to the terminal. For debugging use.
 */
void print_matrix_to_terminal(MATRIX_T* matrix, MOMO_OPTIONS_T* options, SUMMARY_T* summary) {
  const char* alph_letters = summary->alph_letters;
  int i;
  int j;
  for (i = 0; i < options->width-1; i++) {
    for (j = 0; j < strlen(alph_letters); j++) {
      printf("%.1f\t", get_matrix_cell_defcheck(i, j, matrix));
    }
    printf("\n");
  }
  printf("\n");
}


/**
 * For each mod inside the mod table that passes all the
 * filters specified by the user, create a motif. This also
 * sets the number of (passing) mods/modtypes.
 */
void create_motifs(MOMO_OPTIONS_T* options,
                          SUMMARY_T* summary,
                          SEQ_T** all_sequences,
                          int num_sequences) {
  
  ALPH_T* alph = summary->alph;
  const char* alph_letters = summary->alph_letters;
  int i;
  
  // Get the mod table
  HASH_TABLE mod_table = summary->mod_table;
  ARRAYLST_T * mod_table_keys = summary->mod_table_keys;
  
  // Initialize statistic counters.
  unsigned long num_mods = 0;
  unsigned long num_types = 0;
  unsigned long pass_mods = 0;
  unsigned long pass_types = 0;
  
  // For each mod inside the mod table
  for (i = 0; i < arraylst_size(mod_table_keys); i++) {
    // Get the mod entry
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    mod_entry->mod_name = hash_get_entry_key(hash_entry);
    
    // Count the number of mods and mod types and add to overall counter
    unsigned long num_occurrences = mod_entry->mod_occurrences;
    num_mods += num_occurrences;
    num_types++;
    
    // For each passing mod, we will add to passing mod and mod type counter and create a motif
    unsigned long num_passing = arraylst_size(mod_entry->seq_list);
    if (num_passing >= options->min_occurrences) {
      // Increment counters
      pass_mods += num_passing;
      pass_types++;
      
      // Create motif
      if (options->algorithm == simple) {
        create_simple_motif(summary, options, mod_entry);
      } else if (options->algorithm == motifx) {
        create_motifx_motifs(summary, options, mod_entry, all_sequences, num_sequences);
      } else if (options->algorithm == modl) {
        create_modl_motifs(summary, options, mod_entry, all_sequences, num_sequences);
      }
    }
  }
  
  summary->num_mod = num_mods;
  summary->num_modtype = num_types;
  summary->num_mod_passing = pass_mods;
  summary->num_modtype_passing = pass_types;
}
