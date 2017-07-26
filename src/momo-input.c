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
#include "fasta-io.h"
#include "momo.h"
#include "array-list.h"
#include "io.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "read_seq_file.h"
#include "momo-input.h"
#include "prealigned-io.h"
#include "momo-algorithm.h"
#include "momo-modl.h"


/* local variables */
#define DATA_HASH_SIZE 100003

/***********************************************************************
 Initialize mod info
 ***********************************************************************/
MOD_INFO_T* init_modinfo(MOMO_OPTIONS_T* options, SUMMARY_T* summary) {
  int i;
  const char* alph_letters = summary->alph_letters;
  
  MOD_INFO_T * modinfo = mm_malloc(sizeof(MOD_INFO_T));
  
  // Create and init an entry inside the mod table
  modinfo->mod_occurrences = 1;
  modinfo->seq_table = (options->eliminate_repeats) ? hash_create(DATA_HASH_SIZE) : NULL;
  modinfo->seq_list = arraylst_create();
  modinfo->motifinfos = arraylst_create();
  modinfo->amino_acids = (BOOLEAN_T*) mm_malloc(sizeof(BOOLEAN) * strlen(alph_letters));
  for (i = 0; i < strlen(alph_letters); ++i) {
    modinfo->amino_acids[i] = false;
  }
  modinfo->mod_name = NULL;
  modinfo->bg_seq_table = NULL;
  modinfo->bg_seq_list = NULL;
  
  modinfo->modl_ops = NULL;
  
  return modinfo;
}

/***********************************************************************
 Clean up mod info
 ***********************************************************************/
void cleanup_modinfo(MOD_INFO_T* modinfo) {
  int i;
  // Each MOD_INFO_T contains: motifinfos, mod_name, mod_occurrences, seq_list/seq_table, bg_sequences, amino_acids
  
  // Free motifinfos
  // Each MOTIF_INFO_T contains: motif, seqs (these point to the inside the seq_list/seq_table)
  for (i = 0; i < arraylst_size(modinfo->motifinfos); ++i) {
    MOTIF_INFO_T* motifinfo = arraylst_get(i, modinfo->motifinfos);
    destroy_motif(motifinfo->motif);
    arraylst_destroy(NULL, motifinfo->seqs);
    myfree(motifinfo);
  }
  arraylst_destroy(NULL, modinfo->motifinfos);
  
  // Free sequences from seq_list/seq_table
  for (i = 0; i < arraylst_size(modinfo->seq_list); ++i) {
    free_seq(modinfo->seq_table == NULL ? arraylst_get(i, modinfo->seq_list) : hash_get_entry_value(arraylst_get(i, modinfo->seq_list)));
  }

  arraylst_destroy(NULL, modinfo->seq_list);
  if (modinfo->seq_table != NULL) {
    hash_destroy(modinfo->seq_table);
  }
  
  // clean up bg sequences
  if (modinfo->bg_seq_list != NULL) {
    for (i=0; i < arraylst_size(modinfo->bg_seq_list); ++i) {
      free_seq(modinfo->bg_seq_table == NULL ? arraylst_get(i, modinfo->bg_seq_list) : hash_get_entry_value(arraylst_get(i, modinfo->bg_seq_list)));
    }
    
    arraylst_destroy(NULL, modinfo->bg_seq_list);
    if (modinfo->bg_seq_table != NULL) {
      hash_destroy(modinfo->bg_seq_table);
    }
  }
  myfree(modinfo->amino_acids);
  
  // clean up modl ops
  if (modinfo->modl_ops != NULL) {
    for (i = 0; i < arraylst_size(modinfo->modl_ops); ++i) {
      MODL_STEP_T* modl_step = arraylst_get(i, modinfo->modl_ops);
      cleanup_modl_step(modl_step);
    }
    arraylst_destroy(NULL, modinfo->modl_ops);
  }
  myfree(modinfo);
}

/**
 * If a protein database is provided, parse the sequences
 */
void read_protein_database_sequences(MOMO_OPTIONS_T* options,
                                     SUMMARY_T* summary,
                                     SEQ_T*** all_sequences,
                                     int* num_sequences) {
  char* protein_database_filename = options->protein_database_filename;
  FILETYPE_T bg_filetype = options->bg_filetype;
  ALPH_T* alph = summary->alph;
  
  if (protein_database_filename) {
    FILE *fasta_fp;
    fasta_fp = fopen(protein_database_filename, "r");
    if (options->bg_filetype == fasta) {
      read_many_fastas(alph, fasta_fp, MAX_SEQ, num_sequences, all_sequences);
    } else { // else prealigned
      read_many_line_delimited_sequences(alph, fasta_fp, MAX_SEQ, num_sequences, options->width, all_sequences);
    }
    fclose(fasta_fp);
  }
}

/**
 * Given a sequence and a start idx, we will change result to contain a peptide
 * of the length k. If the idx is is negative or greater than the length
 * of the sequence, we will replace with X.
 */
void set_kmer(char** result, char* sequence, int start_idx, int k) {
  int i;
  (*result)[k] = '\0';
  for (i = 0; i < k; ++i) {
    int curr_idx = start_idx + i; // start_idx + i
    if (curr_idx < 0 || curr_idx > strlen(sequence) - 1) {
      (*result)[i] = 'X';
    } else {
      (*result)[i] = sequence[curr_idx];
    }
  }
}

/**
 * Sets the O(1) lookup table. This hashes from every unique kmer within
 * the protein database to an arraylist of KMER_LOC_T objects.
 */
void create_hash_fasta_preprocess_table(MOMO_OPTIONS_T* options,
                                        SUMMARY_T* summary,
                                        SEQ_T** all_sequences,
                                        int num_sequences) {
  int i;
  int j;
  int k;
  
  if (options->hash_fasta && options->protein_database_filename) {
    HASH_TABLE hash_fasta_table = hash_create(DATA_HASH_SIZE);	// hash table of mods
    ARRAYLST_T * hash_fasta_table_keys = arraylst_create(); // hash table keys
    int hash_fasta_width = options->hash_fasta_width;
    char* kmer = mm_malloc(hash_fasta_width + 1);
    for (i = 0; i < num_sequences; i++) {
      char* raw_sequence = get_raw_sequence(all_sequences[i]);
      int raw_sequence_length = strlen(raw_sequence);
      for (j = 0; j < raw_sequence_length - hash_fasta_width + 1; ++j) {
        set_kmer(&kmer, raw_sequence, j, hash_fasta_width);
        KMER_LOC_T * kmerloc = mm_malloc(sizeof(KMER_LOC_T));
        kmerloc->seqidx = i;
        kmerloc->kmeridx = j;
        if (!hash_lookup_str(kmer, hash_fasta_table)) {
          ARRAYLST_T * locations = arraylst_create();
          arraylst_add(kmerloc, locations);
          hash_insert_str_value(kmer, locations, hash_fasta_table);
          HASH_TABLE_ENTRY * hash_entry = hash_lookup_str(kmer, hash_fasta_table);
          arraylst_add(hash_entry, hash_fasta_table_keys);
        } else {
          HASH_TABLE_ENTRY * hash_entry = hash_lookup_str(kmer, hash_fasta_table);
          ARRAYLST_T * locations = hash_get_entry_value(hash_entry);
          arraylst_add(kmerloc, locations);
        }
      }
    }
    summary->hash_fasta_table = hash_fasta_table;
    summary->hash_fasta_table_keys = hash_fasta_table_keys;
    
    // Clean up
    myfree(kmer);
  }
}

/**
 * Returns whether the given value in the filter field passes
 * the filter threshold specified by the user
 */
BOOLEAN_T check_filter_threshold(char* value_str,
                                 MOMO_OPTIONS_T* options) {
  if (!options->filter) {
    return true;
  } else {
    FILTERTYPE_T filter_type = options->filter_type;
    double filter_threshold = options->filter_threshold;
    double value = atof(value_str);
    if (filter_type == le) {
      return (value <= filter_threshold);
    } else if (filter_type == lt) {
      return (value < filter_threshold);
    } else if (filter_type == eq) {
      return (value == filter_threshold);
    } else if (filter_type == gt) {
      return (value > filter_threshold);
    } else { // filter_type = ge
      return (value >= filter_threshold);
    }
  }
}

/**
 * Given a peptide (modless), we want to see if we can find the protein
 * it originated from. If not, we will return the peptide.
 */
void find_fasta(SUMMARY_T* summary,
                 MOMO_OPTIONS_T* options,
                 SEQ_T **all_sequences,
                 int num_sequences,
                 char* modless,
                 KMER_LOC_T* result) {
  int i;
  int j;
  int hash_fasta_width = options->hash_fasta_width;
  HASH_TABLE hash_fasta_table = summary->hash_fasta_table;
  
  // Use the O(1) lookup table to find the location of the peptide within the protein database.
  // Returns a kmer info object if found.
  if (options->hash_fasta && options->protein_database_filename && strlen(modless) >= hash_fasta_width) {
    char* subseq = mm_malloc(hash_fasta_width + 1);
    for (i = 0; i < strlen(modless) - hash_fasta_width + 1; ++i) {
      set_kmer(&subseq, modless, 0, hash_fasta_width);
      HASH_TABLE_ENTRY * hash_entry = hash_lookup_str(subseq, hash_fasta_table);
      if (hash_entry) {
        ARRAYLST_T* kmerloclist = (ARRAYLST_T*) hash_get_entry_value(hash_entry);
        for (j = 0; j < arraylst_size(kmerloclist); ++j) {
          KMER_LOC_T* kmerloc = arraylst_get(j, kmerloclist);
          char* protein = get_raw_sequence(all_sequences[kmerloc->seqidx]);
          char* peptide = protein + kmerloc->kmeridx;
          if (strncmp(modless, peptide, strlen(modless)) == 0) {
            myfree(subseq);
            result->seqidx = kmerloc->seqidx;
            result->kmeridx = kmerloc->kmeridx;
            return;
          }
        }
      }
    }
    myfree(subseq);
  }
      
  // Use linear search to find the location of a peptide within
  // the protein database. Returns a kmer info object if found.
  for (i = 0; i < num_sequences; i++) {
    char* raw_sequence = get_raw_sequence(all_sequences[i]);
    char* found = strstr(raw_sequence, modless);
    if (found) {
      result->seqidx = i;
      result->kmeridx = found - raw_sequence;
      return;
    }
  }
}

/**
 * Adds a kmer to a sequence list (and table if eliminate repeats) can refer to either bg or fg sequences
 */
void add_kmer_to_sequence_list(HASH_TABLE* seq_table, ARRAYLST_T** seq_list, char* sequence, int idx, MOMO_OPTIONS_T* options) {
  // Set up variables
  int motif_width = options->width;
  int eliminate_repeat_width = options->eliminate_repeat_width;
  
  // Create eliminate_repeat_str and motif_str. This is the mod
  // along with its flanking amino acids for the respective length.
  char* eliminate_repeat_str = mm_malloc(eliminate_repeat_width + 1);
  set_kmer(&eliminate_repeat_str, sequence, idx - eliminate_repeat_width/2, eliminate_repeat_width);
  char* motif_str = mm_malloc(motif_width + 1);
  set_kmer(&motif_str, sequence, idx - motif_width/2, motif_width);
  
  BOOLEAN_T contains_unknowns = strchr(motif_str, 'X');
  
  // Store the sequence.
  if (!(options->remove_unknowns && contains_unknowns)) {
    SEQ_T* newsequence = allocate_seq("tempname", "tempdescription", 0, motif_str);
    if (options->eliminate_repeats) {
      if (!hash_lookup_str(eliminate_repeat_str, *seq_table)) {
        hash_insert_str_value(eliminate_repeat_str, newsequence, *seq_table);
        HASH_TABLE_ENTRY* hash_entry = hash_lookup_str(eliminate_repeat_str, *seq_table);
        arraylst_add(hash_entry, *seq_list);
      } else {
        free_seq(newsequence);
      }
    } else {
      arraylst_add(newsequence, *seq_list);
    }
  }
  
  // Clean Up
  myfree(eliminate_repeat_str);
  myfree(motif_str);
}

/**
 * Add a single mod with its flanking amino acids to the mod table
 */
void add_mod_to_mod_list(HASH_TABLE mod_table,
                    char* modless,
                    int mod_index,
                    char* mod_name,
                    char* protein_id,
                    char* scan,
                    SEQ_T** all_sequences,
                    int num_sequences,
                    MOMO_OPTIONS_T* options,
                    SUMMARY_T* summary,
                    BOOLEAN_T passes_threshold) {
  int i;
  
  // Set up variables
  MOD_INFO_T * modinfo = hash_get_entry_value(hash_lookup_str(mod_name, mod_table));
  int motif_width = options->width;
  int eliminate_repeat_width = options->eliminate_repeat_width;
  
  // For this mod, we will create a kmerinfo object.
  KMER_LOC_T* kmerloc = mm_malloc(sizeof(kmerloc));
  kmerloc->seqidx = -1;
  kmerloc->kmeridx = -1;
  find_fasta(summary, options, all_sequences, num_sequences, modless, kmerloc);
  
  char* protein = (kmerloc->seqidx >= 0) ? get_raw_sequence(all_sequences[kmerloc->seqidx]) : modless;
  int start_of_mod = (kmerloc->kmeridx >= 0) ? kmerloc->kmeridx + mod_index : mod_index;
  
  if (passes_threshold) {
    add_kmer_to_sequence_list(&modinfo->seq_table, &modinfo->seq_list, protein, start_of_mod, options);
  }
  free(kmerloc);
}

/**
 * Parse a PSM header line. We want to locat the idx of filter_value, sequence, scan, and protein_id
 */
void parse_tsv_header(MOMO_OPTIONS_T* options,
                      int* index_of_filter,
                      int* index_of_sequence,
                      int* index_of_scan,
                      int* index_of_protein,
                      FILE* fp, char** header,
                      size_t* len) {
  int i = 0;
  
  getline(header, len, fp);            // read the first line, this should be the header
  *header = strtok(*header, "\n");     // strip new line
  
  char *token;
  token = strtok(*header, "\t");
  while(token != NULL) {
    if (strcmp(token, options->sequence_column) == 0) {
      *index_of_sequence = i;
    }
    if (options->filter && strcmp(token, options->filter_field) == 0) {
      *index_of_filter = i;
    }
    if (strcmp(token, "scan") == 0) {
      *index_of_scan = i;
    }
    if (strcmp(token, "protein id") == 0) {
      *index_of_protein = i;
    }
    token = strtok(NULL, "\t");
    i++;
  }
  
  // ensure that input file has required columns
  if (*index_of_sequence == -1) {
    die("Could not find a sequence column inside the psm file.");
  }
/*  if (*index_of_scan == -1) {
    die("Could not find scan number inside the psm file.");
  }
  if (*index_of_protein == -1) {
    die("Could not find a protein id column inside the psm file.");
  }*/
  if (options->filter && *index_of_filter == -1) {
    die("Could not find the filter column inside the psm file.");
  }
}

/**
 * Parse a tsv body line into filter_value, sequence, scan, and protein_id
 */
void parse_tsv_bodyline(MOMO_OPTIONS_T* options,
                        char** filter_value,
                        char** sequence,
                        char** scan,
                        char** protein_id,
                        int index_of_filter,
                        int index_of_sequence,
                        int index_of_scan,
                        int index_of_protein,
                        char* line) {
  char* token;
  token = strtok(line, "\t");
  int i = 0;
  while (token != NULL) {
    if (i == index_of_sequence) {
      *sequence = strdup(token);
    }
    if (i == index_of_scan) {
      *scan = strdup(token);
    }
    if (i == index_of_protein) {
      *protein_id = strdup(token);
    }
    if (options->filter && i == index_of_filter) {
      *filter_value = strdup(token);
    }
    token = strtok(NULL, "\t");
    i++;
  }
}

/**
 * Given a sequence, duplicate it into modless sequence and fill
 * an array with each position indicating whether or not it has a mod. 
 * If yes, then will contain the mod name. Otherwise, it is set to NULL.
 */
void get_prealigned_mod_array(ARRAYLST_T* mod_array, char** modless,  char* sequence) {
  *modless = strdup(sequence);
  char* mod = mm_malloc(2);
  mod[0] = sequence[strlen(sequence)/2];
  mod[1] = '\0';
  arraylst_put(strlen(sequence)/2, mod, mod_array);
}

/**
 * Given a modified sequence, set the modless sequence and fill
 * an array with each position indicating whether or not it has a mod.
 * If yes, then will contain the mod name. Otherwise, it is set to NULL.
 */
void get_tsv_mod_array(ARRAYLST_T* mod_array, char** modless, char* sequence) {
  size_t seq_len = strlen(sequence);
  *modless = mm_malloc(seq_len + 1);
  *modless[0] = '\0';
  
  char* sequence_mod_copy = strdup(sequence);
  char* sequence_mod_iterator = sequence_mod_copy;
  char* ptr1 = strchr(sequence_mod_iterator, '[');
  char* ptr2 = strchr(sequence_mod_iterator, ']');
  
  // while mod has been found
  while (ptr1 && ptr2) {
    // find the length of the mod. For example, S[79.0] = 5 since we disregard the brackets.
    int mod_len = ptr2-ptr1;
    
    // create a copy of the mod
    char* mod = mm_malloc(mod_len + 1);
    strncpy(mod, ptr1-1, 1);
    strncpy(mod+1, ptr1+1, mod_len - 1);
    mod[mod_len] = '\0';
    
    // update the modless sequence
    strncat(*modless, sequence_mod_iterator, ptr1-sequence_mod_iterator);
    
    // store the mod at its proper index.
    arraylst_put(strlen(*modless) - 1, mod, mod_array);
    
    //update
    sequence_mod_iterator = ptr2+1;
    ptr1 = strchr(sequence_mod_iterator, '[');
    ptr2 = strchr(sequence_mod_iterator, ']');
  }
  // add the remaining amino acids to the end of modless
  strncat(*modless, sequence_mod_iterator, strlen(sequence_mod_iterator));
  
  myfree(sequence_mod_copy);
}

/**
 * Given a single phosphorylation file, add all the modifications and respective flanking
 * amino acid information into the mod table.
 */
void add_phospho_file_to_table(char* phospho_filename,
                                  MOMO_OPTIONS_T* options,
                                  SUMMARY_T* summary,
                                  SEQ_T** all_sequences,
                                  int num_sequences) {
  
  int i = 0;
  int j = 0;
  const char* alph_letters = summary->alph_letters;
  ALPH_T* alph = summary->alph;
  
  // open and prepare file for reading
  FILE* fp;
  char* line = NULL;
  size_t len = 0;
  ssize_t read;
  fp = fopen(phospho_filename, "r");
  if (fp == NULL) {
    die("Failed to read filename: %s", phospho_filename);
  }
  
  // parse header to find the indices of required data
  int index_of_filter = -1;           // the tab-delimited index of filter;
  int index_of_sequence = -1;         // the tab-delimited index of column sequence
  int index_of_scan = -1;             // the tab-delimited index of scan;
  int index_of_protein = -1;          // the tab-delimited index of protein id;

  if (options->fg_filetype == psm) {
    parse_tsv_header(options, &index_of_filter, &index_of_sequence, &index_of_scan, &index_of_protein, fp, &line, &len);
  }
  
  HASH_TABLE mod_table = summary->mod_table;
  ARRAYLST_T * mod_table_keys = summary->mod_table_keys;
  
  // for each body line (PSM) inside the file:
  while ((read = getline(&line, &len, fp)) != -1) {
    // strip new line and make a copy
    line = strtok(line, "\n");
    
    char* unprocessed_sequence = NULL;
    char* sequence = NULL;
    char* filter_value = NULL;
    
    // These are currently not used.
    char* scan = NULL;
    char* protein_id = NULL;
    
    if (options->fg_filetype == psm) {
      parse_tsv_bodyline(options, &filter_value, &unprocessed_sequence, &scan, &protein_id, index_of_filter, index_of_sequence, index_of_scan, index_of_protein, line);
    } else { // prealigned or prealigned fasta
      if (!(options->fg_filetype == fasta && line != NULL && line[0] == '>')) {
        parse_prealigned_bodyline(&unprocessed_sequence, line, options->width);
      }
    }
    // this means that this line does not have a sequence (e.g. fasta id lines)
    if (unprocessed_sequence == NULL) {
      continue;
    }
    
    if (options->fg_filetype != psm) {
      sequence = strdup(unprocessed_sequence);
    } else {
      // PSMs may appear in different formats, we will modify them to be in standard crux
      if (strlen(unprocessed_sequence) > 2 && unprocessed_sequence[1] == '.') {
        // sequences are in A.AS[79.9]A.- format, remove the dots and '-' from ends.
        sequence = mm_malloc(strlen(unprocessed_sequence) + 1);
        sequence[0] = '\0';
        if (unprocessed_sequence[0] != '-') {
          strncat(sequence, unprocessed_sequence, 1);
        }
        int pep_length = strlen(unprocessed_sequence) - 4;
        strncat(sequence, unprocessed_sequence + 2, pep_length);
        if (unprocessed_sequence[strlen(unprocessed_sequence) - 1] != '-') {
          strncat(sequence, unprocessed_sequence + strlen(unprocessed_sequence) - 1, 1);
        }
      } else {
        int count = 0;
        int plusminus_idx = -1;
        for (i = 0; i < strlen(unprocessed_sequence); ++i) {
          if (unprocessed_sequence[i] == '+' || unprocessed_sequence[i] == '-') {
            count++;
            if (plusminus_idx == -1) {
              plusminus_idx = i;
            }
          }
        }
        if (plusminus_idx != -1 && unprocessed_sequence[plusminus_idx-1] != '[') {
          // sequences are in AS+79.9A format, sequences are missing brackets, create bracket rep
          sequence = mm_malloc(strlen(unprocessed_sequence) + 2*count + 1);
          sequence[0] = '\0';
          for (i = 0; i < strlen(unprocessed_sequence); ++i) {
            if (unprocessed_sequence[i] == '+' || unprocessed_sequence[i] == '-') {
              strncat(sequence, "[", 1);
            }
            strncat(sequence, unprocessed_sequence + i, 1);
            if (!isalpha(unprocessed_sequence[i]) && isalpha(unprocessed_sequence[i+1])) {
              strncat(sequence, "]", 1);
            }
          }
        } else {
          // sequence is fine as it is
          sequence = strdup(unprocessed_sequence);
        }
      }
    }
    
    // check whether or not this psm passes the filter threshold.
    BOOLEAN_T passes_threshold = options->fg_filetype == prealigned || options->fg_filetype == fasta || check_filter_threshold(filter_value, options);
    
    // NOTE: mod_name = X[mod_mass] when we create a mod per aminoacid + mod pair,
    // and simply [mod_mass] when we create a single mod per mass.
    
    // we will now create a modless version of sequence, while:
    // 1. update mod_table / mod_table_keys. This is a hash table
    // mapping from mod_name => MOD_INFO_T.
    // 2. create an arraylist: This stores the mod_name of mods that pass the threshold
    // at the index they are located inside the modless string.
    
    size_t seq_len = strlen(sequence);
    // This is supposed to store the locations of the mods within the sequence.
    ARRAYLST_T * mod_array = arraylst_create_sized(seq_len);
    for (i = 0; i < seq_len; i++) {
      arraylst_add(NULL, mod_array);
    }
    // modless sequence
    char* modless = NULL;
    
    // note that after we obtain the mods, we need to free them.
    if (options->fg_filetype == psm) {
      // sequence contains multiple mods
      get_tsv_mod_array(mod_array, &modless, sequence);
    } else { // fg_filetype = prealigned or fasta
      // sequence only contains a single mod
      get_prealigned_mod_array(mod_array, &modless, sequence);
    }
    
    // update the tables for each valid mod:
    for (i = 0; i < seq_len; i++) {
      char* mod = arraylst_get(i, mod_array);
      if (mod != NULL) {
        char* mod_name = (options->single_motif_per_mass) ? mod + 1 : mod; // e.g. 79.9 or S79.9
        HASH_TABLE_ENTRY* hash_entry = hash_lookup_str(mod_name, mod_table);
        // update the mod_table for each mod.
        if (!hash_entry) {
          MOD_INFO_T* newmodinfo = init_modinfo(options, summary);
          
          // we are adding an entry to the hash table
          hash_insert_str_value(mod_name, (void *) newmodinfo, mod_table);
          hash_entry = hash_lookup_str(mod_name, mod_table);
          arraylst_add(hash_entry, mod_table_keys);
        } else {
          // Update the entry within the mod info
          MOD_INFO_T * currmodinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
          currmodinfo->mod_occurrences++;
        }
        // update the list of amino acids associated with this mod
        MOD_INFO_T * currmodinfo = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
        currmodinfo->amino_acids[strchr(alph_letters, mod[0]) - alph_letters] = true;
        
        // find the flanks of this mod and add it to the list of sequences associated with the mod name
        add_mod_to_mod_list(mod_table, modless, i, mod_name, protein_id, scan, all_sequences, num_sequences, options, summary, passes_threshold);
      }
    }
    
    // clean up
    arraylst_destroy(free, mod_array);
    myfree(modless);
    
    myfree(unprocessed_sequence);
    myfree(sequence);
    myfree(filter_value);
    
    // These are currently not used.
    myfree(scan);
    myfree(protein_id);
    
  }
  
  // clean up
  myfree(line);
  fclose(fp);
}

/**
 * For each phosphorylation file, add all the modifications and respective flanking
 * amino acid information into the mod table.
 */
void add_phospho_files_to_table(MOMO_OPTIONS_T* options,
                                        SUMMARY_T* summary,
                                        SEQ_T** all_sequences,
                                        int num_sequences) {
  int i;
  for (i = 0; i < arraylst_size(options->phospho_filenames); ++i) {
      add_phospho_file_to_table((char*) arraylst_get(i, options->phospho_filenames), options, summary, all_sequences, num_sequences);
  }
}

/**
 * Given a list of sequences, set the background sequences for all modifications. For example,
 * if a mod is associated with amino acids S,T,Y, then find all the kmers centralized around
 * these letters and add these to the background sequences.
 */
void add_background_sequences_to_table(MOMO_OPTIONS_T* options, SUMMARY_T* summary, SEQ_T** protein_database_sequences, int num_sequences) {
  if (options->algorithm == motifx || options->algorithm == modl) {
    int i;
    int j;
    int k;
    
    ARRAYLST_T* mod_table_keys = summary->mod_table_keys;
    const char* alph_letters = summary->alph_letters;
    
    // for each mod
    for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
      // initialize the background table and list
      MOD_INFO_T* modinfo = (MOD_INFO_T*) hash_get_entry_value(arraylst_get(i, mod_table_keys));
      ARRAYLST_T* bg_seq_list = arraylst_create();
      HASH_TABLE bg_seq_table = (options->eliminate_repeats) ? hash_create(DATA_HASH_SIZE) : NULL;
      modinfo->bg_seq_list = bg_seq_list;
      modinfo->bg_seq_table = bg_seq_table;
      
      // amino acids for this mod
      BOOLEAN_T* amino_acids = modinfo->amino_acids;
      
      // for each sequence
      for (j = 0; j < num_sequences; ++j) {
        char* protein = get_raw_sequence(protein_database_sequences[j]);
        
        // if fasta, we need to search the entire protein
        for (k = 0; k < strlen(protein); ++k) {
          if (options->bg_filetype == fasta || (options->bg_filetype == prealigned && k == options->width/2)) {
            char* alph_letters_substring = strchr(alph_letters, protein[k]);
            if (alph_letters_substring) {
              int aa_idx = alph_letters_substring - alph_letters;
              if (amino_acids[aa_idx]) {
                add_kmer_to_sequence_list(&bg_seq_table, &bg_seq_list, protein, k, options);
              }
            }
          }
        }
      }
    }
  }
}
