#include "array-list.h"
#include "momo.h"

void cleanup_modinfo(MOD_INFO_T* modinfo);

/**
 * If a protein database is provided, parse the sequences
 */
void read_protein_database_sequences(
  MOMO_OPTIONS_T* options,
  SUMMARY_T* summary,
  SEQ_T*** all_sequences,
  int* num_sequences);

/**
 * Sets the O(1) lookup table. This hashes from every unique kmer within
 * the protein database to an arraylist of KMER_LOC_T objects.
 */
void create_hash_fasta_preprocess_table(
  MOMO_OPTIONS_T* options,
  SUMMARY_T* summary,
  SEQ_T** all_sequences,
  int num_sequences);

/**
 * For each phosphorylation file, add all the modifications and respective flanking
 * amino acid information into the mod table.
 */
void add_phospho_files_to_table(
  MOMO_OPTIONS_T* options,
  SUMMARY_T* summary,
  SEQ_T** all_sequences,
  int num_sequences);

/**
 * Given a list of sequences, we want to set the background sequences for all modifications. For example,
 * if a mod is associated with amino acids S,T,Y, then we want to find all the kmers centralized around
 * these letters add these to the background sequences.
 */
void add_background_sequences_to_table(
  MOMO_OPTIONS_T* options,
  SUMMARY_T* summary,
  SEQ_T** protein_database_sequences,
  int num_sequences);


// Structure for hash_fasta. seqidx represents the index of the the sequence within the protein database. kmeridx
// represents the index of the kmer within the sequence.
typedef struct kmerloc {
  int seqidx;                                     // Index of within FASTA file
  int kmeridx;                                    // Index of the kmer within sequence
} KMER_LOC_T;
