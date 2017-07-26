/*
 * $Id: read_seq_file.h 776 2006-05-10 17:27:13Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 19:09:14  nadya
 * Initial revision
 *
 */

#ifndef READ_DATA_FILE_H
#define READ_DATA_FILE_H

extern DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  ALPH_T *alph,			/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  BOOLEAN ignore_dup 		// allow duplicate names
);

extern DATASET *create_meme_dataset_from_momo(
  ARRAYLST_T* seq_array,		/* name of file to open */
  ALPH_T *alph,		/* alphabet used in sequences */
  int width,
  BOOLEAN_T eliminate_repeats
);

extern SAMPLE *get_sample_by_name(
  char *sample_name
);
extern void add_control_samples (
  DATASET *dataset,                     // primary dataset
  DATASET *control,                     // control dataset
  int c1,                               // size of control group 1
  int c2,                               // size of control group 2
  BOOLEAN use_comp                      // consider both strands
);
extern void shuffle_dataset_order(
  DATASET *dataset              // dataset to shuffle
);
extern void shuffle_dataset_letters(
  DATASET *dataset,	// dataset to have sequences shuffled and replaced
  int kmer,		// preserve k-mer counts in each sequence
  int revcomp		// include reverse complements of sequences
);
extern void destroy_sample(SAMPLE *sample);
#endif
