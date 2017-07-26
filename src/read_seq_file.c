/***********************************************************************
*								       *
*	MEME							       *
*	Copyright 1994-2015, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 9-30-99 tlb; add sequence to SAMPLE */
/* 6-29-99 tlb; add resic */
/* 9-22-08 pm; add posprior */
/*
	Supports FASTA and SWISSPROT formats.

Note:
	Assumes SWISSPROT format if "ID" appears at beginning of first line.
	Ignores duplicate (same <id>) sequences.

All formats:
	<id>		:= <alpha>+
			   The *unique* name of the sequence.  May *not*
			   contain whitespace!

			   NOTE: if <id> is "WEIGHTS", the <description>
			   field is parsed as a string of sequence weights 
			   and any <sequence> field (if any) is ignored.
			   All weights are assigned in order to the 
			   sequences in the file. If there are more sequences
			   than weights, the remainder are given weight one.
			   Weights must be greater than zero and less than
			   or equal to one.  Weights may be specified by
			   more than one "WEIGHTS" entry which may appear
			   anywhere in the file.

	<description>	:= <alpha>+
			   A description of the sequence;
		           may contain whitespace.

	<sequence>	:= <alpha>+
			   The DNA or protein sequence data; 
                           may contain whitespace.

	<text>		:= <alpha>*
			   Text which is ignored.


Pearson FASTA format:
	INPUT DATA FILE FORMAT:
		[
		<header>
		<sequence>+
		]*

	<header>	:= ">"<id>" "<description>
			    Everything up to the first space is 
			    assumed to be the unique name of the sequence.
	
SWISSPROT Format
	[
	<ID>
	<DE>+
	<SQ>
	<sequence>+
	//
	]*

	<ID>		:= "ID   "<id>" "<text>
	<DE>		:= "DE   "<description>
	<SQ>		:= "SQ   "<text>
*/

#include "meme.h"
#include "read_sequence.h"
#include "ushuffle.h"
#include "seq.h"

static void parse_sequence_weights(
  DATASET *dataset,
  char *sample_name,
  char *sample_de,
  char *sequence
);

void set_database_freq_and_entropy(
  DATASET *dataset,
  BOOLEAN use_comp  		/* use complementary strands, too */
);

static inline int compare_res(
  const uint8_t *res1,
  int len1,
  const uint8_t *res2,
  int len2
);

extern inline int compare_samples_by_res (
  const void *elem1,
  const void *elem2
);

static void free_sample(
  SAMPLE *s			// sample to free
);

static SAMPLE *create_sample(
  ALPH_T *alph,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence,  	/* the sequence */
  char *descript,	// Description of the sample
  BOOLEAN revcomp	// make space for reverse complement
);

/* local variables */
// JIJ: this was 100 000 but I changed it to 100 003 to make it prime because
// apparently prime numbers are better for avoiding hash collisions.
#define DATA_HASH_SIZE 100003
// The hash table will be appended to if read_seq_file is called more
// than once.  This is so that positive sequences the negative sequence file
// will be ignored.
static HASH_TABLE ht_seq_names = NULL;	/* hash of dataset seq names */

/**********************************************************************/
/*
        get_sample_by_name

        Check if sample name is defined in ht_seq_names hash table
*/
/**********************************************************************/
SAMPLE *get_sample_by_name(
  char *sample_name
)
{
  HASH_TABLE_ENTRY * hash_entry = hash_lookup_str(sample_name, ht_seq_names);
  return(hash_entry != NULL ? (SAMPLE *) hash_get_entry_value(hash_entry) : NULL);
} // get_sample_by_name


DATASET *create_meme_dataset_from_momo(
                                       ARRAYLST_T* seq_array,		/* name of file to open */
                                       ALPH_T *alph,		/* alphabet used in sequences */
                                       int width,
                                       BOOLEAN_T eliminate_repeats) {
  int i, j;
  char *sample_name;                      /* name of sample read */
  char *sample_de;                        /* descriptor text for sample */
  char *sequence;                         /* sequence read */
  long length = (long) width;			/* length of sequence */
  BOOLEAN error=FALSE;		/* none yet */
  
  /* create a dataset */
  DATASET *dataset = (DATASET *) mymalloc(sizeof(DATASET));
  dataset->alph = alph_hold(alph);
  dataset->psp_w = 0;			// indicates no PSP was read
  dataset->log_psp_w = 0;		// so log_psp will get initialized
  dataset->samples = NULL;		/* no samples */
  dataset->n_samples = 0;		/* no samples yet */
  dataset->n_group[0] = dataset->n_group[1] = dataset->n_group[2] = 0;
  dataset->seq_weights = NULL;		// user-proved sequence weights
  dataset->n_wgts = 0;
  /* initialize maximum length of sequences */
  dataset->max_slength = 0;
  dataset->min_slength = 10000000;
  
  printf("here1\n");
  for (i = 0; i < arraylst_size(seq_array); ++i) {
    printf("here\n");
    SEQ_T* seqobject = eliminate_repeats ? hash_get_entry_value(arraylst_get(i, seq_array)) : arraylst_get(i, seq_array);
    sample_name = get_seq_name(seqobject);
    sample_de = get_seq_description(seqobject);
    sequence = get_raw_sequence(seqobject);
    
    //    HASH_TABLE_ENTRY *hash_entry; // needed to add pointer to sample
    
    /* create the sample */
    SAMPLE *sample = create_sample(alph, length, sample_name, sequence, sample_de, FALSE);
    if (sample == NULL) {error = TRUE; continue;}
    
    // Record the original index of the sample so we can
    // access it after the dataset is shuffled.
    sample->orig_index = dataset->n_samples;
    
    // All samples are in Group 0 initially.
    sample->group = 0;
    
    /* record maximum length of actual sequences */
    dataset->max_slength = length;
    dataset->min_slength = length;
    
    printf("min_slength: %lu\n", dataset->min_slength);
    
    /* put the sample in the array of samples */
    if ((dataset->n_samples % RCHUNK) == 0) {
      Resize(dataset->samples, dataset->n_samples + RCHUNK, SAMPLE *);
    }
    dataset->samples[dataset->n_samples++] = sample;
    //    hash_entry = hash_lookup_str(sample_name, ht_seq_names);
    //    if (!hash_entry) {
    //      fprintf(stderr, "hash error: added sample ID then failed to find it in read_seq_file ()\n");
    //      error = TRUE;
    //      break;
    //    }
    //    hash_set_entry_value(sample, hash_entry);
    
    /* cleanup sequence (create_sample copies it) */
    myfree(sequence);
    
  } /* sequences */
  if (length < 0) error = TRUE;			/* read_sequence error */
  
  /* check that datafile contained at least one sample */
  if (!error) {
    if (dataset->n_samples == 0) {
      fprintf(stderr, "No sequences found in files.  Check file format.\n");
      error = TRUE;
    }
  }
  
  /* exit if there was an error */
  if (error) exit(1);
  
  /* resize the array of samples and intialize the size of Group 0*/
  if (dataset->n_samples) Resize(dataset->samples, dataset->n_samples, SAMPLE*);
  dataset->n_group[0] = dataset->n_samples;
  
  dataset->input_order = NULL;
  
  // Set the frequencies and entropies of the database.
  set_database_freq_and_entropy(dataset, FALSE);
  
  return dataset;
}

/**********************************************************************/
/*
	read_seq_file

	Open a sequence file and read in the sequences. 
	Returns a dataset->

	setup_hash_alphabet must have been called prior to first call to this.
*/
/**********************************************************************/
DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  ALPH_T *alph,		/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  BOOLEAN ignore_dup 		// allow duplicate names
)
{
  int i, j;
  char *sample_name;		/* name of sample read */
  char *sample_de;		/* descriptor text for sample */
  char *sequence;		/* sequence read */
  long length;			/* length of sequence */
  BOOLEAN error=FALSE;		/* none yet */

  /* create a hash table of sequence names */
  if (!ht_seq_names) ht_seq_names = hash_create(DATA_HASH_SIZE);

  /* create a dataset */
  DATASET *dataset = (DATASET *) mymalloc(sizeof(DATASET));
  dataset->alph = alph_hold(alph);
  dataset->psp_w = 0;			// indicates no PSP was read
  dataset->log_psp_w = 0;		// so log_psp will get initialized
  dataset->samples = NULL;		/* no samples */
  dataset->n_samples = 0;		/* no samples yet */
  dataset->n_group[0] = dataset->n_group[1] = dataset->n_group[2] = 0;
  dataset->seq_weights = NULL;		// user-proved sequence weights
  dataset->n_wgts = 0;
  /* initialize maximum length of sequences */
  dataset->max_slength = 0;
  dataset->min_slength = 10000000;

  /* open data file */
  FILE *data_file;			/* file pointer */
  if (file_name == NULL) {
    fprintf(stderr, "You must specify a data file or `stdin'.\n");
    exit(1);
  } else if (strcmp(file_name, "stdin")) {
    data_file = fopen(file_name, "r"); 
    if (data_file == NULL) {
      fprintf(stderr, "Cannot open file `%s'.\n", file_name);
      exit(1);
    }
  } else {
    data_file = stdin;
  }

  while (read_sequence(alph, data_file, &sample_name, &sample_de, &sequence, 
    &length)) {

    /* skip sequence if an error occurred */
    if (length < 0) continue;

    /* parse weights if given; make (more than enough) room in array */
    if (strcmp(sample_name, "WEIGHTS")==0) {
      parse_sequence_weights(dataset, sample_name, sample_de, sequence);
      continue;
    }

    /* ignore duplicate (same sample name) sequences */
    if (hash_lookup_str(sample_name, ht_seq_names) != NULL) {
      if (! ignore_dup) {
	fprintf(stderr, "Duplicate sequence name. Skipping '%s'.\n", sample_name);
	myfree(sample_name);
	myfree(sample_de);
	myfree(sequence);
	continue;
      }
    } else {
      // Don't free sample_name since the hash table doesn't copy it.
      hash_insert_str(sample_name, ht_seq_names);  /* put name in hash table */
    }

    HASH_TABLE_ENTRY *hash_entry; // needed to add pointer to sample

    /* create the sample */
    SAMPLE *sample = create_sample(alph, length, sample_name, sequence, sample_de, use_comp);
    if (sample == NULL) {error = TRUE; continue;}

    // Record the original index of the sample so we can
    // access it after the dataset is shuffled.
    sample->orig_index = dataset->n_samples;

    // All samples are in Group 0 initially.
    sample->group = 0;

    /* record maximum length of actual sequences */
    dataset->max_slength = MAX(sample->length, dataset->max_slength);
    dataset->min_slength = MIN(sample->length, dataset->min_slength);

    /* put the sample in the array of samples */
    if ((dataset->n_samples % RCHUNK) == 0) {
      Resize(dataset->samples, dataset->n_samples + RCHUNK, SAMPLE *);
    }
    dataset->samples[dataset->n_samples++] = sample;
    hash_entry = hash_lookup_str(sample_name, ht_seq_names);
    if (!hash_entry) {
      fprintf(stderr, "hash error: added sample ID then failed to find it in read_seq_file ()\n");
      error = TRUE;
      break;
    }
    hash_set_entry_value(sample, hash_entry);

    /* cleanup sequence (create_sample copies it) */
    myfree(sequence);
    
  } /* sequences */
  if (length < 0) error = TRUE;			/* read_sequence error */
  
  /* check that datafile contained at least one sample */
  if (!error) {
    if (dataset->n_samples == 0) {
      fprintf(stderr, "No sequences found in file `%s'.  Check file format.\n",
	      file_name);
      error = TRUE;
    }
  }

  /* exit if there was an error */
  if (error) exit(1);

  /* resize the array of samples and intialize the size of Group 0*/
  if (dataset->n_samples) Resize(dataset->samples, dataset->n_samples, SAMPLE*);
  dataset->n_group[0] = dataset->n_samples;

  dataset->input_order = NULL;

  // Set the frequencies and entropies of the database.
  set_database_freq_and_entropy(dataset, use_comp);

  return dataset;
} /* read_seq_file */

/**********************************************************************/
/*
	shuffle_dataset_order

	Shuffle the order of the sequences in a dataset after
	first ordering them alphabetically by integer-encoding.
	(This makes sure the shuffled order
	is independent of strand and input order.)

*/
/**********************************************************************/
void shuffle_dataset_order(
  DATASET *dataset		// dataset to shuffle
)
{
  int i, j;
  int n = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  // Save an array containing pointers to the samples in their
  // input order.  
  Resize(dataset->input_order, n, SAMPLE*);
  memcpy(dataset->input_order, samples, n*sizeof(SAMPLE*));

  // With shuffling turned on, Classic sometimes gives different results
  // than <= v4.10.2 due to its previous (unintended) dependency on sequence order.
  if (1) {
    // Sort the samples alphabetically by their integer-encoding,
    // accounting for the possibility of two strands.
    qsort(samples, dataset->n_samples, sizeof(SAMPLE*), compare_samples_by_res);

    // Shuffle the order of the sorted SAMPLE pointers using Fisher-Yates.
    for (i=n-1; i>0; i--) {
      j = drand_mt() * (i + 1);
      SAMPLE *tmp = samples[i]; samples[i] = samples[j]; samples[j] = tmp; // swap 
    } 
  } else {
    fprintf(stderr, "Shuffling of sequence order is is SUPPRESSED.  FIXME!!!\n");
  }

} // shuffle_dataset_order

/**********************************************************************/
/*
	add_control_samples

	Add ctfrac fraction of the provided control samples 
	to groups 1 and 2, split evenly.

	Do not free the control dataset after calling this.
*/
/**********************************************************************/
extern void add_control_samples (
  DATASET *dataset,			// primary dataset
  DATASET *control,			// control dataset
  int c1,				// size of control group 1
  int c2,				// size of control group 2
  BOOLEAN use_comp			// consider both strands
)
{
  int i;

  // Resize the samples array to make room.
  int n_control = c1 + c2;
  Resize(dataset->samples, dataset->n_samples + n_control, SAMPLE*);

  // Add each control sample to group 1 or 2.
  for (i=0; i<n_control; i++) {
    // Copy the sample
    SAMPLE *s = control->samples[i];
    // FIXME: Because create_sample doesn't initialize log_cum_back
    // we need to do a direct link.  That means never free the control dataset.
    //SAMPLE *sample = create_sample(dataset->alphabet, s->length, 
    //  s->sample_name, s->seq, s->descript, use_comp);
    SAMPLE *sample = s;
    
    if (i < c1) {
      sample->group = 1;
      dataset->n_group[1]++;
    } else {
      sample->group = 2;
      dataset->n_group[2]++;
    }

    /* record maximum length of actual sequences */
    // Not needed: these are only used in status reports.
    //dataset->max_slength = MAX(sample->length, dataset->max_slength);
    //dataset->min_slength = MIN(sample->length, dataset->min_slength);

    /* put the sample in the array of samples */
    dataset->samples[dataset->n_samples++] = sample;
  }

  // Reset the frequencies and entropies of the database.
  myfree(dataset->res_freq);
  set_database_freq_and_entropy(dataset, use_comp);

//FIXME
//printf("ADDED: group0 %d group1 %d group2 %d n %d\n", 
//  dataset->n_group[0], dataset->n_group[1], dataset->n_group[2], dataset->n_samples);
} // add_control_samples

/**********************************************************************/
/*
	parse_sequence_weights

	Read sequence weights from a FASTA header:
		>WEIGHTS w_1 w_2 ... w_n
*/
/**********************************************************************/
static void parse_sequence_weights(
  DATASET *dataset,
  char *sample_name,
  char *sample_de,
  char *sequence
)
{
  double wgt; 
  char *wgt_str = sample_de;
  Resize(dataset->seq_weights, dataset->n_wgts+(int)strlen(wgt_str), double);
  while (sscanf(wgt_str, "%lf", &wgt) == 1) {
    if (wgt <= 0 || wgt > 1) {
      fprintf(stderr, 
	"Weights must be larger than zero and no greater than 1.\n");
      exit(1);
    }
    dataset->seq_weights[dataset->n_wgts++] = wgt;	/* save weight */
    wgt_str += strspn(wgt_str, "      ");		/* skip white */
    wgt_str += strcspn(wgt_str, "     ");		/* skip token */
  }
  myfree(sample_name);
  myfree(sample_de);
  myfree(sequence);
} // parse_sequence_weights

/**********************************************************************/
/*
	set_database_freq_and_entropy

  	Calculate the prior residue frequencies and entropies 
     	and |D|, size of the dataset
  	tlb; 5/9/97 wgt_total_res and weighted res_freq
*/
/**********************************************************************/
void set_database_freq_and_entropy(
  DATASET *dataset,
  BOOLEAN use_comp  		/* use complementary strands, too */
)
{
  int i, j;
  double *seq_weights = dataset->seq_weights;
  int n_wgts = dataset->n_wgts;
  ALPH_T *alph = dataset->alph;	// alphabet

  // count residues
  dataset->res_freq = NULL;
  Resize(dataset->res_freq, alph_size_core(alph), double);
  for (i=0; i < alph_size_core(alph); i++) { dataset->res_freq[i] = 0; }
  dataset->total_res = 0;
  dataset->wgt_total_res = 0;
  for (i=0; i<dataset->n_samples; i++) {		/* sequence */
    long slen = dataset->samples[i]->length;
    double sw = dataset->samples[i]->sw = (n_wgts > i ? seq_weights[i] : 1);
    dataset->total_res += slen;
    dataset->wgt_total_res += slen*sw;
    for (j=0; j < alph_size_core(alph); j++) {
      if (use_comp) { 			/* using complementary strand as well */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j]/2.0;
        dataset->res_freq[alph_complement(alph, j)] +=
          sw * dataset->samples[i]->counts[j]/2.0;
      } else {				/* not using complementary strands */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j];
      }
    }
  }
  /* convert counts to frequencies */
  for (i=0; i < alph_size_core(alph); i++) {
    dataset->res_freq[i] /= dataset->wgt_total_res;
  }
} // set_database_freq_and_entropy

/**********************************************************************/
/*
        shuffle_dataset_letters

        Shuffle the letters in each sequence in a dataset, preserving
        the k-mer counts in each sequence.
        Renames the sequences appending "_shuf_copy#".
        Modifies the dataset in place.
*/
/**********************************************************************/
void shuffle_dataset_letters(
  DATASET *dataset,	// dataset to have sequences shuffled and replaced
  int kmer,		// preserve k-mer counts in each sequence
  int revcomp		// include reverse complements of sequences
)
{
  int i,j,k;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  int n_copies = 1;
  int pad = log(n_copies)/log(10) + 1;
  int i_copy = 1;

  if (kmer > 0) {
    for (i=0; i<n_samples; i++){		/* sequence */
      SAMPLE *s = samples[i];
      int length = s->length;
      // Replace the sequence with its reverse complement if it is smaller.
      // Need to do this on the integer-encoded sequence due to aliases.
      if (revcomp && compare_res(s->resic, s->length, s->res, s->length) < 0) {
	char *rc = NULL;
	Resize(rc, length+1, char);
        for (j=0, k=length-1; j<length; j++, k--) rc[j] = comp_sym(dataset->alph, s->seq[k]);
        rc[length] = '\0';
        SAMPLE *new = create_sample(dataset->alph, s->length, s->sample_name, rc, s->descript, revcomp);
        myfree(s);				// free larger sample
        s = new;				// replace the sample
      }
      // Create the new sample name.
      char *name = NULL; Resize(name, strlen(s->sample_name) + pad + 14, char);
      sprintf(name, "%s_shuf_%d", s->sample_name, i_copy);
      // Create a shuffled copy of the sequence.
      char *shuf = NULL; Resize(shuf, length+1, char);
      shuf[length] = '\0';
      ushuffle1(s->seq, length, kmer);
      ushuffle2(shuf);
      // Create a new sample.
      SAMPLE *new = create_sample(dataset->alph, length, name, shuf, s->descript, revcomp);
      // Replace sample with a shuffled version.
      samples[i] = new;
      free_sample(s);				// free old sample
    }
  }
} // shuffle_dataset_letters

/**********************************************************************/
/*
	compare_res

	Compare integer-encoded sequences.
*/
/**********************************************************************/
static inline int compare_res(
  const uint8_t *res1,
  int len1,
  const uint8_t *res2,
  int len2
)
{
  int i;
  int len = MIN(len1, len2);
  for (i=0; i<len; i++) {
    if (res1[i] < res2[i]) {
      return(-1);
    } else if (res1[i] > res2[i]) {
      return(+1);
    }
  }
  // res1 == res2
  if (len1 < len2) {
    return(-1);
  } else if (len1 > len2) {
    return(+1);
  } else {
    return(res1 < res2);		// Break ties using address.
  }
} // compare_res

/**********************************************************************/
/*
	compare_samples_by_res

	Compare two samples based on their integer-coded sequences.
*/
/**********************************************************************/
inline int compare_samples_by_res (
  const void *elem1,
  const void *elem2
)
{
  const uint8_t *res1 = (*(SAMPLE **)elem1)->res;
  const uint8_t *resic1 = (*(SAMPLE **)elem1)->resic;
  const int len1 = (*(SAMPLE **)elem1)->length;

  const uint8_t *res2 = (*(SAMPLE **)elem2)->res;
  const uint8_t *resic2 = (*(SAMPLE **)elem2)->resic;
  const int len2 = (*(SAMPLE **)elem2)->length;

  // Replace res with MIN(res, resic) to make 
  // comparison insensitive to strand given.
  if (resic1 && compare_res(resic1, len1, res1, len1) < 0) res1 = resic1;
  if (resic2 && compare_res(resic2, len2, res2, len2) < 0) res2 = resic2;

  // Now compare the integer-encoded sequences.
  return (compare_res(res1, len1, res2, len2));
} // compare_samples_by_res

/**********************************************************************/
/*
	free_sample
*/
/**********************************************************************/
static void free_sample(
  SAMPLE *s			// sample to free
)
{
  BOOLEAN revcomp = (s->resic != NULL);
  myfree(s->sample_name);
  myfree(s->descript);
  myfree(s->seq);
  myfree(s->res);
  myfree(s->resic);
  myfree(s->weights);
  myfree(s->not_o);
  myfree(s->log_not_o);
  myfree(s->pY);
  myfree(s->pYic);
  if (revcomp) {
    s->z -= s->length;
    s->log_psp -= s->length;
  }
  myfree(s->z);
  myfree(s->psp_original);
  myfree(s->log_psp);
  myfree(s->counts);
  myfree(s->logcumback);
  myfree(s->sites);

  myfree(s);
} // free_sample

/**********************************************************************/
/*
	create_sample

	Create a sample.

	Returns the sample or NULL on error.

*/
/**********************************************************************/
static SAMPLE *create_sample(
  ALPH_T *alph,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence,       /* the sequence */
  char *descript,	// Description of the sample
  BOOLEAN revcomp	// make space for reverse complement
)
{
  long i, j; 
  SAMPLE *new1;

  /* disallow zero length samples */
  if (length == 0) {
    fprintf(stderr, "\nZero length sequences not allowed. (Sequence `%s').\n",
      name);
    return NULL;
  }

  /* create the record to hold the sample and its associated data */
  new1 = (SAMPLE *) mymalloc(sizeof(SAMPLE));

  /* assign the name and sequence data */
  new1->sample_name = strdup(name);
  new1->descript = strdup(descript);
  new1->length = length;
  new1->seq = strdup(sequence);
  new1->sw = 1;

  /* set up encoded version of sequence and weights */
  new1->res = (uint8_t *) mymalloc(length * sizeof(char));
  if (revcomp) {
    new1->resic = (uint8_t *) mymalloc(length * sizeof(char));
  } else {
    new1->resic = NULL;
  }
  new1->weights = (double *) mymalloc(length * sizeof(double));
  for (i=0; i<length; i++) { 
    // convert anything that's not a core alphabet symbol into wildcard
    int e = alph_encodec(alph, sequence[i]);
    new1->res[i] = e;
    // FIXME: make resic obsolete
    if (new1->resic) new1->resic[length-i-1] = alph_complement(alph, e);
    new1->weights[i] = 1.0;
  }

  /* set up erasing arrays */
  new1->not_o = (double *) mymalloc(length * sizeof(double));
  new1->log_not_o = (int *) mymalloc(length * sizeof(int));

  /* set up arrays to hold posterior probabilities */
  create_2array(new1->pY, int, 3, length);
  new1->pYic = (char *) mymalloc(length * sizeof(char));

  /* Set up Z array and log_psp arrays */
  if (revcomp) {
    // Offset pointer -lseq so that z[j], j in [-lseq...+lseq]
    new1->z_buf = mm_malloc(sizeof(double) * (2*length+1));
    new1->z = (new1->z_buf)+length;
    new1->log_psp_buf = mm_malloc(sizeof(double) * (2*length+1));
    new1->log_psp = (new1->log_psp_buf)+length;
  } else {
    new1->z_buf = mm_malloc(sizeof(double) * (length+1));
    new1->z = new1->z_buf;
    new1->log_psp_buf = mm_malloc(sizeof(double) * (length+1));
    new1->log_psp = new1->log_psp_buf;
  }
  new1->psp_original = NULL;
  new1->max_log_psp = LOG(0);

  /* set up array to hold character counts and get character counts */
  new1->counts = (double *) calloc(alph_size_core(alph), sizeof(double));
  double f = 1.0 / alph_size_core(alph); // uniform frequency
  for (i=0; i<length; i++) {
    int e = new1->res[i];
    if (e != alph_wild(alph)) { // normal letter
      new1->counts[e]++;
    } else { // wildcard
      for (j = 0; j < alph_size_core(alph); j++) new1->counts[j] += f;
    }
  }

  new1->logcumback = (double *) mymalloc((length+1) * sizeof(double));

  new1->nsites = 0;
  new1->sites = NULL;
  new1->minpv = 1;
  new1->group = 0;

  return new1;
} /* create_sample */

void destroy_sample(SAMPLE *sample) {
  free(sample->sample_name);
  free(sample->seq);
  free(sample->descript);
  free(sample->res);
  free(sample->resic);
  free(sample->pYic);
  free(sample->weights);
  free(sample->not_o);
  free(sample->log_not_o);
  free(sample->logcumback);
  free(sample->sites);
  //free(sample->psp_original); // seems to be freed elsewhere
  free_2array(sample->pY, 3);
  free(sample->counts);
  free(sample->z_buf);
  free(sample->log_psp_buf);
  memset(sample, 0, sizeof(SAMPLE));
  free(sample);
} // destroy_sample
