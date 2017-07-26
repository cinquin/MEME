#ifndef SPAMO_MATCHES_H
#define SPAMO_MATCHES_H

#include "red-black-tree.h"
#include "linked-list.h"
#include "motif.h"

// constants that can be ORd together to create the 4 quadrants
#define LEFT 0
#define RIGHT 1
#define SAME 0
#define OPPO 2

// constants for the orientations that combine quadrants
#define UP_SEC_PAL 4
#define UP_PRI_PAL 5
#define DOWN_PRI_PAL 6
#define DOWN_SEC_PAL 7
#define BOTH_PAL 8

// constants used for constructing an quadrant set from an orient
#define Q_TOP_LEFT 1
#define Q_TOP_RIGHT 2
#define Q_BOTTOM_LEFT 4
#define Q_BOTTOM_RIGHT 8
extern const uint8_t orient2set[];

#define QUADRANT_CONTRIBUTES(quad, orient) ((orient2set[orient] & (1 << (quad)))  != 0)

// These macros work for quads or orients. Note that a orient can be both left and right!
#define LEFT_SIDE(orient) ((orient2set[orient] & (Q_TOP_LEFT | Q_BOTTOM_LEFT)) != 0)
#define RIGHT_SIDE(orient) ((orient2set[orient] & (Q_TOP_RIGHT | Q_BOTTOM_RIGHT)) != 0)

// These macros only work with quads, do not use them for orients
#define STRAND(quad) ((quad) & OPPO)
#define SIDE(quad) ((quad) & RIGHT)

#define NQUADS 4        // four basic orientations
#define NORIENTS 9      // 9 biologically interesting combinations of 4 quadrants
#define BLOCKSIZE 100   // for allocating chunks of space for sequence indexes

/*
 * To conserve space the different values have the following meanings:
 * 0            : there is no match
 * value > 0    : there is a match on the positive strand at value
 * value < 0    : there is a match on the negative strand at abs(value)
 */

typedef struct {
  char *source;
  char *name;
  time_t last_mod;
  int loaded;
  int excluded_duplicate_id;
  int excluded_tooshort;
  int excluded_nomatch;
  int excluded_ambigs;
  int excluded_similar;
  int erased_primary_matches;
} SEQUENCE_DB_T;

typedef struct {
  int idx; // output index (only given if contributes to a significant spacing)
  int index; // internal index
  int length;
  char *data;
  char *name;
  ARRAY_T *primary_matches; // Position [0] contains location of primary match.
  ARRAY_T *trimmed_primary_matches;  // Re-scan of trimmed sequences for erasing extra primary sites.
  bool contributes; // contributes to a significant spacing
} SEQUENCE_T;

typedef struct {
  int id; //an id assigned to the database based on the order it was listed
  char *source;
  char *name;
  time_t last_mod;
  int loaded;
  int excluded;
} MOTIF_DB_T;

typedef struct {
  int bin_count;        // number of possible bins
  int *count;           // count[i] is number of occurences in bin i
  double *pvalue;       // pvalue[i] is (adjusted) p-value for bin i
  LINKLST_T **sequences;// sequences[i] is the list of sequences contributing to count[i]
} SPACING_T;

typedef struct {
  double pvalue;        // adjusted p-value of significant spacing
  int bin;              // bin of spacing
  int orient;           // orientation of motifs
  MOTIF_T *alignment_motif[NQUADS];     // alignment motif of contributing sequences in each quadrant
  MOTIF_T *inferred_motif;      // inferred secondary motif based on alignment
} SIGSPACE_T;

typedef struct {
  int idx; // an index assigned based on the order the motifs were written out to file
  MOTIF_DB_T *db;
  MOTIF_T *motif;
  SPACING_T spacings[NORIENTS];
  int total_spacings;
  int max_in_one_bin;
  SIGSPACE_T *sigs; //significant spacings
  int sig_count;
  double min_pvalue; //minimum p-value of spacings
  int best_orient; // keep track of the best orient (in case no spacings pass threshold)
  int best_bin; // keep track of the best bin (in case no spacings pass threshold)
  BOOLEAN passes_evalue_cutoff;
  int *seqs;       // the ids of the best sequences
  int seq_count;
} SECONDARY_MOTIF_T;

typedef struct {
  int db_id;
  char *motif_id;
} SECONDARY_KEY_T;

typedef struct GROUPED_MOTIF {
  SECONDARY_MOTIF_T *best;
  LINKLST_T *others;
} GROUPED_MOTIF_T;

/**************************************************************************
 * Compares two keys of a secondary motif for equality.
 * First compares based on the motif database and second compares
 * based on the name of the motif.
 **************************************************************************/
int secondary_key_compare(const void *p1, const void *p2);

/**************************************************************************
 * Copies a key of a secondary motif
 **************************************************************************/
void* secondary_key_copy(void *p);

/**************************************************************************
 * Creates a structure to hold the motif database information.
 * All strings are copied.
 **************************************************************************/
SEQUENCE_DB_T* create_sequence_db(char *file);

/**************************************************************************
 * Destroy the sequence db
 **************************************************************************/
void destroy_sequence_db(SEQUENCE_DB_T *db);

/**************************************************************************
 * Creates a structure to hold the secondary database information.
 **************************************************************************/
MOTIF_DB_T* create_motif_db(int id, char *file);

/**************************************************************************
 * Destroy the secondary db
 **************************************************************************/
void destroy_motif_db(void *secondary_db);

/**************************************************************************
 * Destroy the sequence (including the name attribute)
 **************************************************************************/
void destroy_sequence(void *sequence);

/**************************************************************************
 * Create a secondary motif. As the number of sequences is unknown at this
 * point the sequence_matches array is left unallocated. All pvalues are
 * initilized to 1.
 **************************************************************************/
SECONDARY_MOTIF_T* create_secondary_motif(int margin, int bin, 
    MOTIF_DB_T *db, MOTIF_T *motif);

/**************************************************************************
 * Destroys a secondary motif. 
 * It takes a void * pointer so it can be used in the collection objects.
 **************************************************************************/
void destroy_secondary_motif(void *);

/**************************************************************************
 * Create a group with a good secondary motif and a group
 * of redundant secondary motifs.
 **************************************************************************/
GROUPED_MOTIF_T* create_grouped_motif(SECONDARY_MOTIF_T* best);

/**************************************************************************
 * Destroys a grouped motif
 * 
 * Relies on the secondary motifs begin destroyed elsewhere.
 **************************************************************************/
void destroy_grouped_motif(void *p);

/**************************************************************************
 * Calculate the total number of pvalue calculations that will be done
 * by the program. This number is used to correct the pvalues for multiple
 * tests using a bonferoni correction.
 **************************************************************************/
int calculate_test_count(
  int margin, 
  int bin, 
  int test_max, 
  RBTREE_T *secondary_motifs
);

/**************************************************************************
 * This does most of the calculation steps after the matching
 * positions of the motifs have been derived from a scan
 * file. The steps undertaken are:
 * 1) bin the matches
 * 2) compute the pvalues of the tested region around the primary
 * 3) compute the sequence id set of the most significant peak
 * 4) set the state of the motif to loaded.
 *
 **************************************************************************/
void process_matches(
  int margin, 
  int bin, 
  double significance_threshold, 
  double motif_evalue_cutoff,
  int test_max,
  MOTIF_T *primary_motif, 
  RBTREE_T *sequences,
  SECONDARY_MOTIF_T *secondary_motif,
  int n_secondary_motifs,
  //int *secondary_matches
  ARRAY_T **secondary_matches
);

/**************************************************************************
 * Copy a (segment of a) string, optionally with reverse complementation.
 *
 * The destination string must be allocated and freed by the caller
 * and must be at least size "length+1".
 * Copying begins at the address "from_str" and the segment of
 * length bytes is copied to "to_str", with reverse complementation
 * if requested.  Then a null is added.
 **************************************************************************/
void copy_string_with_rc(
  ALPH_T *alph,         // Alphabet
  char *from_str,       // Start of segment to copy.
  char *to_str,         // Destination for copy.
  int length,           // The number of characters to copy.
  BOOLEAN_T rc          // Copy the reverse complement of the segment.
);

#endif
