/***********************************************************************
 * FILE: motif-spec.h
 * AUTHOR: James Johnson
 * CREATE DATE: 1 September 2011
 * PROJECT: All
 * COPYRIGHT: 2011, UQ
 * DESCRIPTION: Specification of data structure for representing one motif.
 * As per the principal of information hiding this should only be included
 * by files that need to create motifs. Other files should rely on the
 * accessor functions in motif.h
 ***********************************************************************/
#ifndef MOTIF_SPEC_H
#define MOTIF_SPEC_H

struct motif_t {
  int idx; // The motif loading index.
  char id[MAX_MOTIF_ID_LENGTH+2]; // The motif ID. First character is strand [+-?]
  char id2[MAX_MOTIF_ID_LENGTH+1]; // The motif secondary ID (for Tomtom).
  char *consensus;    // Single letter consensus computed from PSSM.
  int length;         // The length of the motif.
  ALPH_T* alph;       // The alphabet of the motif.
  int flags;          // A bitfield of flags
  double evalue;      // E-value associated with this motif.
  double log_evalue;  // log (base 10) of E-value associated with this motif.
  double num_sites;   // Number of sites associated with this motif.
  double complexity;  // Complexity of this motif.
  MATRIX_T* freqs;    // The frequency matrix. [length by (alph_size + ambigs)] 
  MATRIX_T* scores;   // Position specific scoring matrix [length by alph_size]
  char *url;          // Optional url, empty string if not set.
  int trim_left;      // count of columns trimmed on left, 0 = untrimmed
  int trim_right;     // count of columns trimmed on right, 0 = untrimmed
};

#define MOTIF_HAS_AMBIGS 1
#define MOTIF_BOTH_STRANDS 2
#define MOTIF_RESERVED_FLAGS 2

/***********************************************************************
 * Allocate memory for a MEME motif and initialize from score and
 * frequency matrices.
 ***********************************************************************/
MOTIF_T* allocate_motif
  (char *id, char *id2, ALPH_T *alph, MATRIX_T* freqs, MATRIX_T* scores);

/***********************************************************************
 * Set the motif id.
 ***********************************************************************/

void set_motif_id
  (const char* id,
   int len,
   MOTIF_T* motif);

void set_motif_id2
  (const char* id,
   int len,
   MOTIF_T* motif);

/***********************************************************************
 * Set the strand of a motif.
 ***********************************************************************/
void set_motif_strand
  (char strand,
   MOTIF_T *motif);

/***********************************************************************
 * Set the site count of a motif
 ***********************************************************************/
void set_motif_nsites(MOTIF_T* motif, int num_sites);

/***********************************************************************
 * Set the url of a motif
 ***********************************************************************/
void set_motif_url
  (char    *url,
   MOTIF_T *motif);

#endif
