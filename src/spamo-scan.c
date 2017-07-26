
#include "matrix.h"
#include "alphabet.h"
#include "assert.h"
#include "pssm.h"
#include "spamo-matches.h"
#include "spamo-scan.h"
#include "binary-search.h"

#include <stdlib.h>

/**************************************************************************
 * Scans specified regions of a sequence with a PSSM and its reverse 
 * complement and returns a list of matches (positions) that score above
 * the specified threshold.  If best_hit_only is TRUE,
 * it returns the  single best match or randomly selects among ties.
 * The strand of each match is encoded in its position: negative
 * means opposite strand.
 **************************************************************************/
static inline ARRAY_T* best_in_region(
  PSSM_T *pssm,             // the motif score matrix
  PSSM_T *pssm_rc,          // the revese complement motif score matrix
  char *seq,                // the sequence to score
  int *regions,             // the regions to scan
  double score_threshold,   // the minimum score considered a hit
  int *hits,                // the pre-allocated hits cache should be large enough to store a hit for every position
  int hits_size,            // the size of the pre-allocatted hits cache.
  BOOLEAN_T best_hit_only   // old SpaMo: used best hit only
) {
  int i, region_index, position_index, motif_offset, motif_length, region_end;
  int aindex, success, hit_count;
  char c;
  double lo, lo_rc; 
  double *hits_lo;
  ATYPE best_lo;
  ALPH_T *alph;
  MATRIX_T *pssm_matrix, *pssm_rc_matrix;
  bool revcomp;
  //basic checks on parameters
  assert(pssm != NULL);
  assert(pssm_rc == NULL || get_pssm_w(pssm) == get_pssm_w(pssm_rc));
  assert(seq != NULL);
  assert(regions != NULL);
  // check if we're scanning revcomp
  revcomp = (pssm_rc != NULL);
  // init hits list
  hits_lo = best_hit_only ? (double *)mm_malloc(hits_size * sizeof(double)) : NULL;
  //get the alphabet
  alph = pssm->alph;
  //initilize the best position: 0 means no match; +ve is match on 5` strand; -ve is match on 3` strand
  // to make this work an origin of 1 is used.
  best_lo = -BIG;
  hit_count = 0;
  // get the motif length
  motif_length = get_pssm_w(pssm);
  // get the pssm matrixes
  pssm_matrix = pssm->matrix;
  pssm_rc_matrix = (revcomp ? pssm_rc->matrix : NULL);
  //scan regions, check for overlaps and those not in sorted order
  for (region_index = 0; regions[region_index] != -1; region_index += 2) {
    //error checks, so nothing is too suprising
    assert(regions[region_index] >= 0);
    assert(region_index == 0 || regions[region_index-1] < regions[region_index]);
    //scan a region
    region_end = regions[region_index+1] - motif_length + 1;
    for (position_index = regions[region_index]; position_index <= region_end; ++position_index) {
      //scan a position
      success = TRUE; // assume success
      lo = 0; //reset sums to zero
      lo_rc = 0;
      for (motif_offset = 0; motif_offset < motif_length; ++motif_offset) {
        c = seq[position_index+motif_offset];

        // Check for gaps and ambiguity codes at this site
        if(!alph_is_core(alph, c)) {
          success = FALSE;
          break;
        }
        //note these scores are scaled
        aindex = alph_index(alph, c);
        lo += get_matrix_cell(motif_offset, aindex, pssm_matrix);
        if (revcomp) lo_rc += get_matrix_cell(motif_offset, aindex, pssm_rc_matrix); 
      }
      if (success) {
        // now revert the scaled lo scores to unscaled
        lo = get_unscaled_pssm_score(lo, pssm);
        if (revcomp) lo_rc = get_unscaled_pssm_score(lo_rc, pssm);

        // Save hit if it meets threshold.
        double hit_lo;
        int position;
        if (lo >= score_threshold || (revcomp && lo_rc >= score_threshold)) {
          // Determine score and the strand; break ties randomly; negative means reverse strand.
          if (revcomp && (lo_rc > lo || (lo_rc == lo && ( ((double)rand()/RAND_MAX) > 0.5) )))  {
            hit_lo = lo_rc;
            position = -(position_index + 1);   // negative strand
          } else {
            hit_lo = lo;
            position = position_index + 1;      // positive strand
          }
          // Save all hits, keeping track of index of best one.
          hits[hit_count] = position;         // Save hit.
          if (best_hit_only) hits_lo[hit_count] = hit_lo;     // Save hit logodds if necessary.
          if (hit_lo > best_lo) {
            best_lo = hit_lo;
          }
          hit_count++;
        } // score above threshold
      } // found hit (no ambigs)
    } // position
  } // region

  // Return an array of all the hits; put a (randomly chosen) best hit first in array.
  if (best_lo >= score_threshold) {
    ARRAY_T* hits_array = allocate_array(hit_count);
    for (i=0; i<hit_count; i++) set_array_item(i, (ATYPE)(hits[i]), hits_array); 
    // If best_hit_only, randomly choose among ties for best hit.
    if (best_hit_only) {        // put best hits first
      int nbest = 0;
      int tail = hit_count-1;
      // Put ties in front of hits_array
      for (i=0; i<hit_count; i++) {
        int index = (hits_lo[i] == best_lo) ? nbest++ : tail--;
        set_array_item(index, (ATYPE)(hits[i]), hits_array);
      }
      // Now randomly choose a best hit and move to front of hits_array if there is more than one.
      if (nbest > 1) {
        int best_index = nbest * ((double)rand()/RAND_MAX);
        ATYPE tmp_best_hit = get_array_item(best_index, hits_array);
        ATYPE tmp_first_hit = get_array_item(0, hits_array);
        set_array_item(0, tmp_best_hit, hits_array);
        set_array_item(best_index, tmp_first_hit, hits_array);
      }
    } else {    // otherwise, order of hits doesn't matter
      for (i=0; i<hit_count; i++) set_array_item(i, (ATYPE)(hits[i]), hits_array); 
    }
    myfree(hits_lo);
    return hits_array;
  } else {
    return NULL;
  }

} // best_in_region

/**************************************************************************
* Compare two ATYPEs and return <0, 0, >0
* if the first value is <, =, > the second value.
**************************************************************************/
int atype_compare(
  const void *v1,
  const void *v2
)
{
  double result;

  const ATYPE* s1 = (const ATYPE*) v1;
  const ATYPE* s2 = (const ATYPE*) v2;

  if ((result = s2 - s1) != 0) {
    return (result<0) ? -1 : +1;
  } else {
    return (int) (s2 - s1);
  }
}

/**************************************************************************
 * Scan all the sequences with the primary motif and store the best match
 * for each in the sequence data structure.
 **************************************************************************/
void scan_spamo_primary(
  int margin, 
  double score_threshold, 
  ARRAY_T *background, 
  MOTIF_T *motif, 
  RBTREE_T *sequences,
  BOOLEAN_T trimmed
) {
  MOTIF_T *motif_rc;
  PSSM_T *pssm, *pssm_rc;
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int regions[3];
  int *hits, hits_size;
  bool revcomp;
  if ((node = rbtree_first(sequences)) == NULL) return;
  revcomp = alph_has_complement(get_motif_alph(motif));
  sequence = (SEQUENCE_T*)rbtree_value(node);
  hits_size = (sequence->length - 2 * margin) * 2;
  hits = mm_malloc(sizeof(int) * hits_size);

  regions[0] = margin; // first valid index to scan
  regions[2] = -1; //terminate list with negative 1
  //prepare a reverse complement of the motif
  motif_rc = (revcomp ? dup_rc_motif(motif) : NULL);
  //convert the motif and the reverse complement into PSSMs
  pssm =    build_motif_pssm(motif,    background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  pssm_rc = (revcomp ? build_motif_pssm(motif_rc, background, background, NULL, 0, PSSM_RANGE, 0, FALSE) : NULL);
  // if the score threshold is between 0 and -1 then compute a threshold from the pssm
  if (score_threshold < 0 && score_threshold >= -1) {
    score_threshold = pssm_best_match_score(pssm) * (-score_threshold);
  }
  //scan each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    if (((sequence->length - 2 * margin) * 2) > hits_size) {
      hits_size = (sequence->length - 2 * margin) * 2;
      hits = mm_realloc(hits, sizeof(int) * hits_size);
    }
    // Scan only the margins to find the best hit there.
    regions[1] = sequence->length - margin - 1; // last valid index to include in scan
    if (trimmed) {
      sequence->trimmed_primary_matches = best_in_region(pssm, pssm_rc, sequence->data, regions, score_threshold, hits, hits_size, FALSE);
    } else {
      sequence->primary_matches = best_in_region(pssm, pssm_rc, sequence->data, regions, score_threshold, hits, hits_size, TRUE);
    }
  }
  //clean up the PSSMs
  free_pssm(pssm);
  free_pssm(pssm_rc);
  //clean up the motif
  destroy_motif(motif_rc);
  //clean up the hits
  free(hits);
}

/**************************************************************************
 * Scan all the sequences with the secondary motif and store the best match
 * for each in matches list.
 **************************************************************************/
void scan_spamo_secondary(
  int margin, 
  double score_threshold,
  BOOLEAN_T use_best_secondary,
  ARRAY_T *background, 
  MOTIF_T *motif,
  RBTREE_T *sequences, 
  ARRAY_T **matches,
  int *hits,
  int hits_size
) {
  MOTIF_T *motif_rc;
  PSSM_T *pssm, *pssm_rc;
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int regions[5];
  bool revcomp;
  revcomp = alph_has_complement(get_motif_alph(motif));
  if ((node = rbtree_first(sequences)) == NULL) return;
  sequence = (SEQUENCE_T*)rbtree_value(node);
  // set up the two regions flanking the primary motif
  regions[0] = 0; // first valid index to scan
  regions[1] = margin - 1;
  regions[2] = sequence->length - margin;
  regions[3] = sequence->length - 1;
  regions[4] = -1; // terminate list with negative 1
  // prepare a reverse complement of the motif
  motif_rc = (revcomp ? dup_rc_motif(motif) : NULL);
  // convert the motif and the reverse complement into PSSMs
  pssm =    build_motif_pssm(motif,    background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  pssm_rc = (revcomp ? build_motif_pssm(motif_rc, background, background, NULL, 0, PSSM_RANGE, 0, FALSE) : NULL);
  // if the score threshold is between 0 and -1 then compute a threshold from the pssm
  if (score_threshold < 0 && score_threshold >= -1) {
    score_threshold = pssm_best_match_score(pssm) * (-score_threshold);
  }
  // scan all the sequences
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    matches[sequence->index] = best_in_region(pssm, pssm_rc, sequence->data, regions, score_threshold, hits, hits_size, use_best_secondary);
  }
  //clean up the PSSMs
  free_pssm(pssm);
  free_pssm(pssm_rc);
  //clean up the motif
  destroy_motif(motif_rc);
}
