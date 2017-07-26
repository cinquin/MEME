/********************************************************************
 * FILE: compute-prior-dist.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy Bailey
 * CREATE DATE: 11/03/2010
 * PROJECT: MEME suite
 * COPYRIGHT: 2010 UW
 ********************************************************************/
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pssm.h"
#include "prior-reader-from-psp.h"
#include "array.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/********************************************************************
 * This program reads a MEME PSP file and computes the binned
 * distribution of priors. The distribution is writen to stdout.
 ********************************************************************/
int main(int argc, char *argv[]) {

  char *usage = "compute-prior-dist <num-bins> <psp-file>";

  if (argc != 3) {
    fprintf(stderr, "Usage: %s\n", usage);
    return -1;
  }

  int num_bins = atoi(argv[1]);
  if (num_bins <= 0) {
    fprintf(stderr, "Usage: %s\n", usage);
    return -1;
  }

  const char *filename = argv[2];

  // Read each prior, find max and min of distribution.
  DATA_BLOCK_READER_T *psp_reader = NULL;
  psp_reader = new_prior_reader_from_psp(FALSE /* Don't try to parse genomic coord.*/, filename);
  DATA_BLOCK_T *psp_block = new_prior_block();

  int prior_array_size = 100;
  ARRAY_T *raw_priors = allocate_array(prior_array_size);
  int num_priors = 0;
  while (psp_reader->go_to_next_sequence(psp_reader) != FALSE) {
    while (psp_reader->get_next_block(psp_reader, psp_block) != FALSE) {
      double prior = get_prior_from_data_block(psp_block);
      if (prior == 0.0) {
        // Skip priors that are exactly 0.0
        continue;
      }
      if (num_priors == INT_MAX) {
        die("Number of priors exceeded maximum allowed value of %d", INT_MAX);
      }
      set_array_item(num_priors, prior, raw_priors);
      ++num_priors;
      if (num_priors >= prior_array_size) {
        resize_array(raw_priors, 2 * prior_array_size);
        prior_array_size = 2 * prior_array_size;
      }
    }
  }
  free_data_block(psp_block);
  free_data_block_reader(psp_reader);

  ARRAY_T *priors = extract_subarray(raw_priors, 0, num_priors);
  free_array(raw_priors);
  double median_prior = compute_median(priors);
  double min_prior = get_array_item(0, priors);
  double max_prior = get_array_item(num_priors - 1, priors);

  // Print min, max, and median
  printf("#min %6.5f\n", min_prior);
  printf("#max %6.5f\n", max_prior);
  printf("#median %6.5f\n", median_prior);

  // Special case if priors are exactly uniform.
  if (min_prior == max_prior) {
    printf("%6.5f\n", 1.0);
    return 0;
  }

  // Create the array of bins, intialized to 0.
  double *prior_dist = mm_calloc(num_bins, sizeof(double));
  double scale = (num_bins - 1) / (max_prior - min_prior);
  double offset = min_prior;
  int dist_index = 0;

  int i;
  for (i = 0; i < num_priors; ++i) {
      double prior = get_array_item(i, priors);
      dist_index = raw_to_scaled(prior, 1, scale, offset);
      ++prior_dist[dist_index];
  }

  for (dist_index = 0; dist_index < num_bins; ++dist_index) {
    // Print normalized bin counts
    prior_dist[dist_index] /= num_priors;
    printf("%6.5f\n", prior_dist[dist_index]);
  }

  return 0;
}
