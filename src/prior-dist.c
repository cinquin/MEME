/********************************************************************
 * FILE: prior-dist.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-11-19
 * COPYRIGHT: 2010 UW
 *
 * This file contains the implmentation for the data structures and
 * functions used to store and manipulate the distribution of priors.
 ********************************************************************/

#include <errno.h>
#include "array.h"
#include "utils.h"

typedef struct psp_dist {
  double min;
  double max;
  double median;
  double offset;
  double scale;
  double range;
  ARRAY_T *dist;
} PRIOR_DIST_T;

/***********************************************************************
 * Read an array of unknown length from a file.
 * Caller is resposbile for freeing array.
 ***********************************************************************/
ARRAY_T *read_priors_from_file(
  const char* filename, 
  double *min,
  double *max,
  double *median
) {

  FILE *array_file = fopen(filename, "r");
  if (array_file == NULL) {
    die(
      "Unable to open file: %s.\nError message: %s.\n", 
      filename, 
      strerror(errno)
    );
  }

  // Prior dist files being with min, max, and median
  int num_read = 0;
  num_read = fscanf(array_file, "#min %lf\n", min);
  if (num_read != 1) {
    die("Minimum prior not found in file %s\n", filename);
  }
  num_read = fscanf(array_file, "#max %lf\n", max);
  if (num_read != 1) {
    die("Maximum prior not found in file %s\n", filename);
  }
  num_read = fscanf(array_file, "#median %lf\n", median);
  if (num_read != 1) {
    die("Median prior not found in file %s\n", filename);
  }

  // The remaining values are the frequencies for each bin.
  ATYPE value;
  int array_size = 100;
  ARRAY_T *array = allocate_array(array_size);
  int i_item = 0;
  while ((num_read = fscanf(array_file, ASCAN, &value)) == 1) {
    set_array_item(i_item, value, array);
    ++i_item;
    if (i_item >= array_size) {
      resize_array(array, 2 * array_size);
      array_size = 2 *array_size;
    }
  }

  if (num_read == 0) {
    die("Error reading array at position %d.\n", i_item);
  }
  fclose(array_file);

  resize_array(array, i_item);

  return array;
}

/***********************************************************************
 * Create a new PRIOR_DIST_T object by reading the distribution 
 * from a file.
 ***********************************************************************/
PRIOR_DIST_T *new_prior_dist(const char *filename) {
    
  PRIOR_DIST_T *prior_dist = mm_malloc(sizeof(PRIOR_DIST_T));
  prior_dist->dist = read_priors_from_file(
    filename,
    &(prior_dist->min),
    &(prior_dist->max),
    &(prior_dist->median)
  );
  // The first two items in the array contain
  // the minimum and maximium priors.
  prior_dist->range = get_array_length(prior_dist->dist);
  prior_dist->offset = prior_dist->min;
  prior_dist->scale = (prior_dist->range - 1) / (prior_dist->max - prior_dist->min);

  return prior_dist;
}

/***********************************************************************
 * Free a PRIOR_DIST_T object.
 ***********************************************************************/
void free_prior_dist(PRIOR_DIST_T *prior_dist) {
  free_array(prior_dist->dist);
  myfree(prior_dist);
}

/***********************************************************************
 * Get minimum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_minimum(PRIOR_DIST_T *prior_dist) {
  return prior_dist->min;
}

/***********************************************************************
 * Get maximum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_maximum(PRIOR_DIST_T *prior_dist) {
  return prior_dist->max;
}

/***********************************************************************
 * Get median prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_median(PRIOR_DIST_T *prior_dist) {
  return prior_dist->median;
}

/***********************************************************************
 * Get the array containing distribution from PRIOR_DIST_T object.
 * The caller should not free this array. It will be freed when
 * the PRIOR_DIST_T is freed.
 ***********************************************************************/
ARRAY_T *get_prior_dist_array(PRIOR_DIST_T *prior_dist) {
  return prior_dist->dist;
}

/***********************************************************************
 * Get the length of the array containing the distribution from 
 * a PRIOR_DIST_T object.
 ***********************************************************************/
int get_prior_dist_length(PRIOR_DIST_T *prior_dist) {
  return get_array_length(prior_dist->dist);
}

/***********************************************************************
 * Get the offset for converting an index into the prior dist array
 * into a raw value and visa versa
 ***********************************************************************/
double get_prior_dist_offset(PRIOR_DIST_T *prior_dist) {
  return prior_dist->offset;
}
/***********************************************************************
 * Get the scale for converting an index into the prior dist array
 * into a raw value and visa verse
 * a PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_scale(PRIOR_DIST_T *prior_dist) {
  return prior_dist->scale;
}

/**************************************************************************
*	get_min_lo_priors
*
*	Returns the minimum log-odds prior from the disbtribution of priors.
*
**************************************************************************/
double get_min_lo_prior(PRIOR_DIST_T *prior_dist, double alpha) {
  double min_prior = prior_dist->min;
  return my_log2((alpha * min_prior) / (1.0L - alpha * min_prior));
}

/**************************************************************************
*	get_max_lo_priors
*
*	Returns the maximum log-odds prior from the disbtribution of priors.
*
**************************************************************************/
double get_max_lo_prior(PRIOR_DIST_T *prior_dist, double alpha) {
  double max_prior = prior_dist->max;
  return my_log2((alpha * max_prior) / (1.0L - alpha * max_prior));
}

