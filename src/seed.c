/*
 * @file seed.c 
 *
 * The functions for the SEED object.
 *
 * $Id$
 *
 */


#include "seed.h"
#include "alphabet.h"
#include <assert.h>

#define trace(Y)     fprintf(stderr,Y)

/**
 * new_seed generates a new SEED object.
 * \return A pointer to the new object.
 */
SEED *new_seed(
  char *str_seed,   ///< An ascii representation of the seed.
  double score,     ///< Score of the seed as a starting point for local search.
  int iseq,         ///< Index of sequence seed is from.
  int ipos          ///< Position in sequence seed is from.
)
{
  // Set aside memory for the seed and its contents:
  SEED *seed = NULL;
  Resize(seed, 1, SEED);

  // Determine w, in order to allocate memory for the seed representation:
  int w = 0;
  while (str_seed[w] != '\0') {
    w++;
  }

  seed->str_seed = NULL;
  Resize(seed->str_seed, w+1, char); 

  // Set the contents of the seed
  set_seed(seed, str_seed, score, iseq, ipos);
  return(seed);
}


/**
 * compare_seed 
 * 
 * This function compares two seeds. The rule is:
 *      Returns 1   if  score1 > score2,
 *      Returns -1  if  score1 < score2,
 *      If score1 == score2 then
 *        Returns 1  if  iseq_1 < iseq_2
 *        Returns -1 if  iseq_1 > iseq_2
 *      else if iseq_1 == iseq_2 then
 *        Returns 1  if  ipos_1 < ipos_2
 *        Returns -1 if  ipos_1 > ipos_2
 *      else
 *        Returns 0
 *
 */
int compare_seed(
  SEED *s1,         ///< pointer to a SEED object
  SEED *s2          ///< pointer to a SEED object
)
{
  double score_1 = get_seed_score(s1);
  double score_2 = get_seed_score(s2);
  int iseq_1 = s1->iseq;
  int iseq_2 = s2->iseq;
  int ipos_1 = s1->ipos;
  int ipos_2 = s2->ipos;
  if (score_1 != score_2) {
    return (score_1 > score_2 ? 1 : -1);
  } else if (iseq_1 != iseq_2) {
    return (iseq_1 < iseq_2 ? 1 : -1);
  } else if (ipos_1 != ipos_2) {
    return (ipos_1 < ipos_2 ? 1 : -1);
  } else if (strcmp(s1->str_seed, s2->str_seed) != 0) {
    return (strcmp(s1->str_seed, s2->str_seed) < 0 ? 1 : -1);
  } else {
    // the strings are the same
    return 0;
  }
} // compare_seed


/**
 * get_seed_score retrieves the score for a given seed.
 * \return The score for the input seed.
 */
double get_seed_score(
  SEED *seed         ///< The seed object whose score is being retrieved.
)
{
  return(seed->score);
}


/**
 * get_e_seed
 *
 * Converts the string representation of the specified seed into an integer
 * encoded representation, and returns a pointer to the start of the resulting
 * array. Memory is allocated for the array in this function. Deallocation must
 * be controlled by the caller of this function. Note that the length of the
 * array can get retrieved via "get_width()" for the same SEED object.
 *
 * \return The e_seed for the input seed.
 */
uint8_t *get_e_seed(
  ALPH_T *alph,
  SEED *seed         ///< The seed object whose e_seed is being retrieved.
)
{
  int w = get_width(seed);
  // Allocate memory:
  uint8_t *e_seed = NULL;
  Resize(e_seed, w, uint8_t);

  // Perform conversion:
  int seed_idx;
  for (seed_idx = 0; seed_idx < w; seed_idx++) {
    char lett = (seed->str_seed)[seed_idx];
    uint8_t e_lett = alph_encode(alph, lett);
    e_seed[seed_idx] = e_lett;
  }

  // Return the result:
  return e_seed;
}


/**
 * set_seed sets the fields for the input pointer to a seed object.
 */
void set_seed(
  SEED *seed,        ///< The seed object whose fields are to be set.
  char *str_seed,    ///< An ascii representation of the seed.
  double score,      ///< Score of the seed as a local search starting point.
  int iseq,          ///< Index of sequence seed is from.
  int ipos           ///< Position in sequence seed is from.
)
{
  // It is valid for a new seed to be constructed with a NULL ascii string:
  if (seed == NULL) {
    seed->str_seed = NULL;
  } else {
    int seed_idx = 0;
    // Initialise the SEED ascii representation:
    while (str_seed[seed_idx] != '\0') {
      seed->str_seed[seed_idx] = str_seed[seed_idx];
      seed_idx++;
    }
    seed->str_seed[seed_idx] = '\0';
  }

  seed->score = score;
  seed->iseq = iseq;
  seed->ipos = ipos;
}


/**
 * copy_seed copies the fields of the input seed object into the fields of
 * the a new seed object.
 * An EXACT copy of the seed is made, such that compare_seed between the
 * original and the copy yields zero. 
 * \return A pointer to a new seed that is a copy of the original.
 */
SEED *copy_seed(
  SEED *orig_object  ///< The existing seed object that will be copied
)
{
  SEED *orig_seed = orig_object;
  SEED *seed = new_seed(
    get_str_seed(orig_seed), 
    get_seed_score(orig_seed),
    orig_seed->iseq,
    orig_seed->ipos
  );

  return(seed);
}


/** 
 * free_seed destroys the input seed object. 
 */
void free_seed(
  SEED *obj         ///< The seed object to be destroyed
)
{
  SEED *seed = obj;
  myfree(seed->str_seed);
  myfree(seed);
}


/**
 * print_seed prints a representation of the seed to standard out, thus
 * providing a way to debug code related to seed objects.
 */
void print_seed(
  FILE *outfile,     ///< The file to print to
  SEED *seed         ///< Seed to be printed.
)
{
  fprintf(outfile, "--------------------------------------------------\n");
  fprintf(outfile, "SEED:\n");
  if (seed == NULL) {
    fprintf(outfile, "NULL\n");
  } else {
    // Print the stored string representation of the seed:
    fprintf(outfile, "str_seed: ");
    if (seed->str_seed != NULL) {
      fprintf(outfile, "%s\n", seed->str_seed);
    } else {
      fprintf(outfile, "NULL\n");
    }

    fprintf(outfile, "score: %f\n", seed->score);
    fprintf(outfile, "iseq: %d\n", seed->iseq);
    fprintf(outfile, "ipos: %d\n", seed->ipos);
  }
  fprintf(outfile, "**************************************************\n\n\n");
}


/**
 * get_str_seed retrieves the ascii string representation of this e_seed.
 * /return ascii representation of the e_seed.
 */
char *get_str_seed(
  SEED *seed         ///< The seed whose representation is being retrieved.
)
{
  assert (seed != NULL);

  return(seed->str_seed);
}


/**
 * get_width
 *
 * Retrieves the length of the seed. Note that this is an 0(w) operation
 * because w is not explicitly stored in the SEED object and is instead
 * determined by analysing str_seed.
 *
 * /return length of str_seed.
 */
int get_width(
  SEED *seed         ///< The seed whose length is being retrieved.
)
{
  // w is equal to the length of the string representation of the seed:
  return strlen(seed->str_seed);

  /*
  int seed_idx = 0;
  // Iterate through str_seed to determine w:
  while(seed->str_seed[seed_idx] != '\0') {
    seed_idx++;
  }
  return seed_idx;
  */
}


/**
 * to_str_seed
 *
 * This function converts an integer encoded representation of a seed into an
 * ascii representation of it. Memory for the string is dynamically allocated
 * here, and it is the caller's responsibility to later free that memory.
 */
char *to_str_seed(
  ALPH_T *alph,      // alphabet
  uint8_t *e_seed,   // Integer encoded representation.
  int w              // The length of the string.
)
{
  char *str_seed = NULL;
  Resize(str_seed, w+1, char);

  int seed_idx;
  for (seed_idx = 0; seed_idx < w; seed_idx++) {
    str_seed[seed_idx] = alph_char(alph, e_seed[seed_idx]);
  }
  str_seed[w] = '\0';
  return str_seed;
}


/**
 * to_e_seed
 *
 * Converts a string representation of a seed into its integer encoding, storing
 * the length of the seed in the given parameter and returning the integer
 * array.
 */
uint8_t *to_e_seed (
  ALPH_T *alph,
  char *str_seed,    ///< ASCII representation
  int *seed_len      ///< The length of the resulting e_seed
)
{
  uint8_t *e_seed = NULL;
  *seed_len = strlen(str_seed);
  Resize(e_seed, *seed_len, uint8_t);

  int seed_idx;
  for (seed_idx = 0; seed_idx < *seed_len; seed_idx++) {
    e_seed[seed_idx] = alph_encode(alph, str_seed[seed_idx]);
  }

  return e_seed;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
