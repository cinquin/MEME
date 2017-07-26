/*
 * @file seed.h 
 *
 * A definition of the interface for SEED objects.
 *
 * $Id$
 *
 * $Log$
 * Revision 1.1.2.12  2006/03/10 06:58:01  twhitington
 * print_heaps and "eval_best_seed" command line switches added.
 *
 * Revision 1.1.2.11  2006/03/02 04:55:15  twhitington
 * Removed old FIXME comment regarding string and integer rep storage.
 *
 * Revision 1.1.2.10  2006/03/02 04:31:12  twhitington
 * Changed get_w to get_width.
 *
 * Revision 1.1.2.9  2006/02/22 01:53:32  twhitington
 * Changed SEED objects so they don't store e_seed, and are instantiated with
 * ascii strings instead of encoded representations.
 *
 * Revision 1.1.2.8  2006/02/01 02:26:04  twhitington
 * Changed formal parameter names in seed.h.
 *
 * Revision 1.1.2.7  2006/02/01 01:38:46  twhitington
 * Modified subseq7 so that Tcm make test is passed;
 * Catered for rejection of seeds by align_top_subsequences.
 *
 * Revision 1.1.2.6  2006/01/31 08:21:19  tbailey
 * Changed way functions are passed to create_heap
 *
 * Revision 1.1.2.5  2006/01/30 06:06:46  twhitington
 * Introduced serial number variable to seed objects.
 * Introduced str_seed variable to seed objects.
 *
 * Revision 1.1.2.4  2006/01/27 07:32:47  twhitington
 * Debugging changes made to subseq7.c, and debugging problems with get_key
 * in seed.c.
 *
 * Revision 1.1.2.3  2006/01/19 04:59:08  twhitington
 * Introduced heap objects to S_POINT objects. Changed get_nsites0 and get_starts accordingly, in starts.c. Corresponding changes were required in the seed interface and the Makefile.
 *
 * Revision 1.1.2.2  2006/01/18 07:09:08  twhitington
 * Implemented SEED object. Also introduced seed-tests program for testing SEED objects.
 *
 * Revision 1.1.2.1  2006/01/16 06:04:36  tbailey
 * Add seed.c and seed.h for branching project.
 *
*/

#ifndef SEED_H
#define SEED_H
#include "macros.h"
#include "stdint.h"
#include "alphabet.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif


/**
 * A SEED type, for encapsulating an alphabet sequence that could could as a
 * starting point for meme.
 */
typedef struct seed SEED;
struct seed{
  char *str_seed;    ///< An ascii representation of the seed.
  double score;      ///< Score of e_seed as a starting point for local search
  int iseq;	     ///< The sequence number index
  int ipos;	     ///< The sequence position index
};


/**
 * new_seed generates a new SEED object.
 * \return A pointer to the new object.
 */
SEED *new_seed(
  char *str_seed,   ///< An ascii representation of the seed.
  double score,     ///< Score of the seed as a starting point for local search.
  int iseq,         ///< Index of sequence seed is from.
  int ipos          ///< Position in sequence seed is from.
);


/**
 * compare_seed 
 * 
 * This function compares two seeds. The rule is:
 * 	Returns 1   if  score1 > score2,
 * 	Returns -1  if  score1 < score2,
 *	If score1 == score2 then	
 * 	  Returns 1   if  iseq_1 < iseq_2
 * 	  Returns -1   if  iseq_1 > iseq_2
 *	else if iseq_1 == iseq_2 then	
 * 	  Returns 1   if  ipos_1 < ipos_2
 * 	  Returns -1   if  ipos_1 > ipos_2
 *	else if ipos_1 == ipos_2 then	
 * 	  Returns -1  if  serial_no_1 > serial_no_2
 * 	  Returns 1   if  serial_no_1 < serial_no_2
 *      else 
 *	  Returns 0
 *
 *
 */
int compare_seed(
  SEED *p1,          ///< pointer to a SEED object
  SEED *p2           ///< pointer to a SEED object
);


/**
 * get_seed_score retrieves the score for a given seed.
 * \return The score for the input seed.
 */ 
double get_seed_score(
  SEED *seed     ///< seed object
);


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
  SEED *seed
);


/**
 * set_seed sets the fields for the input pointer to a seed object.
 */
void set_seed(
  SEED *seed,        ///< The seed object whose fields are to be set
  char *str_seed,    ///< An ascii representation of the seed.
  double score,      ///< Score of the seed as a local search starting point.
  int iseq,          ///< Index of sequence seed is from.
  int ipos           ///< Position in sequence seed is from.
);


/**
 * copy_seed copies the fields of the input seed object into the fields of
 * the a new seed object.
 * An EXACT copy of the seed is made, such that compare_seed between the
 * original and the copy yields zero. In particular, note that the
 * serial_number of the new seed is equal to that of the old seed.
 * \return A pointer to a new seed that is a copy of the original.
 */
SEED *copy_seed(
  SEED *orig_object  ///< The existing seed object that will be copied
);


/** 
 * free_seed destroyes the input seed object. 
 */ 
void free_seed(
  SEED *obj          ///< The seed object to be destroyed
);


/**
 * print_seed prints a representation of the seed to standard out, thus
 * providing a way to debug code related to seed objects.
 */
void print_seed(
  FILE *outfile,     ///< The file to print to
  SEED *seed         ///< Seed to be printed.
);


/**
 * get_str_seed retrieves the ascii string representation of this e_seed.
 * /return ascii representation of the e_seed.
 */
char *get_str_seed(
  SEED *seed          ///< The seed whose representation is being retrieved.
); 


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
);


/**
 * to_str_seed
 *
 * This function converts an integer encoded representation of a seed into an
 * ascii representation of it. Memory for the string is dynamically allocated
 * here, and it is the caller's responsibility to later free that memory.
 */
char *to_str_seed(
  ALPH_T *alph,
  uint8_t *e_seed,      ///< Integer encoded representation.
  int w              ///< The length of the string.
);


uint8_t *to_e_seed (
  ALPH_T *alph,
  char *str_seed,    ///< ASCII representation
  int *seed_len      ///< The length of the resulting e_seed
);

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
