/********************************************************************
 * FILE: alphabet.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4-17-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Define the amino acid and nucleotide alphabets.
 ********************************************************************/
#ifndef ALPHABET_H
#define ALPHABET_H

#include "array.h"
#include "json-writer.h"
#include "string-builder.h"
#include "utils.h"

#include <assert.h>
#include <limits.h>

/*
 * Defines a alphabet
 *
 * IMPORTANT
 * If you need access to this object then please write your accessor
 * either as a macro in alphabet.h or a function in alphabet.c
 * DO NOT just access the structure as while it is exposed here that
 * is purely so macros have access. You should treat it as a black box
 * everywhere else.
 *
 * The only other place that should access the guts is alph-in.c/h which is
 * responsible for loading the alphabet from files and the unit-tests.
 */
typedef struct {
  // as this structure is shared widely it uses reference counting
  // to manage the memory rather than making copies.
  uint64_t references;
  // I need a quick way to check if it is a built in alphabet for code that depends
  // very strongly on the alphabet
  int flags;
  // the human readable name of the alphabet (ie DNA, Protein or something else)
  // if not provided then it will be set to the list of core symbols
  char* alph_name;
  // the number of core symbols (which start at index 1)
  int ncore;
  // the full number of symbols (core + ambiguous which start at index 1)
  int nfull;
  // all the symbols where index 0 is a special error placeholder
  // this is intended to allow the inner loop to avoid branching
  char* symbols;
  // the aliases for all the symbols (starting at index 1)
  char** aliases;
  // the names for all the symbols (starting at index 1, 0 is reserved for error state)
  char** names;
  // the colours for all the symbols in the first 24 bits
  // red = bits 17-24, green = bits 9-16, blue = bits 1-8, unassigned = bits 25-32
  uint32_t* colours;
  // the number of comprising symbols, for core symbols this is 1
  uint8_t* ncomprise;
  // Zero terminated list of core symbols for each symbol (including core symbols)
  uint8_t** comprise;
  // The complements of all symbols. Symbols without complements will return 0.
  uint8_t* complement;
  // Encode a unsigned char to the index in the symbol array for all symbols and their aliases.
  // If a character is not a symbol/alias then it will encode to 0 (error state).
  uint8_t encode[UCHAR_MAX + 1];
  // Encode a unsigned char to the index in the symbol array for core symbols and their aliases.
  // If a charater is not a core symbol/alias then it will encode to 0 (error state).
  uint8_t encode2core[UCHAR_MAX + 1];
  // Same as (encode[sym] - 1) but errors are encoded as wildcards.
  uint8_t encodesafe[UCHAR_MAX + 1];
  // Same as (encode2core[sym] - 1) but errors are encoded as wildcards.
  uint8_t encodesafe2core[UCHAR_MAX + 1];
} ALPH_T;

// structure for translating sequence in one alphabet to another
typedef struct {
  // the source alphabet
  ALPH_T *src_alph;
  // the destination alphabet
  ALPH_T *dest_alph;
  // the number of symbols in the source alphabet to group for a translation
  uint8_t src_nsym;
  // the number of symbols in the destination alphabet output by a translation
  uint8_t dest_nsym;
  // Translate src_nsym symbols in the src_alph into
  // dest_nsym symbols in the dest_alph.
  uint32_t *xlate;
} XLATE_T;

// a background calculation object
typedef struct bgcalc BGCALC_T;

#define ALPH_FLAG_BUILTIN 1

#define ALPH_FLAG_EXTENDS 6
#define ALPH_FLAG_EXTENDS_RNA 2
#define ALPH_FLAG_EXTENDS_DNA 4
#define ALPH_FLAG_EXTENDS_PROTEIN 6

#define ALPH_CASE_INSENSITIVE 8

/*
 * Compare 2 pointers to single symbols.
 * Order: letters < numbers < symbols
 */
int alph_sym_cmp(const void *v1, const void *v2);

/*
 * Compare two NUL terminated symbol strings.
 * Order: length=1 < longest < shortest,  letters < numbers < symbols
 */
int alph_str_cmp(const void *v1, const void *v2);

/*
 * alph_dna
 * Get the built-in DNA alphabet.
 */
ALPH_T* alph_dna();

/*
 * alph_rna
 * Get the built-in RNA alphabet.
 */
ALPH_T* alph_rna();

/*
 * alph_protein
 * Get the built-in protein alphabet.
 */
ALPH_T* alph_protein();

/*
 * alph_generic
 * Get a alphabet with all the passed characters as core symbols.
 */
ALPH_T* alph_generic(const char *core_symbols);

/*
 * alph_load
 * load the alphabet from a file.
 */
ALPH_T* alph_load(const char *filename, bool verbose);

/*
 * alph_pick
 * Evaluate a list of alphabets against a list of symbols and counts and
 * return the index of the alphabet that best represents the counts.
 */
int alph_pick(int nalphs, ALPH_T **alphs, char* symbols, int64_t* counts);

/*
 * alph_guess
 * Evaluate the standard alphabets against a list of symbols and counts and
 * return the alphabet that best represents the counts.
 *
 * Note that this method requires initilizing all the standard alphabets
 * which means it should not be used in inner loops - consider using
 * alph_pick instead.
 */
ALPH_T* alph_guess(char* symbols, int64_t* counts);

/*
 * alph_hold
 * Increment the reference count.
 * This should be paired with a call to alph_release to enable
 * freeing of the memory when it is unused.
 */
ALPH_T* alph_hold(ALPH_T *alphabet);

/*
 * alph_release
 * Deincrement the reference count. 
 * If the reference count reaches zero then the alphabet is destroyed.
 * There should be a call to alph_release to pair the initial creation
 * of the alphabet and every subsequent call to alph_hold.
 */
void alph_release(ALPH_T *alphabet);

/*
 * alph_print_header
 * Print the alphabet header to a file.
 * This is separate mainly so MEME can make it look pretty
 * by wrapping in lines of stars.
 */
void alph_print_header(ALPH_T *alphabet, FILE *out);

/*
 * alph_print
 * print the alphabet to a file.
 */
void alph_print(ALPH_T *alphabet, bool header, FILE *out);

/*
 * alph_print_xml
 * print the alphabet to a XML file.
 */
void alph_print_xml(ALPH_T *alphabet, char *tag, char *pad, char *indent, FILE *out);

/*
 * alph_print_json
 * print the alphabet to a json writer.
 */
void alph_print_json(ALPH_T *alph, JSONWR_T* jsonwr);

/*
 * alph_equal
 * test if two alphabets are the same.
 */
bool alph_equal(const ALPH_T *a1, const ALPH_T *a2);

/*
 * alph_core_subset
 * test if the core of one alphabet is the subset of another
 *  0: not subset
 * -1: subset but complements don't match
 *  1: subset with matching complements
 */
int alph_core_subset(const ALPH_T *sub_alph, const ALPH_T *super_alph);

/*
 * alph_extends_x
 * Allows an alphabet to specify one standard alphabet that it extends.
 * The alphabet must contain all core symbols and complements of the alphabet
 * being extended.
 */
#define alph_extends(alph) (((alph)->flags & ALPH_FLAG_EXTENDS) != 0)
#define alph_extends_rna(alph) (((alph)->flags & ALPH_FLAG_EXTENDS) == ALPH_FLAG_EXTENDS_RNA)
#define alph_extends_dna(alph) (((alph)->flags & ALPH_FLAG_EXTENDS) == ALPH_FLAG_EXTENDS_DNA)
#define alph_extends_protein(alph) (((alph)->flags & ALPH_FLAG_EXTENDS) == ALPH_FLAG_EXTENDS_PROTEIN)

/*
 * alph_is_builtin_x
 * When you really need to know if the alphabet is exactly one of the built-in
 * alphabets then you can use these tests.
 */
#define alph_is_builtin(alph) (((alph)->flags & ALPH_FLAG_BUILTIN) != 0)
#define alph_is_builtin_rna(alph) (alph_is_builtin(alph) && alph_extends_rna(alph))
#define alph_is_builtin_dna(alph) (alph_is_builtin(alph) && alph_extends_dna(alph))
#define alph_is_builtin_protein(alph) (alph_is_builtin(alph) && alph_extends_protein(alph))

/*
 * Get the name of the alphabet - useful for error messages.
 * Allows NULL alphabets.
 */
static inline const char* alph_name(const ALPH_T *alph) {
  return (alph != NULL ? alph->alph_name : "undefined");
}

/*
 * Get the size of the alphabet
 * pairs = the count of complementary pairs
 * core = just the core letters
 * wild = core letters + wildcard
 * all = core letter + ambiguous letter (without gap)
 * ambig = ambiguous letters
 * full = core letters + ambigous letters (including gap)
 */
int alph_size_pairs(const ALPH_T *a);
#define alph_size_core(alph) ((alph)->ncore)
#define alph_size_wild(alph) ((alph)->ncore + 1)
#define alph_size_all(alph) ((alph)->nfull - ((alph)->ncomprise[(alph)->nfull] ? 0 : 1))
#define alph_size_ambig(alph) ((alph)->nfull - (alph)->ncore)
#define alph_size_full(alph) ((alph)->nfull)

/*
 * Check if the alphabet has a complement
 * Assume that if the first symbol has a complement then the alphabet
 * is complementable (need to refine this!).
 */
#define alph_has_complement(alph) ((alph)->ncore > 0 && (alph)->complement[1] != 0)

/*
 * Check if the alphabet is case insensitive
 */
#define alph_is_case_insensitive(alph) (((alph)->flags & ALPH_CASE_INSENSITIVE) != 0)

/*
 * Get the alphabet index of a letter. 
 *
 * Assumes that the alphabet is valid and that the letter is
 * represented by one byte and does not exceed 255 (if treated as unsigned)
 */
#define alph_index(alph, letter) ((alph)->encode[(uint8_t)(letter)] - 1)

/*
 * Get the alphabet index of a letter if it is in the core alphabet. 
 *
 * Assumes that the alphabet is valid and that the letter is
 * represented by one byte and does not exceed 255 (if treated as unsigned)
 */
#define alph_indexc(alph, letter) ((alph)->encode2core[(uint8_t)(letter)] - 1)

/*
 * Get the alphabet index of a letter or the index of the wildcard when
 * the letter is unknown.
 *
 * Assumes that the alphabet is valid and that the letter is
 * represented by one byte and does not exceed 255 (if treated as unsigned)
 */
#define alph_encode(alph, letter) ((alph)->encodesafe[(uint8_t)(letter)])

/*
 * Get the alphabet index of a letter if it is in the core alphabet or the
 * index of the wildcard when the letter is not core. 
 *
 * Assumes that the alphabet is valid and that the letter is
 * represented by one byte and does not exceed 255 (if treated as unsigned)
 */
#define alph_encodec(alph, letter) ((alph)->encodesafe2core[(uint8_t)(letter)])


/*
 * Get the alphabet index of the wildcard
 */
#define alph_wild(alph) ((alph)->ncore)

/*
 * Get the complement index of a index. 
 */
#define alph_complement(alph, index) ((alph)->complement[(index) + 1] - 1)

/*
 * Get the letter at an index of the alphabet
 */
#define alph_char(alph, index) ((alph)->symbols[(index) + 1])

/*
 * Checks to see if the symbol name is actually interesting or
 * just a stringified version of the symbol.
 */
static inline bool alph_has_sym_name(ALPH_T *alph, int index) {
  char *name;
  name = alph->names[index + 1];
  return !(name[0] == alph->symbols[index + 1] && name[1] == '\0');
}

/*
 * Get the symbol name at an index of the alphabet
 */
#define alph_sym_name(alph, index) ((alph)->names[(index) + 1])

/*
 * Get an ID that can represent the symbol in XML as an ID field
 * or even as an attribute.
 */
static inline char* alph_xml_id(ALPH_T *alph, int index, STR_T* buffer) {
  char s;
  s = alph->symbols[index + 1];
  str_clear(buffer);
  if ((s >= 'A' && s <= 'Z') || (s >= 'a' && s <= 'z')) {
    str_appendf(buffer, "%c", s);
  } else if (s >= '0' && s <= '9') {
    str_appendf(buffer, "n%c", s);
  } else {
    // any other symbols get encoded as hexadecimal
    str_appendf(buffer, "x%2X", (uint8_t)s);
  }
  return str_internal(buffer);
}

/*
 * Get the number of comprising core symbols at an index of the alphabet
 */
#define alph_ncomprise(alph, index) ((alph)->ncomprise[(index) + 1])

/*
 * Get the index of a comprising core symbol for an index of the alphabet
 */
#define alph_comprise(alph, index, cindex) ((alph)->comprise[(index) + 1][cindex] - 1)

/*
 * Get the colour at an index of the alphabet
 */
#define alph_colour(alph, index) ((alph)->colours[(index) + 1])
#define alph_colour_r(alph, index) (((alph)->colours[(index) + 1] >> 16) & 0xFF)
#define alph_colour_g(alph, index) (((alph)->colours[(index) + 1] >> 8) & 0xFF)
#define alph_colour_b(alph, index) ((alph)->colours[(index) + 1] & 0xFF)

/*
 * Get the aliases for a symbol index.
 */
#define alph_aliases(alph, index) ((alph)->aliases[(index) + 1])

/*
 * Get the alphabet string (no ambiguous characters)
 */
const char* alph_string(ALPH_T *alph, STR_T *buffer);

/*
 * Get the wildcard letter of the alphabet
 */
#define alph_wildcard(alph) ((alph)->symbols[(alph)->ncore + 1])

/*
 * Determine if a letter is known in this alphabet
 */
#define alph_is_known(alph, letter) ((alph)->encode[(uint8_t)(letter)] != 0)

/*
 * Determine if a letter is a concrete representation
 */
#define alph_is_concrete(alph, letter) ((alph)->encode2core[(uint8_t)(letter)] != 0)
#define alph_is_core(alph, letter) ((alph)->encode2core[(uint8_t)(letter)] != 0)

/*
 * Determine if a letter is ambiguous
 */
#define alph_is_ambiguous(alph, letter) ((alph)->encode[(uint8_t)(letter)] > (alph)->ncore)

/*
 * Determine if a letter is a wildcard for this alphabet
 */
#define alph_is_wildcard(alph, letter) ((alph)->encode[(uint8_t)(letter)] == ((alph)->ncore + 1))


/*
 * Determine if a letter is the main way of representing the symbol.
 */
bool alph_is_prime(ALPH_T *alph, char letter);

/*
 * Get the number of ambiguous letters in a string.
*/
int get_num_ambiguous_letters(ALPH_T *alph, char *string, int length);

/*
 *  Tests the letter against the alphabet. If the alphabet is unknown
 *  it attempts to find a predefined alphabet that matches the letter.
 *  For simplicy this assumes you will pass indexes in asscending order.
 *  Returns false if the letter is unacceptable
 */
bool alph_test(ALPH_T **alpha, int index, char letter);

/*
 *  Tests the alphabet string and attempts to return the ALPH_T*
 *  for one of the predefined alphabets.
 *  If the alphabet is from some buffer a max size can be set
 *  however it will still only test until the first null byte.
 *  If the string is null terminated just set a max > 20 
 */
ALPH_T* alph_type(const char *alphabet, int max);

/*
 * Calculate the values of the ambiguous letters from sums of
 * the normal letter values.
 */
void calc_ambigs(ALPH_T *alph, bool log_space, ARRAY_T* array);

/*
 * Replace the elements an array of frequences with the average
 * over complementary bases.
 */
void average_freq_with_complement(ALPH_T *alph, ARRAY_T *freqs);

/*
 * Load the frequencies that used to be in hash_alph.h
 *
 */
ARRAY_T* get_mast_frequencies(ALPH_T *alph, bool has_ambigs, bool translate);

/*
 * Load the non-redundant database frequencies into the array.
 * FIXME broken function, currently aliased to get_uniform_frequencies
 */
ARRAY_T* get_nrdb_frequencies(ALPH_T *alph, ARRAY_T* freqs);

/*
 * Load uniform frequencies into the array.
 */
ARRAY_T* get_uniform_frequencies(ALPH_T *alph, ARRAY_T* freqs);

/*
 * Add pseudocount and renormalize a frequency array
*/
void normalize_frequencies(ALPH_T *alph, ARRAY_T *freqs, double pseudo);

/*
 * Load background file frequencies into the array.
 */
ARRAY_T* get_file_frequencies(ALPH_T *alph, char *bg_filename);

/*
 * Loads a markov model from a background file.
 */
ARRAY_T* load_markov_model(ALPH_T *alph, int* order, const char *bg_filename);

/*
 * Load background file frequencies into the array.
 */
ARRAY_T* load_markov_model_without_alph(const char *bg_filename, int *order, char **syms);

/*
 * Calculates a markov model from sequence data.
 *
 * When *save is NULL then it will be initialised.
 *
 * When seq is specified then the counts for the contained chains will be added
 * to the calculation and the function will return NULL.
 *
 * If join_seq is true then there will be considered to be no gap between the last
 * passed sequence and this one.
 *
 * When seq is not specified then the result of the background calculation will
 * be returned in the array and the save object will be deallocated.
 *
 */
ARRAY_T* calculate_markov_model(ALPH_T *alph, int order,
    double epsilon, bool join_seq, const char *seq, BGCALC_T** save);

/*
 * Averages background probability tuples with their reverse complement.
 * Expects only the core alphabet - this must be applied before extending.
 */
void average_rc_markov_model(ALPH_T *alph, int order, ARRAY_T *bg);

/*
 * resize_markov_model
 * Allow extending a markov model to allow for a larger alphabet - for example
 * to add ambiguous symbols which can be derived from the existing ones.
 */
void resize_markov_model(int asize0, int asize1, ARRAY_T *tuples, int *order_p);

/*
 * Method used to extend the markov model to includ the ambigous characters
 */
typedef enum ambig_calc {
  SUM_FREQS,
  SUM_LOGS,
  AVG_FREQS
} AMBIG_CALC_EN;

/*
 * Given a markov model with frequencies for the core letters, calculate the
 * values of the ambiguous letters from sums of the core letter values.
 * As an alterntive to having all the ambiguous values it can be limited
 * to the wildcard only which is assumed to be the first symbol after the core
 * symbols.
 */
void extend_markov_model(ALPH_T* alph, bool wildcard_only, AMBIG_CALC_EN method, ARRAY_T* tuples);

/*
 * Given a frequency of seeing any ambiguous symbol calculate frequencies for
 * all entries in the markov model.
 *
 * Zero order ambigous entries are filled in as:
 *    Pr(ambig) = ambig_fraction / (asize1 - asize0)
 * then they are normalized.
 *
 * Higher order entries are calculated based on the previous order frequencies and zero order frequencies as:
 * Pr(word | symbol) = Pr(word) . Pr(symbol)
 * where "word" is a list of symbols.
 * After entries for a order have been filled in they are normalized.
 */
void extrapolate_markov_model(int asize0, int asize1, double ambig_fraction, ARRAY_T* tuples);

/*
 * Get log of the probability of each position in a sequence given a 
 * Markov background model.
 *
 * Note that the Markov model must be normalized for each prefix.
 * 
 * Returns the cumulative background as an array:
 *    logcumback_i = 0, i=0
 *        = log Pr(s_{0,i-1} | H_0), otherwise.
 * and the total (log) cumulative background probability.
 * 
 * The probability of any length-w substring starting at position i in the
 * sequence (in the context of the string) can then be computed as:
 *    last_p = i+w-1;
 *    log_p = logcumback[last_p+1] - logcumback[i];
 * 
 * The background model, H_0, is a Markov model defined by the
 * order, n, and the conditional probabilities, a_cp, where
 *    a_cp[s2i(wa)] = Pr(a | w).
 */
double calculate_log_cumulative_background(ALPH_T *alph, bool wildcard_only,
    int order, ARRAY_T *a_cp, const char *seq, double *logcumback);

/*
 * Get a background distribution
 *  - by reading values from a file if filename is given, or
 *  - equal to the NRDB frequencies if filename is NULL.
 */
ARRAY_T* get_background(ALPH_T *alph, char* bg_filename);

/*
 * Take the counts from each ambiguous character and evenly distribute
 * them among the corresponding concrete characters.
 *
 * This function operates in log space.
 */
void dist_ambigs(ALPH_T *alph, ARRAY_T* freqs);

/*
 * complement_swap_freqs
 * swap complementary positions between arrays and recalulate ambiguous values
 * when the arrays have extra allocated space.
 * You may pass the same array twice to just complement it.
 */
void complement_swap_freqs(ALPH_T *alph, ARRAY_T* a1, ARRAY_T* a2);

/*
 * translate_seq
 * Converts a sequence in place to a more restricted form:
 * ALPH_NO_ALIASES - convert any alias symbols into the primary representation.
 * ALPH_NO_AMBIGS - convert any ambiguous symbols into the main wildcard symbol.
 * ALPH_NO_UNKNOWN - convert any unknown symbols into the main wildcard symbol.
 */
#define ALPH_NO_ALIASES 1
#define ALPH_NO_AMBIGS 2
#define ALPH_NO_UNKNOWN 4
void translate_seq(ALPH_T *alph, char *sequence, int flags);

/*
 * invcomp_seq
 * Converts a sequence in place to its reverse complement.
 */
void invcomp_seq(ALPH_T *alph, char *sequence, long length);

/*
 * comp_sym
 * Converts a symbol to its reverse complement.
 */
// TODO We might have to add another lookup table to the ALPH_T to make this faster
#define comp_sym(alph, sym) ((alph)->symbols[(alph->complement[(alph)->encode[(uint8_t)(sym)]])])

/*
 * dhash_seq
 * Hash sequence to index in two-letter hashed alphabet.
 * If translating, hash two translated letters (eg codons) to two-letter index.
 *
 */
int* dhash_seq(ALPH_T *alph, XLATE_T *trans,  bool full_alph, char *sequence, long length);

/*
 * xlate_dna2protein
 * Get a translation of DNA to protein.
 */
XLATE_T* xlate_dna2protein();

/*
 * xlate_destroy
 * Destroys the passed alphabet translator object.
 */
void xlate_destroy(XLATE_T *translator);

/*
 * xlate_print
 * Print out all translations that don't equal the wildcard.
 */
void xlate_print(XLATE_T *translator, FILE *out);

/*
 * xlate_src_alph
 * Get the source alphabet of the translator.
 */
#define xlate_src_alph(translator) ((translator)->src_alph)

/*
 * xlate_src_nsyms
 * Get the number of symbols that will be translated.
 */
#define xlate_src_nsyms(translator) ((translator)->src_nsym)

/*
 * xlate_dest_alph
 * Get the destination alphabet of the translator.
 */
#define xlate_dest_alph(translator) ((translator)->dest_alph)

/*
 * xlate_pos
 * Convert a subsequence into a position for translating
 */
static inline __attribute__((always_inline)) int xlate_pos(XLATE_T *xlate, bool invcomp, const char *sequence) {
  int i, index;
  ALPH_T *a;
  a = xlate->src_alph;
  if (!invcomp) {
    index = a->encode[(uint8_t)sequence[0]];
    for (i = 1; i < xlate->src_nsym; i++) {
      index = (index * (a->nfull + 1)) + a->encode[(uint8_t)sequence[i]];
    }
  } else {
    index = a->complement[a->encode[(uint8_t)sequence[xlate->src_nsym - 1]]];
    for (i = xlate->src_nsym - 2; i >= 0; i--) {
      index = (index * (a->nfull + 1)) + a->complement[a->encode[(uint8_t)sequence[i]]];
    }
  }
  return index;
}

/*
 * xlate_index
 * Get the index in the destination alphabet for some letters in the source alphabet.
 * The xlate_index2 allows you to pass in the translation postion separately.
 */
#define xlate_index(translator, invcomp, sequence) ((translator)->xlate[xlate_pos(translator, invcomp, sequence)] - 1)
#define xlate_index2(translator, position) ((translator)->xlate[position] - 1)


#endif
