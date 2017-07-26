/********************************************************************
 * FILE: alphabet.c
 * AUTHOR: James Johnson
 * CREATE DATE: 21 July 2011
 * PROJECT: all
 * COPYRIGHT: 2011 UQ
 * DESCRIPTION: Define the amino acid and nucleotide alphabets as
 * well as any utility functions. Many of the functions are rewrites
 * of the original ones which made use of mutable global variables.
 ********************************************************************/

#include "alphabet.h"
#include "alph-in.h"
#include "string-builder.h"
#include "regex-utils.h"
#include "red-black-tree.h"
#include "xlate-in.h"

#include <assert.h>
#include <ctype.h>
#include <limits.h>

// The number of bytes to read per chunk of background file
#define BG_CHUNK_SIZE 100

/*
 * The regular expression to match a background frequency
 * from a background file
 */
const char *BGFREQ_RE = "^[[:space:]]*([a-zA-Z])[[:space:]]+"
    "([+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)[[:space:]]*$";

/*
 * Compare two symbol characters.
 * Takes in pointers to the symbols characters.
 * Order: letters < numbers < symbols
 */
int alph_sym_cmp(const void *v1, const void *v2) {
  char c1, c2;
  bool is_letter_1, is_letter_2, is_num_1, is_num_2;
  c1 = *((const char*)v1);
  c2 = *((const char*)v2);
  // Are they letters? We want letters listed first!
  is_letter_1 = ((c1 >= 'A' && c1 <= 'Z') || (c1 >= 'a' && c1 <= 'z'));
  is_letter_2 = ((c2 >= 'A' && c2 <= 'Z') || (c2 >= 'a' && c2 <= 'z'));
  if (is_letter_1) {
    if (is_letter_2) {
      return c1 - c2; // both are letters
    } else {
      return -1; // v1 is a letter (v2 is not)
    }
  } else if (is_letter_2) {
    return 1; // v2 is a letter (v1 is not)
  }
  // Neither are letters, are they numbers? We want numbers listed before symbols!
  is_num_1 = (c1 >= '0' && c1 <= '9');
  is_num_2 = (c2 >= '0' && c2 <= '9');
  if (is_num_1) {
    if (is_num_2) {
      return c1 - c2; // both are numbers
    } else {
      return -1; // v1 is a number (v2 must be a symbol)
    }
  } else if (is_num_2) {
    return 1; // v2 is a number (v1 must be a symbol)
  }
  // both must be a symbol
  return c1 - c2;
}

/*
 * Compare two NUL terminated symbol strings.
 * Order: length=1 < longest < shortest,  letters < numbers < symbols
 */
int alph_str_cmp(const void *v1, const void *v2) {
  int len1, len2, i;
  len1 = strlen((const char*)v1);
  len2 = strlen((const char*)v2);
  if (len1 != len2) {
    // we want core symbols to go first, followed by ambiguous, longest to shortest
    if (len1 == 1) {
      return -1;
    } else if (len2 == 1) {
      return 1;
    } else {
      return len2 - len1;
    }
  }
  for (i = 0; i < len1; i++) {
    if (((const char*)v1)[i] != ((const char*)v2)[i]) {
      return alph_sym_cmp(((const char*)v1)+i, ((const char*)v2)+i);
    }
  }
  return 0;
}

/*
 * alph_dna
 * Get the built-in DNA alphabet.
 */
ALPH_T* alph_dna() {
  ALPH_READER_T *factory;
  ALPH_T *alphabet;
  factory = alph_reader_create();
  alph_reader_header(factory, 1, "DNA", ALPH_FLAG_BUILTIN | ALPH_FLAG_EXTENDS_DNA);
  alph_reader_core(factory, 'A', NULL, "Adenine", 0xCC0000, 'T');
  alph_reader_core(factory, 'C', NULL, "Cytosine", 0x0000CC, 'G');
  alph_reader_core(factory, 'G', NULL, "Guanine", 0xFFB300, 'C');
  alph_reader_core(factory, 'T', "U", "Thymine", 0x008000, 'A');
  // Note that -1 signifies no assigned colour which allows the factory to choose.
  // Generally this means that ambiguous symbols will be black unless they
  // alias a core symbol like U which will take on the same colour as T.
  alph_reader_ambig(factory, 'W', NULL, "Weak", -1, "AT");
  alph_reader_ambig(factory, 'S', NULL, "Strong", -1, "CG");
  alph_reader_ambig(factory, 'M', NULL, "Amino", -1, "AC");
  alph_reader_ambig(factory, 'K', NULL, "Keto", -1, "GT");
  alph_reader_ambig(factory, 'R', NULL, "Purine", -1, "AG");
  alph_reader_ambig(factory, 'Y', NULL, "Pyrimidine", -1, "CT");
  alph_reader_ambig(factory, 'B', NULL, "Not A", -1, "CGT");
  alph_reader_ambig(factory, 'D', NULL, "Not C", -1, "AGT");
  alph_reader_ambig(factory, 'H', NULL, "Not G", -1, "ACT");
  alph_reader_ambig(factory, 'V', NULL, "Not T", -1, "ACG");
  alph_reader_ambig(factory, 'N', "X.", "Any base", -1, "ACGT");
  // for now treat '.' (aka a gap) as wildcard (though this is a bit odd)
  alph_reader_done(factory);
  if (alph_reader_had_warning(factory) || alph_reader_had_error(factory)) {
    while (alph_reader_has_message(factory)) {
      PARMSG_T *msg;
      msg = alph_reader_next_message(factory);
      parmsg_print(msg, stderr);
      parmsg_destroy(msg);
    }
    fprintf(stderr, "Standard DNA alphabet should not produce warnings or errors!");
    abort();
  }
  alphabet = alph_reader_alphabet(factory);
  alph_reader_destroy(factory);
  return alphabet;
}

/*
 * alph_rna
 * Get the built-in RNA alphabet.
 */
ALPH_T* alph_rna() {
  ALPH_READER_T *factory;
  ALPH_T *alphabet;
  factory = alph_reader_create();
  alph_reader_header(factory, 1, "RNA", ALPH_FLAG_BUILTIN | ALPH_FLAG_EXTENDS_RNA);
  // core symbols
  alph_reader_core(factory, 'A', NULL, "Adenine", 0xCC0000, '\0');
  alph_reader_core(factory, 'C', NULL, "Cytosine", 0x0000CC, '\0');
  alph_reader_core(factory, 'G', NULL, "Guanine", 0xFFB300, '\0');
  alph_reader_core(factory, 'U', "T", "Uracil", 0x008000, '\0');
  // ambiguous symbols
  alph_reader_ambig(factory, 'W', NULL, "Weak", -1, "AU");
  alph_reader_ambig(factory, 'S', NULL, "Strong", -1, "CG");
  alph_reader_ambig(factory, 'M', NULL, "Amino", -1, "AC");
  alph_reader_ambig(factory, 'K', NULL, "Keto", -1, "GU");
  alph_reader_ambig(factory, 'R', NULL, "Purine", -1, "AG");
  alph_reader_ambig(factory, 'Y', NULL, "Pyrimidine", -1, "CU");
  alph_reader_ambig(factory, 'B', NULL, "Not A", -1, "CGU");
  alph_reader_ambig(factory, 'D', NULL, "Not C", -1, "AGU");
  alph_reader_ambig(factory, 'H', NULL, "Not G", -1, "ACU");
  alph_reader_ambig(factory, 'V', NULL, "Not U", -1, "ACG");
  alph_reader_ambig(factory, 'N', "X.", "Any base", -1, "ACGU");
  // for now treat '.' (aka a gap) as wildcard (though this is a bit odd)
  alph_reader_done(factory);
  if (alph_reader_had_warning(factory) || alph_reader_had_error(factory)) {
    while (alph_reader_has_message(factory)) {
      PARMSG_T *msg;
      msg = alph_reader_next_message(factory);
      parmsg_print(msg, stderr);
      parmsg_destroy(msg);
    }
    fprintf(stderr, "Standard RNA alphabet should not produce warnings or errors!");
    abort();
  }
  alphabet = alph_reader_alphabet(factory);
  alph_reader_destroy(factory);
  return alphabet;
}

/*
 * alph_protein
 * Get the built-in protein alphabet.
 */
ALPH_T* alph_protein() {
  ALPH_READER_T *factory;
  ALPH_T *alphabet;
  factory = alph_reader_create();
  alph_reader_header(factory, 1, "Protein", ALPH_FLAG_BUILTIN | ALPH_FLAG_EXTENDS_PROTEIN);
  alph_reader_core(factory, 'A', NULL, "Alanine", 0x0000CC, '\0');
  alph_reader_core(factory, 'R', NULL, "Arginine", 0xCC0000, '\0');
  alph_reader_core(factory, 'N', NULL, "Asparagine", 0x008000, '\0');
  alph_reader_core(factory, 'D', NULL, "Aspartic acid", 0xFF00FF, '\0');
  alph_reader_core(factory, 'C', NULL, "Cysteine", 0x0000CC, '\0');
  alph_reader_core(factory, 'E', NULL, "Glutamic acid", 0xFF00FF, '\0');
  alph_reader_core(factory, 'Q', NULL, "Glutamine", 0x008000, '\0');
  alph_reader_core(factory, 'G', NULL, "Glycine", 0xFFB300, '\0');
  alph_reader_core(factory, 'H', NULL, "Histidine", 0xFFCCCC, '\0');
  alph_reader_core(factory, 'I', NULL, "Isoleucine", 0x0000CC, '\0');
  alph_reader_core(factory, 'L', NULL, "Leucine", 0x0000CC, '\0');
  alph_reader_core(factory, 'K', NULL, "Lysine", 0xCC0000, '\0');
  alph_reader_core(factory, 'M', NULL, "Methionine", 0x0000CC, '\0');
  alph_reader_core(factory, 'F', NULL, "Phenylalanine", 0x0000CC, '\0');
  alph_reader_core(factory, 'P', NULL, "Proline", 0xFFFF00, '\0');
  alph_reader_core(factory, 'S', NULL, "Serine", 0x008000, '\0');
  alph_reader_core(factory, 'T', NULL, "Threonine", 0x008000, '\0');
  alph_reader_core(factory, 'W', NULL, "Tryptophan", 0x0000CC, '\0');
  alph_reader_core(factory, 'Y', NULL, "Tyrosine", 0x33E6CC, '\0');
  alph_reader_core(factory, 'V', NULL, "Valine", 0x0000CC, '\0');
  alph_reader_ambig(factory, 'B', NULL, "Asparagine or Aspartic acid", -1, "ND");
  alph_reader_ambig(factory, 'Z', NULL, "Glutamine or Glutamic acid", -1, "QE");
  alph_reader_ambig(factory, 'J', NULL, "Leucine or Isoleucine", -1, "LI");
  alph_reader_ambig(factory, 'X', "*.", "Any amino acid", -1, "ARNDCEQGHILKMFPSTWYV");
  // for now treat '*' (stop codon) and '.' (gap) as a wildcard
  alph_reader_done(factory);
  if (alph_reader_had_warning(factory) || alph_reader_had_error(factory)) {
    while (alph_reader_has_message(factory)) {
      PARMSG_T *msg;
      msg = alph_reader_next_message(factory);
      parmsg_print(msg, stderr);
      parmsg_destroy(msg);
    }
    fprintf(stderr, "Standard protein alphabet should not produce warnings or errors!");
    abort();
  }
  alphabet = alph_reader_alphabet(factory);
  alph_reader_destroy(factory);
  return alphabet;
}

/*
 * alph_pick
 * Evaluate a list of alphabets against a list of symbols and counts and
 * return the index of the alphabet that best represents the counts.
 */
int alph_pick(int nalphs, ALPH_T **alphs, char* symbols, int64_t* counts) {
  int i, j, idx, prime_seen, all_seen, best_index;
  int64_t prime_count, alt_count, other_count, unknown_count, total_count;
  uint32_t prime_set[4], all_set[4];
  ALPH_T *alph;
  double best_score, best_prime_score, score, prime_score;
  best_index = 0;
  best_score = 0.0;
  best_prime_score = 0.0;
  for (i = 0; i < nalphs; i++) {
    alph = alphs[i];
    prime_count = 0;
    alt_count = 0;
    other_count = 0;
    unknown_count = 0;
    prime_seen = 0; prime_set[0] = 0; prime_set[1] = 0; prime_set[2] = 0; prime_set[3] = 0;
    all_seen = 0; all_set[0] = 0; all_set[1] = 0; all_set[2] = 0; all_set[3] = 0;
    for (j = 0; symbols[j] != '\0'; j++) {
      if (alph_is_known(alph, symbols[j])) {
        idx = alph_index(alph, symbols[j]);
        if ((all_set[idx / 32] & (1 << (idx % 32))) == 0) all_seen++;
        all_set[idx / 32] |= (1 << (idx % 32));
        if (alph_is_core(alph, symbols[j])) {
          if (alph_is_prime(alph, symbols[j])) {
            if ((prime_set[idx / 32] & (1 << (idx % 32))) == 0) prime_seen++;
            prime_set[idx / 32] |= (1 << (idx % 32));
            prime_count += counts[j];
          } else {
            alt_count += counts[j];
          }
        } else if (alph_is_wildcard(alph, symbols[j])) {
          if (alph_is_prime(alph, symbols[j])) {
            prime_count += counts[j];
          } else {
            alt_count += counts[j];
          }
        } else {
          other_count += counts[j];
        }
      } else {
        unknown_count += counts[j];
      }
    }
    total_count = prime_count + alt_count + other_count + unknown_count;
    score = ((double)(prime_count + alt_count) / total_count) * ((double)all_seen / alph_size_core(alph));
    prime_score = ((double)prime_count / total_count) * ((double)prime_seen / alph_size_core(alph));
    if (score > best_score || (score == best_score && prime_score > best_prime_score)) {
      best_index = i;
      best_score = score;
      best_prime_score = prime_score;
    }
  }
  return best_index;
}

/*
 * alph_guess
 * Evaluate the standard alphabets against a list of symbols and counts and
 * return the alphabet that best represents the counts.
 *
 * Note that this method requires initilizing all the standard alphabets
 * which means it should not be used in inner loops - consider using
 * alph_pick instead.
 */
ALPH_T* alph_guess(char* symbols, int64_t* counts) {
  int index;
  ALPH_T* alphs[3];
  alphs[0] = alph_rna();
  alphs[1] = alph_dna();
  alphs[2] = alph_protein();
  index = alph_pick(3, alphs, symbols, counts);
  if (index != 0) alph_release(alphs[0]);
  if (index != 1) alph_release(alphs[1]);
  if (index != 2) alph_release(alphs[2]);
  return alphs[index];
}

/*
 * alph_generic
 * Get a alphabet with all the passed characters as core symbols.
 */
ALPH_T* alph_generic(const char *core_symbols) {
  const char *s;
  ALPH_READER_T *factory;
  ALPH_T *alphabet;
  factory = alph_reader_create();
  alph_reader_header(factory, 1, "", 0);
  for (s = core_symbols; *s; s++) {
    alph_reader_core(factory, *s, NULL, "", 0x000000, '\0');
  }
  // we must have a wildcard even if we don't intend to use it
  alph_reader_ambig(factory, '?', NULL, "wildcard", -1, core_symbols);
  alph_reader_done(factory);
  alphabet = alph_reader_alphabet(factory);
  alph_reader_destroy(factory);
  return alphabet;
}

/*
 * alph_load
 * load the alphabet from a file.
 */
ALPH_T* alph_load(const char *filename, bool verbose) {
  ALPH_READER_T* reader;
  FILE *fp;
  char chunk[100];
  size_t count;
  ALPH_T *alphabet;
  PARMSG_T *message;
  fp = fopen(filename, "r");
  if (fp == NULL) {
    if (verbose) fprintf(stderr, "Failed to open alphabet file: %s\n", filename);
    return NULL;
  }
  reader = alph_reader_create();
  while (!feof(fp) && !ferror(fp)) {
    count = fread(chunk, sizeof(char), sizeof(chunk)/sizeof(char), fp);
    alph_reader_update(reader, chunk, count, feof(fp));
    if (verbose) {
      while (alph_reader_has_message(reader)) {
        message = alph_reader_next_message(reader);
        parmsg_print(message, stderr);
        parmsg_destroy(message);
      }
    }
  }
  if (ferror(fp)) {
    if (verbose) fprintf(stderr, "IO error reading alphabet file: %s\n", filename);
    fclose(fp);
    alph_reader_destroy(reader);
    return NULL;
  }
  fclose(fp);
  alphabet = alph_reader_alphabet(reader);
  alph_reader_destroy(reader);
  return alphabet;
}

/*
 * alph_destroy
 * THIS FUNCTION SHOULD NOT BE CALLED DIRECTLY!
 * Destroys the passed alphabet object. 
 * Prints a warning if the reference count is not zero.
 */
void alph_destroy(ALPH_T *alphabet) {
  int i;
  if (alphabet->references != 0) {
    fputs("WARNING: alphabet destroyed when the reference count was non-zero.\n", stderr);
  }
  for (i = 0; i <= alphabet->nfull; i++) {
    free(alphabet->names[i]);
    free(alphabet->aliases[i]);
    free(alphabet->comprise[i]);
  }
  free(alphabet->alph_name);
  free(alphabet->symbols);
  free(alphabet->names);
  free(alphabet->aliases);
  free(alphabet->colours);
  free(alphabet->ncomprise);
  free(alphabet->comprise);
  free(alphabet->complement);
  // guard against accidental use after free
  // the builtin alphabets rely on this line
  memset(alphabet, 0, sizeof(ALPH_T));
  free(alphabet);
}

/*
 * alph_hold
 * Increment the reference count.
 * This should be paired with a call to alph_release to enable
 * freeing of the memory when it is unused.
 */
ALPH_T* alph_hold(ALPH_T *alphabet) {
  if (alphabet == NULL) return NULL;
  assert(alphabet->references > 0);
  alphabet->references++;
  return alphabet;
}

/*
 * alph_release
 * Deincrement the reference count. 
 * If the reference count reaches zero then the alphabet is destroyed.
 * There should be a call to alph_release to pair the initial creation
 * of the alphabet and every subsequent call to alph_hold.
 */
void alph_release(ALPH_T *alphabet) {
  assert(alphabet->references > 0);
  alphabet->references--;
  if (alphabet->references <= 0) alph_destroy(alphabet);
}

/*
 * print_name
 * Print out the name of a symbol
 */
static void print_name(FILE *out, const char *name) {
  const char *c;
  fputc('"', out);
  for (c = name; *c != '\0'; c++) {
    switch (*c) {
      case '"': fputs("\\\"", out); break;
      case '/': fputs("\\/", out); break;
      case '\\': fputs("\\\\", out); break;
      default: fputc(*c, out);
    }
  }
  fputc('"', out);
}

/*
 * print_name
 * Print out the symbol details
 */
static void print_symbol(FILE *out, const char symbol, const char *name, int32_t colour) {
  fprintf(out, "%c", symbol);
  if (name[0] != symbol || name[1] != '\0') {
    fputc(' ', out);
    print_name(out, name);
  }
  if (colour != 0) fprintf(out, " %06X", colour);
}

/*
 * comprise_cmp
 * Unlike normal strcmp we want to test the length first
 * so ambiguous symbols with fewer comprising symbols are first.
 */
static int comprise_cmp(const void *s1, const void *s2) {
  int result;
  result = strlen((char*)s1) - strlen((char*)s2);
  if (result == 0) result = strcmp((char*)s1, (char*)s2);
  return result;
}

/*
 * alph_print_header
 * Print the alphabet header to a file.
 * This is separate mainly so MEME can make it look pretty
 * by wrapping in lines of stars.
 */
void alph_print_header(ALPH_T *alphabet, FILE *out) {
  fputs("ALPHABET ", out);
  print_name(out, alphabet->alph_name);
  if (alph_extends_rna(alphabet)) {
    fputs(" RNA-LIKE", out);
  } else if (alph_extends_dna(alphabet)) {
    fputs(" DNA-LIKE", out);
  } else if (alph_extends_protein(alphabet)) {
    fputs(" PROTEIN-LIKE", out);
  }
  fputs("\n", out);
}

/*
 * alph_print
 * print the alphabet to a file.
 */
void alph_print(ALPH_T *alphabet, bool header, FILE *out) {
  int i, j, c;
  char *comprise;
  // Print header
  if (header) alph_print_header(alphabet, out);
  // Print core symbols with 2 way complements
  for (i = 1; i <= alphabet->ncore; i++) {
    c = alphabet->complement[i];
    if (alphabet->complement[c] == i && i < c) {
      print_symbol(out, alphabet->symbols[i], alphabet->names[i], alphabet->colours[i]);
      fputs(" ~ ", out);
      print_symbol(out, alphabet->symbols[c], alphabet->names[c], alphabet->colours[c]);
      fputs("\n", out);
    }
  }
  // Print core symbols with no complement
  for (i = 1; i <= alphabet->ncore; i++) {
    if (alphabet->complement[i] == 0) {
      print_symbol(out, alphabet->symbols[i], alphabet->names[i], alphabet->colours[i]);
      fprintf(out, "\n");
    }
  }
  // print ambiguous symbols that have comprising characters
  comprise = mm_malloc(sizeof(char) * (alphabet->ncore + 1));
  for (; i <= alphabet->nfull; i++) {
    if (alphabet->ncomprise[i] == 0) break;
    for (j = 0; alphabet->comprise[i][j] != 0; j++) comprise[j] = alphabet->symbols[alphabet->comprise[i][j]];
    comprise[j] = '\0';
    print_symbol(out, alphabet->symbols[i], alphabet->names[i], alphabet->colours[i]);
    fprintf(out, " = %s\n", comprise);
    for (j = 0; alphabet->aliases[i][j] != '\0'; j++) {
      fprintf(out, "%c = %s\n", alphabet->aliases[i][j], comprise);
    }
  }
  free(comprise);
  // Print aliases of core symbols
  for (i = 1; i <= alphabet->ncore; i++) {
    for (j = 0; alphabet->aliases[i][j] != '\0'; j++) {
      fprintf(out, "%c = %c\n", alphabet->aliases[i][j], alphabet->symbols[i]);
    }
  }
  // Print gap symbols
  i = alphabet->nfull;
  if (alphabet->ncomprise[i] == 0) {
    print_symbol(out, alphabet->symbols[i], alphabet->names[i], alphabet->colours[i]);
    fprintf(out, " =\n");
    for (j = 0; alphabet->aliases[i][j] != '\0'; j++) {
      fprintf(out, "%c =\n", alphabet->aliases[i][j]);
    }
  }
}

static const char* alph_like(ALPH_T *alphabet) {
  switch (alphabet->flags & ALPH_FLAG_EXTENDS) {
    case ALPH_FLAG_EXTENDS_RNA: return "rna";
    case ALPH_FLAG_EXTENDS_DNA: return "dna";
    case ALPH_FLAG_EXTENDS_PROTEIN: return "protein";
  }
  return "";
}

/*
 * alph_print_xml
 * print the alphabet to a xml file.
 */
void alph_print_xml(ALPH_T *alphabet, char *tag, char *pad, char *indent, FILE *out) {
  int i, j;
  STR_T *b;
  b = str_create(10);
  fprintf(out, "%s<%s name=\"%s\"", pad, tag, alph_name(alphabet));
  if (alph_extends(alphabet)) fprintf(out, " like=\"%s\"", alph_like(alphabet));
  fprintf(out, ">\n");
  // Print alphabet
  for (i = 0; i < alph_size_full(alphabet); i++) {
    fprintf(out, "%s%s<letter id=\"%s\" symbol=\"%c\"", pad, indent,
        alph_xml_id(alphabet, i, b), alph_char(alphabet, i));
    if (alph_aliases(alphabet, i)[0] != '\0') {
      fprintf(out, " aliases=\"%s\"", alph_aliases(alphabet, i));
    }
    if (alph_ncomprise(alphabet, i) == 1) { // core symbol
      if (alph_has_complement(alphabet)) {
        fprintf(out, " complement=\"%c\"", 
            alph_char(alphabet, alph_complement(alphabet, i)));
      }
    } else { // ambiguous symbol
      str_clear(b);
      for (j = 0; j < alph_ncomprise(alphabet, i); j++) {
        str_appendf(b, "%c", alph_char(alphabet, alph_comprise(alphabet, i, j)));
      }
      fprintf(out, " equals=\"%s\"", str_internal(b));
    }
    if (alph_has_sym_name(alphabet, i)) {
      fprintf(out, " name=\"%s\"", alph_sym_name(alphabet, i));
    }
    if (alph_colour(alphabet, i)) {
      fprintf(out, " colour=\"%06X\"", alph_colour(alphabet, i));
    }
    fprintf(out, "/>\n");
  }
  fprintf(out, "%s</%s>\n", pad, tag);
  str_destroy(b, false);
}

/*
 * alph_print_json
 * print the alphabet to a json writer.
 */
void alph_print_json(ALPH_T *alph, JSONWR_T* jsonwr) {
  int i, j;
  char symbol[2], colour[7];
  char *comprise, *aliases;
  symbol[1] = '\0';
  comprise = mm_malloc(sizeof(char) * (alph_size_core(alph) + 1));
  comprise[alph_size_core(alph)] = '\0';
  jsonwr_start_object_value(jsonwr);
  jsonwr_str_prop(jsonwr, "name", alph_name(alph));
  if (alph_extends(alph)) jsonwr_str_prop(jsonwr, "like", alph_like(alph));
  jsonwr_lng_prop(jsonwr, "ncore", alph_size_core(alph));
  jsonwr_property(jsonwr, "symbols");
  jsonwr_start_array_value(jsonwr);
  for (i = 0; i < alph_size_full(alph); i++) {
    jsonwr_start_object_value(jsonwr);
    symbol[0] = alph_char(alph, i);
    jsonwr_str_prop(jsonwr, "symbol", symbol);
    aliases = alph_aliases(alph, i);
    if (aliases[0] != '\0') jsonwr_str_prop(jsonwr, "aliases", aliases);
    if (alph_has_sym_name(alph, i)) jsonwr_str_prop(jsonwr, "name", alph_sym_name(alph, i));
    if (alph_colour(alph, i) != 0) {
      snprintf(colour, 7, "%06X", alph_colour(alph, i));
      jsonwr_str_prop(jsonwr, "colour", colour);
    }
    if (i < alph_size_core(alph)) {
      if (alph_complement(alph, i) != -1) {
        symbol[0] = alph_char(alph, alph_complement(alph, i));
        jsonwr_str_prop(jsonwr, "complement", symbol);
      }
    } else {
      for (j = 0; j < alph_ncomprise(alph, i); j++) {
        comprise[j] = alph_char(alph, alph_comprise(alph, i, j));
      }
      comprise[j] = '\0';
      jsonwr_str_prop(jsonwr, "equals", comprise);
    }
    jsonwr_end_object_value(jsonwr);
  }
  jsonwr_end_array_value(jsonwr);
  jsonwr_end_object_value(jsonwr);
  free(comprise);
}

/*
 * return the count of complementary pairs
 */
int alph_size_pairs(const ALPH_T *a) {
  int i, c, pairs;
  for (i = 1, pairs = 0; i <= a->ncore; i++) {
    c = a->complement[i];
    if (i < c && a->complement[c] == i) pairs++;
  }
  return pairs;
}

/*
 * alph_equal
 * test if two alphabets are the same.
 */
bool alph_equal(const ALPH_T *a1, const ALPH_T *a2) {
  int i, j;
  // both must be defined or they can not be considered equal
  if (a1 == NULL || a2 == NULL) return false;
  // An alphabet is equal to itself!
  if (a1 == a2) return true;
  // check that the flags are the same while ignoring the built-in flag
  if ((a1->flags & ~ALPH_FLAG_BUILTIN) != (a2->flags & ~ALPH_FLAG_BUILTIN)) return false;
  // check the alphabet name
  if (strcmp(a1->alph_name, a2->alph_name) != 0) return false;
  // now check that all the symbols match exactly
  if (a1->ncore != a2->ncore) return false;
  if (a1->nfull != a2->nfull) return false;
  if (strcmp(a1->symbols, a2->symbols) != 0) return false;
  for (i = 0; i <= a1->nfull; i++) {
    // check the aliases match
    if (strcmp(a1->aliases[i], a2->aliases[i]) != 0) return false;
    // check the names match
    if (strcmp(a1->names[i], a2->names[i]) != 0) return false;
    // check the colours match
    if (a1->colours[i] != a2->colours[i]) return false;
    // check the number of comprising symbols match
    if (a1->ncomprise[i] != a2->ncomprise[i]) return false;
    // check the comprising symbols match
    for (j = 0; j < a1->ncomprise[i]; j++) {
      if (a1->ncomprise[j] != a2->ncomprise[j]) return false;
    }
    // check the complement matches
    if (a1->complement[i] != a2->complement[i]) return false;
  }
  return true;
}

/*
 * alph_core_subset
 * test if the core of one alphabet is the subset of another
 *  0: not subset
 * -1: subset but complements don't match
 *  1: subset with matching complements
 */
int alph_core_subset(const ALPH_T *sub_alph, const ALPH_T *super_alph) {
  bool complement_same;
  int sub_i, super_i;
  uint32_t bitset[4] = {0, 0, 0, 0};
  complement_same = true;
  for (sub_i = 1; sub_i <= sub_alph->ncore; sub_i++) {
    super_i = super_alph->encode2core[(uint8_t)sub_alph->symbols[sub_i]];
    if (super_i < 1) return 0; // missing symbol so not a superset!
    assert(super_i < 128); // there's only 127 ASCII letters
    // check different sub_alph core symbols map to different super_alph core symbols
    if (bitset[(super_i - 1) / 32] & (1 << ((super_i - 1) % 32))) return 0;
    bitset[(super_i - 1) / 32] |= (1 << ((super_i - 1) % 32));
    // check complement
    if (sub_alph->complement[sub_i] && super_alph->complement[super_i]) {
      if (super_alph->encode2core[(uint8_t)sub_alph->symbols[sub_alph->complement[sub_i]]] != super_alph->complement[super_i]) {
        complement_same = false;
      }
    } else if (sub_alph->complement[sub_i] || super_alph->complement[super_i]) {
      complement_same = false;
    }
  }
  return (complement_same ? 1 : -1);
}

/*
 * Test if the symbol is the prime representation.
 * When the alphabet is case insensitive the either case can be considered prime.
 */
bool alph_is_prime(ALPH_T *alph, char letter) {
  int idx;
  char sym;
  idx = alph->encode[(uint8_t)(letter)];
  if (idx == 0) return false;
  sym = alph->symbols[idx];
  if (alph_is_case_insensitive(alph)) {
    return toupper(sym) == toupper(letter);
  } else {
    return sym == letter;
  }
}

/*
 * Get the alphabet string
 */
const char* alph_string(ALPH_T *alph, STR_T *buffer) {
  str_clear(buffer);
  str_append(buffer, alph->symbols+1, alph->ncore);
  return str_internal(buffer);
}


/*
 * Get the number of ambiguous letters in a string.
*/
int get_num_ambiguous_letters(ALPH_T *alph, char *string, int length) {
  int i;
  int count = 0;
  for (i=0; i<length; i++) if (alph_is_ambiguous(alph, string[i])) count++;
  return(count);
}

/*
 *  Tests the letter against the alphabet. If the alphabet is unknown
 *  it attempts to work it out and set it from the letter.
 *  For simplicy this assumes you will pass positions in asscending order.
 *  Returns false if the letter is unacceptable.
 */
bool alph_test(ALPH_T **alpha, int position, char letter) {
  char uc_letter;
  if (*alpha == NULL) {
    uc_letter = toupper(letter);
    switch (position) {
      case 0:
        return (uc_letter == 'A');
      case 1:
        return (uc_letter == 'C');
      case 2:
        if (uc_letter == 'D') {
          *alpha = alph_protein();
          return true;
        }
        return (uc_letter == 'G'); // DNA or RNA
      case 3:
        if (uc_letter == 'T') {
          *alpha = alph_dna();
        } else if (uc_letter == 'U') {
          *alpha = alph_rna();
        } else {
          return false;
        }
        return true;
      default:// Bad state!
        die("Should not still be attempting to guess by the 5th letter "
            "(position = %d).", position);
        return false;
    }
  } else {
    position++; // adjust to alphabet indexes which use 0 as a error
    if (position > (*alpha)->ncore) return false; // position too big
    return ((*alpha)->encode2core[(uint8_t)letter] == position);
  }
}

/*
 *  Checks the symbols against the core alphabet.
 *  Returns true if it matches.
 */
bool alph_check(ALPH_T *alph, char *syms) {
  int i, alen;
  // test the alphabet
  alen = strlen(syms);
  if (alph_size_core(alph) != alen) return false;
  for (i = 0; i < alen; i++) {
    if (alph_index(alph, syms[i]) != i) return false;
  }
  return true;
}

/*
 *  Tests the alphabet string and attempts to return the associated
 *  built-in alphabet.
 *  If the alphabet string is from some buffer a max size can be set
 *  however it will still only test until the first null byte.
 *  If the string is null terminated just set a max > 20 
 */
ALPH_T* alph_type(const char *alphabet, int max) {
  int i;
  ALPH_T *alph;
  alph = NULL;
  for (i = 0; i < max && alphabet[i] != '\0'; ++i) {
    if (!alph_test(&alph, i, alphabet[i])) {
      if (alph != NULL) alph_release(alph);
      return NULL;
    }
  }
  if (alph == NULL) return NULL;
  if (i != alph->ncore) {
    alph_release(alph);
    return NULL;
  }
  return alph;
}

/*
 * resize_markov_model
 * Allow extending a markov model to allow for a larger alphabet - for example
 * to add ambiguous symbols which can be derived from the existing ones.
 */
void resize_markov_model(int asize0, int asize1, ARRAY_T *tuples, int *order_p) {
  int i, j, ntuples, order, *indexes, len, index0, index1, last;
  assert(asize0 > 0 && asize0 < asize1);
  // given that we know the alphabet size intially we can work out the order of the model
  i = 0;
  ntuples = 0;
  while (ntuples < get_array_length(tuples)) {
    ntuples += pow(asize0, ++i);
  }
  if (ntuples != get_array_length(tuples)) die("Markov model resize failed due to incorrect specified initial alphabet size");
  order = i - 1;
  // calculate the new tuple count needed to fit the new alphabet size
  for (i = 0, ntuples = 0; i <= order; i++) ntuples += pow(asize1, i + 1);
  // resize the array
  resize_array(tuples, ntuples);
  // now start moving entries working in reverse order so collisions don't occur
  indexes = mm_malloc(sizeof(int) * (order + 1));
  last = get_array_length(tuples);
  for (len = order + 1; len > 0; len--) { // length
    // set indexes to last value in alphabet
    for (i = 0; i < len; i++) indexes[i] = asize0;
    i = len - 1;
    while (i >= 0) {
      // convert indexes into positions in the old and new background
      index0 = indexes[0];
      index1 = indexes[0];
      for (j = 1; j < len; j++) {
        index0 = (index0 * asize0) + indexes[j];
        index1 = (index1 * asize1) + indexes[j];
      }
      index0 -= 1;
      index1 -= 1;
      // zero any missed tuples
      for (j = last - 1; j > index1; j--) {
        set_array_item(j, 0.0, tuples);
      }
      // set the moved item
      set_array_item(index1, get_array_item(index0, tuples), tuples);
      // update the last set index
      last = index1;
      // deincrement
      for (i = len - 1; i >= 0; i--) {
        indexes[i]--;
        if (indexes[i] > 0) {
          break;
        } else {
          indexes[i] = asize0;
        }
      }
    }
  }
  free(indexes);
  if (order_p != NULL) *order_p = order;
}

/*
 * Given a markov model with frequencies for the core letters, calculate the
 * values of the ambiguous letters from sums of the core letter values.
 * As an alternative to having all the ambiguous values it can be limited
 * to the wildcard only which is assumed to be the first symbol after the core
 * symbols.
 */
void extend_markov_model(ALPH_T* alph, bool wildcard_only, AMBIG_CALC_EN method, ARRAY_T* tuples) {
  int order, asize, asize_gapless, index1, index2, i, j, k, len, parts;
  int *indexes, *pos;
  double value, sum;
  bool first, has_gap_sym;
  asize = (wildcard_only ? alph_size_wild(alph) : alph_size_full(alph));
  // resize the array to fit the new tuples
  resize_markov_model(alph_size_core(alph), asize, tuples, &order);
  // test to see if need to skip gap symbols
  has_gap_sym = (alph_ncomprise(alph, asize - 1) == 0);
  asize_gapless = (has_gap_sym ? asize - 1 : asize);
  // now iterate over all 
  indexes = mm_malloc(sizeof(int) * (order + 1));
  pos = mm_malloc(sizeof(int) * (order + 1));
  for (len = 1; len <= (order + 1); len++) { // length
    // set indexes to last value in alphabet, that is not a gap
    for (i = 0; i < len; i++) indexes[i] = asize_gapless;
    i = len;
    while (i >= 0) {
      // skip any index that doesn't include ambiguous characters
      for (j = 0; j < len; j++) {
        if (indexes[j] > alph_size_core(alph)) break;
      }
      if (j == len) goto next_index;
      // convert index into position
      index1 = indexes[0];
      for (j = 1; j < len; j++) {
        // Note that indexes is effectively alph_index(alph, sym) + 1
        // hence why we aren't adding one here. It is never zero!
        index1 = (index1 * asize) + indexes[j];
      }
      index1 -= 1;
      // init value calcuation
      first = true;
      sum = 0.0;
      parts = 1; // don't want to divide by zero!
      // determine indexes of all comprising
      for (j = 0; j < len; j++) pos[j] = 0;
      j = len - 1;
      while (j >= 0) {
        // calculate the index substituting comprising characters
        index2 = alph_comprise(alph, indexes[0] - 1, pos[0]) + 1;
        for (k = 1; k < len; k++) {
          index2 = (index2 * asize) + alph_comprise(alph, indexes[k] - 1, pos[k]) + 1;
        }
        index2 -= 1;
        value = get_array_item(index2, tuples);
        // now sum values
        if (first) {
          sum = value;
          first = false;
        } else {
          switch(method) {
            case SUM_LOGS:
              sum = LOG_SUM(sum, value);
              break;
            case AVG_FREQS:
            case SUM_FREQS:
              sum += value;
              break;
          }
          parts++;
        }
        // increment
        for (j = len - 1; j >= 0; j--) {
          pos[j]++;
          if (pos[j] < alph_ncomprise(alph, indexes[j] - 1)) {
            break;
          } else {
            pos[j] = 0;
          }
        }
      }
      // set value
      if (method == AVG_FREQS) {
        set_array_item(index1, sum / parts, tuples);
      } else {
        set_array_item(index1, sum, tuples);
      }
next_index:
      // deincrement
      for (i = len - 1; i >= 0; i--) {
        indexes[i]--;
        if (indexes[i] > 0) {
          break;
        } else {
          indexes[i] = asize_gapless;
        }
      }
    }
  }
  // cleanup
  free(indexes);
  free(pos);
}

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
void extrapolate_markov_model(int asize0, int asize1, double ambig_fraction, ARRAY_T* tuples) {
  int order, index, index_word, i, j, len, range_start, range_length;
  int *indexes;
  double value;
  bool has_gap_sym;
  // 1) resize the model
  // asize1 could be alph_size_wild, alph_size_full or alph_size_all (no gap)
  assert(asize1 > asize0); // must include at least one ambiguous symbol
  // resize the array to fit the new tuples
  resize_markov_model(asize0, asize1, tuples, &order);

  // 2) generate ambigous 0-order frequencies as (ambig_fraction / ambig_count)
  value = ambig_fraction / (asize1 - asize0);
  for (i = asize0; i < asize1; i++) {
    set_array_item(i, value, tuples);
  }
  // normalise
  range_start = 0;
  range_length = asize1;
  normalize_subarray(range_start, range_length, 0.0, tuples);

  // 3) generate n-order ambig frequencies Pr(word . ambig) = Pr(word) * Pr(ambig)
  // now iterate over all 
  index_word = 0;
  indexes = mm_malloc(sizeof(int) * (order + 1));
  for (len = 2; len <= (order + 1); len++) { // length
    range_start += range_length;
    range_length *= asize1;
    // set indexes to last value in alphabet
    for (i = 0; i < len; i++) indexes[i] = asize1;
    i = len;
    while (i >= 0) {
      // skip any index that doesn't include ambiguous characters
      for (j = 0; j < len; j++) {
        if (indexes[j] > asize0) break;
      }
      if (j == len) goto next_index;
      // convert index into position
      index = indexes[0];
      for (j = 1; j < len; j++) {
        index_word = index;
        // Note that indexes is effectively alph_index(alph, sym) + 1
        // hence why we aren't adding one here. It is never zero!
        index = (index * asize1) + indexes[j];
      }
      index -= 1;
      index_word -= 1;
      // calculate value
      value = get_array_item(index_word, tuples) * get_array_item(indexes[len-1] - 1, tuples);
      set_array_item(index, value, tuples);
next_index:
      // deincrement
      for (i = len - 1; i >= 0; i--) {
        indexes[i]--;
        if (indexes[i] > 0) {
          break;
        } else {
          indexes[i] = asize1;
        }
      }
    }
    // normalise
    normalize_subarray(range_start, range_length, 0.0, tuples);
  }
  // cleanup
  free(indexes);
}

/*
 * Calculate the values of the ambiguous letters from sums of
 * the normal letter values (0-order only).
 * Assumes the array is already large enough for the ambiguous values.
 */
void calc_ambigs(ALPH_T* alph, bool log_space, ARRAY_T* array) {
  int i;
  uint8_t *c;
  PROB_T sum;
  PROB_T value;
  if (alph == NULL) die("Alphabet uninitialized.\n");
  if (array == NULL) die("Array unitialized.\n");
  if (get_array_length(array) < alph->nfull) resize_array(array, alph->nfull);
  for (i = alph->ncore + 1; i <= alph->nfull; i++) {
    sum = 0.0;
    for (c = alph->comprise[i]; *c != 0; c++) {
      value = get_array_item(*c - 1, array);
      if (log_space) {
        sum = LOG_SUM(sum, value);
      } else {
        sum += value;
      }
    }
    set_array_item(i - 1, sum, array);
  }
}


/*
 * Replace the elements an array of frequences with the average
 * over complementary bases.
 * This only handles two-way complementary symbols, ignores one-way.
 */
void average_freq_with_complement(ALPH_T *alph, ARRAY_T *freqs) {
  int i, c;
  PROB_T avg_freq;
  for (i = 1; i <= alph->ncore; i++) {
    c = alph->complement[i];
    if (alph->complement[c] == i && i < c) { // two-way complementary pair
      avg_freq = (get_array_item(i - 1, freqs) + get_array_item(c - 1, freqs)) / 2.0;
      set_array_item(i - 1, avg_freq, freqs);
      set_array_item(c - 1, avg_freq, freqs);
    }
  }
}

/*
 * Averages background probability tuples with their reverse complement.
 */
static void average_rc(ALPH_T *alph, int order, int len, int *tuple, ARRAY_T *bg) {
  int i, j, ti, rci;
  double avg;
  for (i = 1; i <= alph->ncore; i++) {
    tuple[len] = i;
    // calculate tuple and reverse complement tuple index
    for (ti = 0, rci = 0, j = 0; j <= len; j++) {
      ti = (alph->ncore * ti) + tuple[j];
      rci = (alph->ncore * rci) + alph->complement[tuple[len - j]];
    }
    // average pair if it is the first time we've seen it
    if (ti < rci) {
      avg = (get_array_item(ti - 1, bg) + get_array_item(rci - 1, bg)) / 2.0;
      set_array_item(ti - 1, avg, bg);
      set_array_item(rci - 1, avg, bg);
    }
    // recursively average larger tuples
    if (len < order) average_rc(alph, order, len + 1, tuple, bg);
  }
}

/*
 * Averages background probability tuples with their reverse complement.
 * Expects only the core alphabet - this must be applied before extending.
 */
void average_rc_markov_model(ALPH_T *alph, int order, ARRAY_T *bg) {
  int *tuple;
  tuple = mm_malloc(sizeof(int) * (order + 1));
  average_rc(alph, order, 0, tuple, bg);
  free(tuple);
}

/*
 * Load the frequencies that used to be in hash_alph.h
 *
 */
ARRAY_T* get_mast_frequencies(ALPH_T *alph, bool has_ambigs, bool translate) {
  const char* MAST_DNA_SYMS = "ABCDGHKMNRSTVWY";
  PROB_T const MAST_DNA_FREQS[] = {
    0.281475655 /* A */, 0.000000649 /* B */, 0.221785822 /* C */,
    0.000001389 /* D */, 0.228634607 /* G */, 0.000001612 /* H */,
    0.000006323 /* K */, 0.000005848 /* M */, 0.000991686 /* N */,
    0.000015904 /* R */, 0.000009514 /* S */, 0.267048106 /* T */,
    0.000000968 /* V */, 0.000005846 /* W */, 0.000016070 /* Y */
  };
  const char* MAST_PROT_SYMS = "ABCDEFGHIKLMNPQRSTVWXYZ";
  PROB_T const MAST_PROT_FREQS[] = {
    0.073091885 /* A */, 0.000021047 /* B */, 0.018145453 /* C */,
    0.051687956 /* D */, 0.062278511 /* E */, 0.040243411 /* F */,
    0.069259642 /* G */, 0.022405456 /* H */, 0.056227000 /* I */,
    0.058435042 /* K */, 0.091621836 /* L */, 0.023044274 /* M */,
    0.046032137 /* N */, 0.050623807 /* P */, 0.040715284 /* Q */,
    0.051846246 /* R */, 0.073729031 /* S */, 0.059352333 /* T */,
    0.064298546 /* V */, 0.013328158 /* W */, 0.000941980 /* X */,
    0.032649745 /* Y */, 0.000021111 /* Z */
  };
  PROB_T const MAST_XPROT_FREQS[] = {
    0.06995463 /* A */, 0.00000228 /* B */, 0.02398853 /* C */,
    0.03950999 /* D */, 0.05199564 /* E */, 0.03955549 /* F */,
    0.07363406 /* G */, 0.02911862 /* H */, 0.04092760 /* I */,
    0.05187635 /* K */, 0.09566732 /* L */, 0.01747826 /* M */,
    0.03216447 /* N */, 0.06506016 /* P */, 0.04180057 /* Q */,
    0.06748940 /* R */, 0.08353402 /* S */, 0.05231376 /* T */,
    0.05927076 /* V */, 0.01585647 /* W */, 0.02720962 /* X */,
    0.02158787 /* Y */, 0.00000413 /* Z */
  };
  
  // vars
  const PROB_T *values;
  const char* syms;
  ARRAY_T *freqs;
  PROB_T value;
  int i, index;
  // allocate the array if it is not given
  freqs = allocate_array(has_ambigs ? alph_size_full(alph) : alph_size_core(alph));
  init_array(0, freqs);
  // get the frequencies
  if (alph_is_builtin_dna(alph)) {
    values = MAST_DNA_FREQS;
    syms = MAST_DNA_SYMS;
  } else if (alph_is_builtin_protein(alph)) {
    values = (translate ? MAST_XPROT_FREQS : MAST_PROT_FREQS);
    syms = MAST_PROT_SYMS;
  } else {
    // TLB 1-Feb-2017 changed default when not standard alphabet to MODEL frequencies
    //value = 1.0 / alph_size_core(alph);
    //for (i = 0; i < alph_size_core(alph); i++) {
    //  set_array_item(i, value, freqs);
    //}
    //return freqs;
    myfree(freqs);
    return NULL;
  }
  for (i = 0; syms[i] != '\0'; i++) {
    index = alph_index(alph, syms[i]);
    assert(index >= 0);
    if (has_ambigs || index < alph_size_core(alph)) {
      set_array_item(index, values[i], freqs);
    }
  }
  return freqs;
}

/*
 * Load the non-redundant database frequencies into the array.
 */
ARRAY_T* get_nrdb_frequencies(ALPH_T *alph, ARRAY_T *freqs) {
  // constants
  PROB_T const NRDB_AA[] = {
    0.073164, /* A */ 0.018163, /* C */ 0.051739, /* D */ 0.062340, /* E */ 
    0.040283, /* F */ 0.069328, /* G */ 0.022428, /* H */ 0.056282, /* I */
    0.058493, /* K */ 0.091712, /* L */ 0.023067, /* M */ 0.046077, /* N */ 
    0.050674, /* P */ 0.040755, /* Q */ 0.051897, /* R */ 0.073802, /* S */ 
    0.059411, /* T */ 0.064362, /* V */ 0.013341, /* W */ 0.032682  /* Y */
  };
  PROB_T const NRDB_DNA[] = {
    0.281774, /* A */ 0.222020, /* C */ 0.228876, /* G */ 0.267330 /* T */
  };
  // vars
  const PROB_T *nrdb_freqs;
  int i;
  // FIXME this is not a great implementation but unless we put the frequencies
  // into the alphabets I'm not sure how to solve it...
  if (alph_is_builtin_dna(alph)) {
    nrdb_freqs = NRDB_DNA;
  } else if (alph_is_builtin_protein(alph)) {
    nrdb_freqs = NRDB_AA;
  } else {
    // well we could make a guess at what alphabet it is and try to use
    // NRDB for dna or protein but I'm not very happy with that option either
    // return get_uniform_frequencies(alph, freqs);
    // TLB 3-Feb-2017; Changed to return NULL so model frequencies will be used for
    // non-standard alphabets.
    return(NULL);
  }
  // allocate the array if it is not given
  if (freqs == NULL) freqs = allocate_array(alph->ncore);
  else if (get_array_length(freqs) < alph->ncore) resize_array(freqs, alph->ncore);
  // copy the frequencies into the array 
  for (i = 0; i < alph->ncore; ++i) {
    set_array_item(i, nrdb_freqs[i], freqs);
  }
  normalize_subarray(0, alph->ncore, 0.0, freqs);
  return freqs;
}

/*
 * Load uniform frequencies into the array.
 */
ARRAY_T* get_uniform_frequencies(ALPH_T *alph, ARRAY_T *freqs) {
  int i;
  PROB_T freq;
  // calculate the uniform frequency
  freq = 1.0/alph->ncore;
  // ensure the array is allocated
  if (freqs == NULL) freqs = allocate_array(alph->ncore);
  else if (get_array_length(freqs) < alph->ncore) resize_array(freqs, alph->ncore);
  // assign core symbols the uniform frequency
  for (i = 0; i < alph->ncore; i++) { 
    set_array_item(i, freq, freqs); 
  }
  return freqs;
}

/*
 * Add pseudocount and renormalize a frequency array
*/
void normalize_frequencies(ALPH_T *alph, ARRAY_T *freqs, double pseudo) {
  int i;
  // Get the sum of the frequencies.
  PROB_T sum = pseudo * alph->ncore;	// sum plus pseudocounts
  for (i = 0; i < alph->ncore; i++) { 
    sum += get_array_item(i, freqs);
  }
  // Add pseudocount and normalize.
  for (i = 0; i < alph->ncore; i++) { 
    PROB_T freq = (get_array_item(i, freqs) + pseudo) / sum;
    set_array_item(i, freq, freqs);
  }
}

/*
 * Load file frequencies into the array.
 */
ARRAY_T* get_file_frequencies(ALPH_T *alph, char *bg_filename) {
  ARRAY_T *freqs;
  int order = 0;
  freqs = load_markov_model(alph, &order, bg_filename);
  return freqs;
}


/*
 * Get a background distribution
 *  - by reading values from a file if filename is given, or
 *  - equal to the NRDB frequencies if filename is NULL.
 */
ARRAY_T* get_background(ALPH_T *alph, char* bg_filename) {
  ARRAY_T* background;

  if ((bg_filename == NULL) || (strcmp(bg_filename, "nrdb") == 0)) {
    background = get_nrdb_frequencies(alph, NULL);
  } else {
    background = get_file_frequencies(alph, bg_filename);
  }
  calc_ambigs(alph, false, background);

  return(background);
}

/*
 * Create a list of symbols from an index.
 */
static inline __attribute__((always_inline)) char* index2chain(const char *alphabet, int index) {
  static char *chain = NULL;
  int i, alen, len;
  char *sym;
  // determine the alphabet size
  alen = strlen(alphabet);
  // count the length
  len = 0;
  for (i = index + 1; i > 0; i = (i - 1) / alen) len++;
  // allocate space for the chain
  chain = mm_realloc(chain, sizeof(char) * (len + 1));
  // fill in the chain
  for (i = index + 1, sym = chain+(len - 1); i > 0; i = ((i - 1) / alen), sym--) {
    *sym = alphabet[(i - 1) % alen];
  }
  chain[len] = '\0';
  return chain;
}

/*
 * Constants and macros for load_markov_model_entry.
 */
typedef enum {
  MKM_OK, MKM_EOF,
  MKM_SYMS_LONG, MKM_SYMS_ERR,
  MKM_NUM_MISSING, MKM_NUM_LONG, MKM_NUM_ERR, MKM_NUM_ZERO,
  MKM_JUNK} MKM_RESULT;

#define MKM_LOAD_MORE(eof_result) \
    size = fread(buffer, sizeof(char), sizeof(buffer) / sizeof(char), fp); \
    text = buffer; \
    offset = 0; \
    if (size == 0) { \
      str_clear(storage); \
      return eof_result; \
    }

/*
 * Load the next entry from the markov model.
 * This allows any amount of whitespace and comments of any length.
 * Both symbols and number have a specifiable maximum buffer size.
 * Symbols are restricted to visible ASCII.
 * Numbers must be valid and satisify 0 <= value <= 1.
 * There must be no trailing junk though comments are allowed at line end.
 *
 * Invariant:
 * The file pointer is at the start of a line or the storage contains any
 * characters needed for the start of the line.
 */
static MKM_RESULT load_markov_model_entry(FILE *fp, STR_T *storage,
    size_t max_syms_len, STR_T *syms, size_t max_num_len, STR_T *num, double *value) {
  char buffer[100];
  char *text, *endptr;
  size_t offset, size;
  MKM_RESULT status;
  // reset the outputs
  str_clear(syms);
  str_clear(num);
  *value = -1;
  status = MKM_OK;
  // initially use the storage as the text source
  size = str_len(storage);
  text = str_internal(storage);
  offset = 0;
skip_whitespace:
  // skip over whitespace
  while (true) {
    for (; offset < size; offset++) {
      if (!isspace(text[offset])) {
        if (text[offset] == '#') {
          goto skip_comment;
        } else {
          goto read_syms;
        }
      }
    }
    MKM_LOAD_MORE(MKM_EOF)
  }
skip_comment:
  // skip the comment
  while (true) {
    for (; offset < size; offset++) {
      if (text[offset] == '\n' || text[offset] == '\r') {
        goto skip_whitespace; // this is above
      }
    }
    MKM_LOAD_MORE(MKM_EOF)
  }
read_syms:
  // read symbols until whitespace, too long or an invalid symbol
  while (true) {
    for (; offset < size; offset++) {
      if (isspace(text[offset])) {
        // found the end of the symbols
        goto skip_gap;
      } else if (str_len(syms) >= max_syms_len) {
        // the symbols go on for too long
        status = MKM_SYMS_LONG;
        goto skip_to_eol;
      } else {
        str_append(syms, text+offset, 1);
        if (!(text[offset] >= '!' && text[offset] <= '~')) {
          // Illegal symbol!
          status = MKM_SYMS_ERR;
          goto skip_to_eol;
        }
      }
    }
    MKM_LOAD_MORE(MKM_NUM_MISSING)
  }
skip_gap:
  // skip over whitespace
  while (true) {
    for (; offset < size; offset++) {
      if (!isspace(text[offset])) {
        goto read_num;
      } else if (text[offset] == '\n' || text[offset] == '\r') {
        // number left unspecified!
        status = MKM_NUM_MISSING;
        goto skip_to_eol;
      }
    }
    MKM_LOAD_MORE(MKM_NUM_MISSING)
  }
read_num:
  // read number until whitespace, too long or invalid number
  while (true) {
    for (; offset < size; offset++) {
      if (isspace(text[offset])) {
        goto parse_num;
      } else if (str_len(num) >= max_num_len) {
        // the symbols go on for too long
        status = MKM_NUM_LONG;
        goto skip_to_eol;
      } else {
        str_append(num, text+offset, 1);
        if (!((text[offset] >= '0' && text[offset] <= '9') ||
              text[offset] == '.' ||
              text[offset] == 'e' || text[offset] == 'E' ||
              text[offset] == '+' || text[offset] == '-')) {
          // Illegal number symbol!
          status = MKM_NUM_ERR;
          goto skip_to_eol;
        }
      }
    }
    // load more
    size = fread(buffer, sizeof(char), sizeof(buffer) / sizeof(char), fp);
    text = buffer;
    offset = 0;
    if (size == 0) {
      str_clear(storage);
      if (str_len(num) > 0) {
        goto parse_num;
      } else {
        return MKM_NUM_MISSING;
      }
    }
  }
parse_num:
  // found the end of the number so parse it
  *value = strtod(str_internal(num), &endptr);
  if (*endptr != '\0' || *value < 0 || *value > 1) {
    status = MKM_NUM_ERR;
    goto skip_to_eol;
  } else if (*value == 0) {
    status = MKM_NUM_ZERO;
    goto skip_to_eol;
  }
  // check that we don't have junk before the end of line
  while (true) {
    for (; offset < size; offset++) {
      if (text[offset] == '#' || text[offset] == '\n' || text[offset] == '\r') {
        // found comment or eol
        goto skip_to_eol;
      } else if (!isspace(text[offset])) {
        // found junk
        status = MKM_JUNK;
        goto skip_to_eol;
      }
    }
    MKM_LOAD_MORE(MKM_OK)
  }
skip_to_eol:
  // skip until the end of the line and return the status
  while (true) {
    for (; offset < size; offset++) {
      if (text[offset] == '\n' || text[offset] == '\r') {
        // store the unprocessed text in storage
        // note that we have to be careful because text could point to storage!
        if (text == buffer) {
          str_set(storage, text+offset, size - offset);
        } else {
          str_delete(storage, 0, offset);
        }
        return status;
      }
    }
    MKM_LOAD_MORE(status)
  }
}

/*
 * Load the next entry from the markov model.
 * Exits if any errors occur.
 * Returns true on success and false on EOF.
 */
static bool load_markov_model_entry2(const char *bg_filename, FILE *fp, STR_T *storage,
    size_t max_syms_len, STR_T *syms, size_t max_num_len, STR_T *num, double *value) {
  switch (load_markov_model_entry(fp, storage, max_syms_len, syms, max_num_len, num, value)) {
      case MKM_OK:
        return true;
      case MKM_SYMS_LONG:
        die("Background file \"%s\" contained a symbol chain that "
            "was excessively long.", bg_filename);
      case MKM_SYMS_ERR:
        die("Background file \"%s\" listed a symbol that is not in "
            "the visible ASCII range.", bg_filename);
      case MKM_NUM_MISSING:
        die("Background file \"%s\" did not list a probability on "
            "the line listing the symbol chain \"%s\".", 
            bg_filename, str_internal(syms));
      case MKM_NUM_LONG:
        die("Background file \"%s\" has a probability string that "
            "was excessively long for the symbol chain \"%s\".",
            bg_filename, str_internal(syms));
      case MKM_NUM_ERR:
        die("Background file \"%s\" has a probability string "
            "that could not be parsed as a probability for "
            "the symbol chain \"%s\".",
            bg_filename, str_internal(num), str_internal(syms));
      case MKM_NUM_ZERO:
        die("Background file \"%s\" has a probability of zero "
            "for the symbol chain \"%s\" which is not supported. "
            "Regenerate your background using pseudocounts to avoid "
            "this problem.", bg_filename, str_internal(syms));
      case MKM_JUNK:
        die("Background file \"%s\" has junk text following the "
            "entry for the symbol chain \"%s\".",
            bg_filename, str_internal(syms));
      case MKM_EOF:
      default:
        return false;
  }
}

/*
 * Load background file frequencies into the array.
 */
ARRAY_T* load_markov_model_without_alph(const char *bg_filename, int *order, char **syms) {
  FILE *fp;
  STR_T *sym_buf, *num_buf, *line_buf;
  double value;
  RBTREE_T *symbols;
  RBNODE_T *node;
  ARRAY_T *tuples;
  int i, alen, index, order_max, order_count;
  uint8_t lookup[UCHAR_MAX + 1];
  char *symbol;
  // get the order value
  order_max = (order != NULL ? *order: INT_MAX);
  // setup buffers
  sym_buf = str_create(10);
  num_buf = str_create(50);
  line_buf = str_create(100);
  // open the background file
  if (!(fp = fopen(bg_filename, "r"))) {
    die("Unable to open background file \"%s\" for reading.\n", bg_filename);
  }
  // setup symbol map for order 0 bg
  symbols = rbtree_create(alph_str_cmp, rbtree_strcpy, free, rbtree_dblcpy, free); // TODO FIXME use true alphabet symbol order
  // read each entry of length 1
  while (load_markov_model_entry2(bg_filename, fp, line_buf, 100, sym_buf, 100, num_buf, &value)) {
    if (str_len(sym_buf) > 1) break;
    if (!rbtree_make(symbols, str_internal(sym_buf), &value)) {
      die("Background file \"%s\" contains more than one entry "
          "for symbol chain \"%s\".", bg_filename, str_internal(sym_buf));
    }
  }
  // rearrange the entries into normal alphabet order and setup the lookup
  memset(lookup, 0, sizeof(lookup));
  alen = rbtree_size(symbols);
  tuples = allocate_array(alen);
  *syms = mm_malloc(sizeof(char) * (alen + 1));
  for (index = 0, node = rbtree_first(symbols); node != NULL; index++, node = rbtree_next(node)) {
    char sym = ((char*)rbtree_key(node))[0];
    (*syms)[index] = sym;
    lookup[(uint8_t)sym] = index + 1;
    set_array_item(index, *((double*)rbtree_value(node)), tuples);
  }
  (*syms)[alen] = '\0';
  // set the order to zero so we can count it
  if (order != NULL) *order = 0;
  // clean up the unneeded order 0 symbol map
  rbtree_destroy(symbols);
  // process other entries
  order_count = alen;
  while (str_len(sym_buf) > 0 && (str_len(sym_buf) - 1) <= order_max) {
    // convert the symbols into an index
    index = 0;
    for (symbol = str_internal(sym_buf); *symbol; symbol++) {
      int symbol_index = lookup[(uint8_t)*symbol];
      if (symbol_index == 0) {
        die("Background file \"%s\" lists the symbol '%c' which is not "
            "a previously defined symbol.", bg_filename, *symbol);
      }
      index = (alen * index) + symbol_index;
    }
    index -= 1;
    // check if we need to expand the tuples to set the chain probability
    if (index >= get_array_length(tuples)) {
      // scan the last order to make sure it was all set
      for (i = get_array_length(tuples) - order_count; i < get_array_length(tuples); i++) {
        if (get_array_item(i, tuples) == -1) {
          die("Background file \"%s\" does not list a probability for the "
              "symbol chain \"%s\".", bg_filename,  index2chain(*syms, i));
        }
      }
      // calculate how much we must expand the tuples by
      order_count *= alen;
      // check that the new item will fit in the expanded array
      if (index >= (get_array_length(tuples) + order_count)) {
        die("Background file \"%s\" does not list all shorter chains "
            "before \"%s\".", bg_filename, str_internal(sym_buf));
      }
      // expand the array
      resize_init_array(tuples, get_array_length(tuples) + order_count, -1);
      // increment the order
      if (order != NULL) *order += 1;
    } else if (get_array_item(index, tuples) != -1) {
      die("Background file \"%s\" has a repeated definition of the "
          "symbol chain \"%s\".", bg_filename, str_internal(sym_buf));
    }
    // store the probability
    set_array_item(index, value, tuples);
    // read the next value
    load_markov_model_entry2(bg_filename, fp, line_buf, 100, sym_buf, 100, num_buf, &value);
  }
  // scan the last order to make sure it was all set
  for (i = get_array_length(tuples) - order_count; i < get_array_length(tuples); i++) {
    if (get_array_item(i, tuples) == -1) {
      die("Background file \"%s\" does not list a probability for the "
          "symbol chain \"%s\".", bg_filename,  index2chain(*syms, i));
    }
  }
  // clean up
  str_destroy(line_buf, false);
  str_destroy(sym_buf, false);
  str_destroy(num_buf, false);
  // return the background
  return tuples;
}

/*
 * Load background file frequencies into the array.
 */
ARRAY_T* load_markov_model(ALPH_T *alph, int* order, const char *bg_filename) {
  char *syms;
  ARRAY_T *tuples;
  syms = NULL;
  tuples = load_markov_model_without_alph(bg_filename, order, &syms);
  if (alph == NULL || alph_check(alph, syms)) {
    free(syms);
    return tuples;
  } else {
    die("Background file \"%s\" is not the %s alphabet as it contained the symbols \"%s\".", bg_filename, alph_name(alph), syms);
    return NULL; // make the compiler happy.
  }
}

/*
 * Structure that holds information for calculating a markov model from sequence data.
 */
struct bgcalc {
  ARRAY_T *tuples;
  double *totals;
  int *history;
};

/*
 * Setup the structure used to store the ongoing background calculation.
 */
static inline BGCALC_T* bgcalc_setup(ALPH_T *alph, int order, double epsilon) {
  BGCALC_T *calc;
  int ntuples, i;
  // calculate how many different counts we need for this order model
  for (i = 0, ntuples = 0; i <= order; i++) ntuples += pow(alph->ncore, i+1);
  // allocate the structure
  calc = mm_malloc(sizeof(BGCALC_T));
  calc->tuples = allocate_array(ntuples);
  // initialise the counts to the pseudocount epsilon
  init_array(epsilon, calc->tuples);
  calc->totals = mm_malloc(sizeof(double) * (order + 1));
  for (i = 0; i <= order; i++) calc->totals[i] = epsilon * pow(alph->ncore, i+1);
  calc->history = mm_malloc(sizeof(int) * (order + 1));
  return calc;
}

/*
 * Finish the calculation of the background probabilities.
 */
static inline ARRAY_T* bgcalc_finish(ALPH_T *alph, int order, BGCALC_T *calc) {
  ARRAY_T *result;
  int index, ncore_pow_w, i, j;
  // convert counts to probabilities
  index = 0;
  for (i = 0; i <= order; i++) {
    ncore_pow_w = pow(alph->ncore, i + 1);
    for (j = 0; j < ncore_pow_w; j++, index++) {
      set_array_item(index, get_array_item(index, calc->tuples) / calc->totals[i], calc->tuples);
    }
  }
  result = calc->tuples;
  // free background calculator
  free(calc->totals);
  free(calc->history);
  memset(calc, 0, sizeof(BGCALC_T));
  free(calc);
  return result;
}

/*
 * Clear the history of the bgcalc object.
 *
 * According to GCC's manual this should be equalivent to a macro in speed and
 * the always_inline attribute should make it inline even when not optimising.
 */
static inline __attribute__((always_inline)) void bgcalc_clear(int order, BGCALC_T *calc) {
  int i;
  for (i = 0; i <= order; i++) calc->history[i] = 0;
}

/*
 * Update the history of the bgcalc object.
 *
 * According to GCC's manual this should be equalivent to a macro in speed and
 * the always_inline attribute should make it inline even when not optimising.
 */
static inline __attribute__((always_inline)) void bgcalc_update(ALPH_T *alph, int order, BGCALC_T *calc, int index) {
  int i;
  if (index) {
    // update history
    for (i = order; i; i--) {
      if (!calc->history[i - 1]) continue;
      calc->history[i] = (alph->ncore * calc->history[i-1]) + index;
      incr_array_item(calc->history[i] - 1, 1, calc->tuples);
      calc->totals[i] += 1;
    }
    calc->history[0] = index;
    incr_array_item(index - 1, 1, calc->tuples);
    calc->totals[0] += 1;
  } else {
    // ambiguous or unknown! clear history
    bgcalc_clear(order, calc);
  }
}

/*
 * Combines the features of get_markov_from_sequence and calculate_background.
 *
 * When *save is NULL then it will be initialised.
 *
 * When seq is specified then the counts for the contained chains will be added
 * to the calculation.
 *
 * If join_seq is true then there will be considered to be no gap between the last
 * passed sequence and this one.
 *
 * When seq is not specified then the result of the background calculation will
 * be returned in the array and the save object will be deallocated.
 *
 */
ARRAY_T* calculate_markov_model(ALPH_T *alph, int order,
    double epsilon, bool join_seq, const char *seq, BGCALC_T** save) {
  BGCALC_T *calc;
  const char *s;
  if (*save == NULL) *save = bgcalc_setup(alph, order, epsilon);
  calc = *save;
  if (seq != NULL) {
    if (!join_seq) bgcalc_clear(order, calc);
    for (s = seq; *s; s++) bgcalc_update(alph, order, calc,
        alph->encode2core[(uint8_t)*s]);
    return NULL;
  } else {
    *save = NULL;
    return bgcalc_finish(alph, order, calc);
  }
}

/*
 * Take the counts from each ambiguous character and evenly distribute
 * them among the corresponding concrete characters.
 *
 * This function operates in log space.
 */
void dist_ambigs(ALPH_T *alph, ARRAY_T* freqs) {
  int i;
  uint8_t *c;
  PROB_T ambig_count, core_count;
  if (alph == NULL) die("Alphabet uninitialized.\n");
  for (i = alph->ncore + 1; i <= alph->nfull; i++) {
    // Get the count to be distributed
    ambig_count = get_array_item(i - 1, freqs);
    // Divide it by the number of comprising core symbols
    ambig_count -= my_log2((PROB_T)(alph->ncomprise[i]));
    // Distribute it in equal portions to the given core symbols
    for (c = alph->comprise[i]; *c != 0; c++) {
      core_count = get_array_item((*c) - 1, freqs);
      // Add the ambiguous counts.
      core_count = LOG_SUM(core_count, ambig_count);
      set_array_item((*c) - 1, core_count, freqs);
    }
    // set the ambiguous count to zero
    set_array_item(i - 1, LOG_ZERO, freqs);
  }
}

/*
 * complement_swap_freqs
 * swap complementary positions between arrays and recalulate ambiguous values
 * when the arrays have extra allocated space.
 * You may pass the same array twice to just complement it.
 */
void complement_swap_freqs(ALPH_T *alph, ARRAY_T* a1, ARRAY_T* a2) {
  int i, c;
  PROB_T temp;
  assert(alph_has_complement(alph));
  assert(get_array_length(a1) >= alph->ncore);
  assert(get_array_length(a2) >= alph->ncore);
  for (i = 1; i <= alph->ncore; i++) {
    c = alph->complement[i];
    if (alph->complement[c] == i && i < c) { // complement pair
      // swap with complement
      temp = get_array_item(i - 1, a1);
      set_array_item(i - 1, get_array_item(c - 1, a2), a1);
      set_array_item(c - 1, temp, a2);
      if (a1 != a2) {
        // only do this if they are different because
        // otherwise we'll just swap them back!
        temp = get_array_item(c - 1, a1);
        set_array_item(c - 1, get_array_item(i - 1, a2), a1);
        set_array_item(i - 1, temp, a2);
      }
    }
  }
  if (get_array_length(a1) >= alph->nfull) calc_ambigs(alph, false, a1);
  if (a1 != a2 && get_array_length(a2) >= alph->nfull) calc_ambigs(alph, false, a2);
}

/*
 * translate_seq
 * Converts a sequence in place to a more restricted form:
 * ALPH_NO_ALIASES - convert any alias symbols into the primary representation.
 * ALPH_NO_AMBIGS - convert any ambiguous symbols into the main wildcard symbol.
 * ALPH_NO_UNKNOWN - convert any unknown symbols into the main wildcard symbol.
 *
 * TODO - this could possibly be made faster by the use of the other encode arrays.
 */
void translate_seq(ALPH_T *alph, char *sequence, int flags) {
  char *s;
  uint8_t idx;
  char wildcard;
  bool no_aliases, no_ambigs, no_unknown;
  no_aliases = ((flags & ALPH_NO_ALIASES) != 0);
  no_ambigs = ((flags & ALPH_NO_AMBIGS) != 0);
  no_unknown = ((flags & ALPH_NO_UNKNOWN) != 0);
  wildcard = alph_wildcard(alph);
  for (s = sequence; *s != '\0'; s++) {
    idx = alph->encode[(uint8_t)(*s)];
    if (idx == 0) {
      if (no_unknown) *s = wildcard;
    } else if (no_ambigs && idx > alph->ncore) {
      *s = wildcard;
    } else if (no_aliases) {
      // simply overwrite with the primary symbol
      *s = alph->symbols[idx];
    }
  }
}

/*
 * invcomp_seq
 * Converts a sequence in place to its reverse complement.
 */
void invcomp_seq(ALPH_T *alph, char *sequence, long length) {
  char *sl, *sr, tmp;
  sl = sequence; 
  sr = sequence+(length - 1);
  for(; sl <= sr; sl++, sr--) {
    tmp = comp_sym(alph, *sl);
    *sl = comp_sym(alph, *sr);
    *sr = tmp;
  }
}

//*************************************************************************
// Get log of the probability of each position in a string given the
// conditional probabilities (esentially a background model that is normalized
// for each prefix). The background must also contain an entry for the wildcard
// at minimum though it can contain entries for all the alphabet.
// 
// Returns the cumulative background as an array:
//    logcumback_i = 0, i=0
//        = log Pr(s_{0,i-1} | H_0), otherwise.
// and the total (log) cumulative background probability.
// 
// The probability of any length-w substring starting at position i in the
// string (in the context of the string) can then be computed as:
//    last_p = i+w-1;
//    log_p = logcumback[last_p+1] - logcumback[i];
// 
// The background model, H_0, is a Markov model defined by the
// order, n, and the conditional probabilities, a_cp, where
//    a_cp[s2i(wa)] = Pr(a | w).
//*************************************************************************
double calculate_log_cumulative_background(ALPH_T *alph, bool wildcard_only, int order, ARRAY_T *a_cp, const char *seq, double *logcumback) {
  int i, j, index, asize, *history;
  double log_pwa;
  uint8_t *encode;
  asize = (wildcard_only ? alph_size_wild(alph) : alph_size_full(alph));
  encode = (wildcard_only ? alph->encode2core : alph->encode);
  history = mm_malloc(sizeof(int) * (order + 1));
  for (i = 0; i <= order; i++) history[i] = 0;
  logcumback[0] = 0.0;
  for (i = 0; seq[i] != '\0'; i++) {
    index = encode[(uint8_t)(seq[i])];
    log_pwa = 0;
    if (index) {
      for (j = order; j; j--) {
        if (!history[j - 1]) continue;
        history[j] = (asize * history[j-1]) + index;
      }
      history[0] = index;
      for (j = order; j >= 0; j--) {
        if (history[j]) {
          log_pwa = log(get_array_item(history[j] - 1, a_cp));
          break;
        }
      }
    } else {
      for (j = 0; j <= order; j++) history[j] = 0;
    }
    logcumback[i+1] = logcumback[i] + log_pwa;
  }
  free(history);
  return logcumback[i]; // total cumulative prob. 
}

/*
 * xdna2protein
 * Get a translation of DNA to protein.
 */
XLATE_T* xlate_dna2protein() {
  ALPH_T *dna_alph, *protein_alph;
  XLATE_READER_T *reader;
  XLATE_T *translator;
  dna_alph = alph_dna();
  protein_alph = alph_protein();
  reader = xlate_reader_create(dna_alph, protein_alph);
  // Phenylalanine
  xlate_reader_add(reader, "TTT", "F");
  xlate_reader_add(reader, "TTC", "F");
  // Leucine
  xlate_reader_add(reader, "TTA", "L");
  xlate_reader_add(reader, "TTG", "L");
  xlate_reader_add(reader, "CTT", "L");
  xlate_reader_add(reader, "CTC", "L");
  xlate_reader_add(reader, "CTA", "L");
  xlate_reader_add(reader, "CTG", "L");
  // Isoleucine
  xlate_reader_add(reader, "ATT", "I");
  xlate_reader_add(reader, "ATC", "I");
  xlate_reader_add(reader, "ATA", "I");
  // Methionine
  xlate_reader_add(reader, "ATG", "M");
  // Valine
  xlate_reader_add(reader, "GTT", "V");
  xlate_reader_add(reader, "GTC", "V");
  xlate_reader_add(reader, "GTA", "V");
  xlate_reader_add(reader, "GTG", "V");
  // Serine
  xlate_reader_add(reader, "TCT", "S");
  xlate_reader_add(reader, "TCC", "S");
  xlate_reader_add(reader, "TCA", "S");
  xlate_reader_add(reader, "TCG", "S");
  // Proline
  xlate_reader_add(reader, "CCT", "P");
  xlate_reader_add(reader, "CCC", "P");
  xlate_reader_add(reader, "CCA", "P");
  xlate_reader_add(reader, "CCG", "P");
  // Threonine
  xlate_reader_add(reader, "ACT", "T");
  xlate_reader_add(reader, "ACC", "T");
  xlate_reader_add(reader, "ACA", "T");
  xlate_reader_add(reader, "ACG", "T");
  // Alanine
  xlate_reader_add(reader, "GCT", "A");
  xlate_reader_add(reader, "GCC", "A");
  xlate_reader_add(reader, "GCA", "A");
  xlate_reader_add(reader, "GCG", "A");
  // Tyrosine
  xlate_reader_add(reader, "TAT", "Y");
  xlate_reader_add(reader, "TAC", "Y");
  // Stop codon (Ochre) "TAA"
  //xlate_reader_add(reader, "TAA", "*");
  // Stop codon (Amber) "TAG"
  //xlate_reader_add(reader, "TAG", "*");
  // Histidine
  xlate_reader_add(reader, "CAT", "H");
  xlate_reader_add(reader, "CAC", "H");
  // Glutamine
  xlate_reader_add(reader, "CAA", "Q");
  xlate_reader_add(reader, "CAG", "Q");
  // Asparagine
  xlate_reader_add(reader, "AAT", "N");
  xlate_reader_add(reader, "AAC", "N");
  // Lysine
  xlate_reader_add(reader, "AAA", "K");
  xlate_reader_add(reader, "AAG", "K");
  // Aspartic acid
  xlate_reader_add(reader, "GAT", "D");
  xlate_reader_add(reader, "GAC", "D");
  // Glutamic acid
  xlate_reader_add(reader, "GAA", "E");
  xlate_reader_add(reader, "GAG", "E");
  // Cysteine
  xlate_reader_add(reader, "TGT", "C");
  xlate_reader_add(reader, "TGC", "C");
  // Stop codon (Opal) "TGA"
  //xlate_reader_add(reader, "TGA", "*");
  // Tryptophan
  xlate_reader_add(reader, "TGG", "W");
  // Arginine
  xlate_reader_add(reader, "CGT", "R");
  xlate_reader_add(reader, "CGC", "R");
  xlate_reader_add(reader, "CGA", "R");
  xlate_reader_add(reader, "CGG", "R");
  // Serine
  xlate_reader_add(reader, "AGT", "S");
  xlate_reader_add(reader, "AGC", "S");
  // Arginine (continued)
  xlate_reader_add(reader, "AGA", "R");
  xlate_reader_add(reader, "AGG", "R");
  // Glycine
  xlate_reader_add(reader, "GGT", "G");
  xlate_reader_add(reader, "GGC", "G");
  xlate_reader_add(reader, "GGA", "G");
  xlate_reader_add(reader, "GGG", "G");
  // do final calcualtions
  xlate_reader_done(reader);
  // print warnings and errors
  while (xlate_reader_has_message(reader)) {
    PARMSG_T *msg;
    msg = xlate_reader_next_message(reader);
    parmsg_print(msg, stderr);
    parmsg_destroy(msg);
  }
  // get the translator
  translator = xlate_reader_translator(reader);
  if (translator == NULL) die("Failed to create DNA to protein translator!");
  xlate_reader_destroy(reader);
  alph_release(dna_alph);
  alph_release(protein_alph);
  return translator;
}

/*
 * xlate_destroy
 * Destroys the passed alphabet translator object.
 */
void xlate_destroy(XLATE_T *translator) {
  alph_release(translator->src_alph);
  alph_release(translator->dest_alph);
  free(translator->xlate);
  memset(translator, 0, sizeof(XLATE_T));
  free(translator);
}

/*
 * xlate_print
 * Print out all translations that don't equal the wildcard.
 */
void xlate_print(XLATE_T *translator, FILE *out) {
  int *symi;
  int i, j, index, wild;
  wild = alph_wild(translator->dest_alph) + 1;
  symi = mm_malloc(sizeof(int) * translator->src_nsym);
  // init loop
  for (i = 0; i < translator->src_nsym; i++) symi[i] = 0;
  i = translator->src_nsym - 1;
  index = 0;
  // test loop
  while (i >= 0) {
    // print if the resulting symbol is not the wildcard
    if (translator->xlate[index] != 0 && translator->xlate[index] != wild) {
      for (j = 0; j < translator->src_nsym; j++) {
        fprintf(out, "%c", alph_char(translator->src_alph, symi[j] - 1));
      }
      fprintf(out, " = %c\n", alph_char(translator->dest_alph, translator->xlate[index] - 1));
    }
    // increment loop
    for (i = translator->src_nsym - 1; i >= 0; i--) {
      symi[i]++;
      if (symi[i] <= alph_size_full(translator->src_alph)) {
        break;
      } else {
        symi[i] = 0;
      }
    }
    index++;
  }
  free(symi);
}

/*
 * alph_dhash
 * Hash sequence to index in two-letter hashed alphabet.
 */
static inline int *alph_dhash(ALPH_T *alph, bool full_alph, char *sequence, long length) {
  long i, alen;
  int *hash_seq, *h;
  char *s;

  alen = (full_alph ? alph_size_full(alph) : alph_size_wild(alph));
  // Hash sequence from letters to positions in alphabet.
  // Pad the hash sequence with wildcards so even the last real position can
  // be double hashed.
  hash_seq = mm_malloc(sizeof(int) * (length + 1));
  if (full_alph) {
    for (i = 0, s = sequence, h = hash_seq; i < length; i++, s++, h++)
      *h = alph_encode(alph, *s);
  } else {
    for (i = 0, s = sequence, h = hash_seq; i < length; i++, s++, h++)
      *h = alph_encodec(alph, *s);
  }
  *h = alen; // pad with special "blank" character

  //Hash sequence to "double letter" logodds alphabet.
  for (i = 0, h = hash_seq; i < length; i++, h++) {
    // alen + 1 is used to give space for the "blank" character
    *h = ((alen + 1) * (*h)) + (*(h+1));
  }

  return hash_seq;
}

/*
 * xlate_dhash
 * Hash sequence to index in two-letter hashed alphabet.
 * If translating, hash two translated letters (eg codons) to two-letter index.
 */
static inline int *xlate_dhash(XLATE_T *translator, bool full_alph, char *sequence, long length) {
  long i, len, alen;
  int *hash_seq, *h, inc;
  char *s;
  ALPH_T *alph;

  alph = xlate_dest_alph(translator);
  alen = (full_alph ? alph_size_full(alph) : alph_size_wild(alph));
  // Hash sequence from letters to positions in alphabet.
  // Pad the hash sequence with wildcards so even the last real position can
  // be double hashed.
  inc = xlate_src_nsyms(translator); // distance to next letter
  len = length - (inc - 1); // last full letter or codon
  hash_seq = mm_malloc(sizeof(int) * (len + inc));
  for (i = 0, s = sequence, h = hash_seq; i < len; i++, s++, h++) {
    *h = xlate_index(translator, false, s);
    if (*h >= alen) *h = alph_wild(alph);
  }
  for ( ; i < len + inc; i++, h++) *h = alen; // pad with special "blank" character

  //Hash sequence to "double letter" logodds alphabet.
  for (i = 0, h = hash_seq; i < len; i++, h++) {
    // alen + 1 is used to give space for the "blank" character
    *h = ((alen + 1) * (*h)) + (*(h+inc));
  }

  return hash_seq;
}

/*
 * dhash_seq
 * Hash sequence to index in two-letter hashed alphabet.
 * If translating, hash two translated letters (eg codons) to two-letter index.
 *
 */
int* dhash_seq(ALPH_T *alph, XLATE_T *trans, bool full_alph, char *sequence, long length) {
  if (trans) return xlate_dhash(trans, full_alph, sequence, length);
  return alph_dhash(alph, full_alph, sequence, length);
}
