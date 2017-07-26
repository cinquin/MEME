
#include "alphabet.h"
#include "linked-list.h"
#include "parser-message.h"
#include "xlate-in.h"

struct xlate_reader {
  XLATE_T translate;
  LINKLST_T *messages;
  bool had_error;
  bool had_warning;
  bool done;
};

/*
 * add_msg
 * Record an error message for reporting to the calling program
 */
static void add_msg(XLATE_READER_T *reader, PARMSG_T *msg) {
  if (msg->severity == SEVERITY_ERROR) {
    reader->had_error = true;
  } else if (msg->severity == SEVERITY_WARNING) {
    reader->had_warning = true;
  }
  linklst_add(msg, reader->messages);
}

/*
 *
 */
static inline __attribute__((always_inline)) uint32_t s2i(ALPH_T *alph, const char *sym) {
  const char *s;
  uint32_t index;
  index = alph_index(alph, *sym) + 1;
  for (s = sym+1; *s; s++) {
    // alphabet size plus 1 for unknown symbols
    index = (index * (alph_size_full(alph) + 1)) + (alph_index(alph, *s) + 1);
  }
  return index;
}

/*
 * Check that there are the correct number of symbols and that they
 * all occur in the core alphabet.
 */
static inline __attribute__((always_inline)) bool check_symbols(
    XLATE_READER_T *reader, const char *name, ALPH_T *alph, const char *syms,
    uint8_t *len) {
  int syms_len, i;
  syms_len = strlen(syms);
  if (*len != 0 && *len != syms_len) {
    add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, -1, -1, 
          "%d %s symbols when expecting %d", syms_len, name, syms_len));
    return false;
  } else {
    *len = syms_len;
  }
  for (i = 0; i < syms_len; i++) {
    if (!alph_is_core(alph, syms[i])) {
      add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, -1, -1, 
            "%s symbol %c is not a core symbol for the %s alphabet",
            name, syms[i], alph_name(alph)));
      return false;
    }
  }
  return true;
}

/*
 *
 */
XLATE_READER_T* xlate_reader_create(ALPH_T *src, ALPH_T *dest) {
  XLATE_READER_T *reader;
  reader = mm_malloc(sizeof(XLATE_READER_T));
  reader->translate.src_alph = alph_hold(src);
  reader->translate.dest_alph = alph_hold(dest);
  reader->translate.src_nsym = 0;
  reader->translate.dest_nsym = 1; // ensure that only one destination symbol is allowed
  reader->translate.xlate = NULL;
  reader->done = false;
  reader->had_error = false;
  reader->had_warning = false;
  reader->messages = linklst_create();
  return reader;
}

/*
 *
 */
void xlate_reader_destroy(XLATE_READER_T* reader) {
  alph_release(reader->translate.src_alph);
  alph_release(reader->translate.dest_alph);
  if (reader->translate.xlate) free(reader->translate.xlate);
  linklst_destroy_all(reader->messages, parmsg_destroy);
  memset(reader, 0, sizeof(XLATE_READER_T));
  free(reader);
}

/*
 * 
 */
void xlate_reader_add(XLATE_READER_T* reader, const char *src_sym, const char *dest_sym) {
  int count;
  XLATE_T *t;
  if (reader->done) die("Reader already done!");
  t = &(reader->translate);
  if (!check_symbols(reader, "source", t->src_alph, src_sym, &(t->src_nsym))) return;
  if (!check_symbols(reader, "destination", t->dest_alph, dest_sym, &(t->dest_nsym))) return;
  if (t->xlate == NULL) {
    count = pow(alph_size_full(t->src_alph) + 1, t->src_nsym);
    t->xlate = mm_malloc(sizeof(uint32_t) * count);
    memset(t->xlate, 0, sizeof(uint32_t) * count);
  }
  t->xlate[s2i(t->src_alph, src_sym)] = alph_index(t->dest_alph, dest_sym[0]) + 1;
}

/*
 * Adds in all missing translations using the wildcard when
 * there is no better option.
 */
void xlate_reader_done(XLATE_READER_T* reader) {
  XLATE_T *t;
  int i, j, k, count, index, index2;
  int *symi, *cmpi;
  bool *targets;
  if (reader->done) return;
  t = &(reader->translate);
  assert(t->dest_nsym == 1);
  // target core symbols
  targets = mm_malloc(sizeof(bool) * alph_size_core(t->dest_alph));
  // symbol indexes
  symi = mm_malloc(sizeof(int) * t->src_nsym);
  // comprising indexes
  cmpi = mm_malloc(sizeof(int) * t->src_nsym);
  // init loop over all inputs
  for (i = 0; i < t->src_nsym; i++) symi[i] = 0;
  i = t->src_nsym - 1;
  index = 0;
  // test input
  while (i >= 0) {
    // check to see if this is set already
    if (t->xlate[index]) goto next_index;
    // Not set! 
    // First check that it's not a error state.
    for (j = 0; j < t->src_nsym; j++) {
      if (symi[j] == 0) {
        // this is an error state because one of the input symbols was unknown
        // so we translate it to the wildcard
        t->xlate[index] = alph_wild(t->dest_alph) + 1;
        goto next_index;
      }
    }
    // This might be because it contains ambiguous symbols in which case we
    // iterate over all matching core symbols and make a list of matching
    // target symbols.
    // Reset list of target core symbols
    for (j = 0; j < alph_size_core(t->dest_alph); j++) targets[j] = false;
    // init comprising symbols loop
    for (j = 0; j < t->src_nsym; j++) cmpi[j] = 0;
    j = t->src_nsym - 1;
    // comprising symbols loop
    while (j >= 0) {
      // convert into a index for the matching string of only core symbols
      index2 = 0;
      for (k = 0; k < t->src_nsym; k++) {
        index2 = (index2 * (alph_size_full(t->src_alph) + 1)) +
          (alph_comprise(t->src_alph, symi[k] - 1, cmpi[k]) + 1);
      }
      // lookup the core symbol
      if (t->xlate[index2]) {
        if (t->xlate[index2] <= alph_size_core(t->dest_alph)) {
          // mark target symbol
          targets[t->xlate[index2] - 1] = true;
        } else {
          // convert ambiguous target symbol into core target symbols
          for (k = 0; k < alph_ncomprise(t->dest_alph, t->xlate[index2] - 1); k++) {
            targets[alph_comprise(t->dest_alph, t->xlate[index2] - 1, k)] = true;
          }
        }
      } else {
        // no assignment so we must assume a wildcard!
        t->xlate[index] = alph_wild(t->dest_alph) + 1;
        goto next_index;
      }
      // increment comprising symbols loop
      for (j = t->src_nsym - 1; j >= 0; j--) {
        cmpi[j]++;
        if (cmpi[j] < alph_ncomprise(t->src_alph, symi[j] - 1)) {
          break;
        } else {
          cmpi[j] = 0;
        }
      }
    }
    // convert list of core targets into single matching target, or wildcard if none match
    count = 0;
    for (j = 0; j < alph_size_core(t->dest_alph); j++) {
      if (targets[j]) count++;
    }
    if (count == 1) {
      // only one core symbol selected so set to that symbol
      for (j = 0; j < alph_size_core(t->dest_alph); j++) {
        if (targets[j]) {
          t->xlate[index] = j + 1;
          break;
        }
      }
    } else if (count == alph_size_core(t->dest_alph)) {
      // all core characters selected so set to wildcard
      t->xlate[index] = alph_wild(t->dest_alph) + 1;
    } else {
      // test ambiguous characters
      for (j = alph_size_core(t->dest_alph); j < alph_size_full(t->dest_alph); j++) {
        if (alph_ncomprise(t->dest_alph, j) == count) {
          for (k = 0; k < count; k++) {
            if (!targets[alph_comprise(t->dest_alph, j, k)]) {
              break;
            }
          }
          if (k == count) {
            // match!
            t->xlate[index] = j + 1;
            break;
          }
        }
      }
      if (j == alph_size_full(t->dest_alph)) {
        // no match so set to wildcard
        t->xlate[index] = alph_wild(t->dest_alph) + 1;
      }
    }
    // increment input symbols
next_index:
    for (i = t->src_nsym - 1; i >= 0; i--) {
      symi[i]++;
      if (symi[i] <= alph_size_full(t->src_alph)) {
        break;
      } else {
        symi[i] = 0;
      }
    }
    index++;
  }
  free(targets);
  free(symi);
  free(cmpi);
  reader->done = true;
}

/*
 *
 */
XLATE_T* xlate_reader_translator(XLATE_READER_T* reader) {
  int count;
  XLATE_T *result, *t;
  if (reader->had_error || !reader->done) return NULL;
  t = &(reader->translate);
  result = mm_malloc(sizeof(XLATE_T));
  result->src_alph = alph_hold(t->src_alph);
  result->src_nsym = t->src_nsym;
  result->dest_alph = alph_hold(t->dest_alph);
  result->dest_nsym = t->dest_nsym;
  count = pow(alph_size_full(t->src_alph) + 1, t->src_nsym);
  result->xlate = mm_malloc(sizeof(uint32_t) * count);
  memcpy(result->xlate, t->xlate, sizeof(uint32_t) * count);
  return result;
}

/*
 * xlate_reader_had_warning
 * Return true if an input caused a recoverable warning.
 */
bool xlate_reader_had_warning(XLATE_READER_T *reader) {
  return reader->had_warning;
}

/*
 * xlate_reader_had_error
 * Return true if an input caused an unrecoverable error.
 */
bool xlate_reader_had_error(XLATE_READER_T *reader) {
  return reader->had_error;
}

/*
 * xlate_reader_has_message
 * Check if there are any pending warning or error mesages.
 */
bool xlate_reader_has_message(XLATE_READER_T *reader) {
  return (linklst_size(reader->messages) > 0);
}

/*
 * xlate_reader_next_message
 * Return the next pending warning or error message.
 * The caller is responsible for freeing memory.
 */
PARMSG_T* xlate_reader_next_message(XLATE_READER_T *reader) {
  return (PARMSG_T*)linklst_pop(reader->messages);
}

#ifdef MAIN
int main(int argc, char **argv) {
  XLATE_T *dna2prot;
  dna2prot = xlate_dna2protein();
  xlate_print(dna2prot, stdout);
  xlate_destroy(dna2prot);
  return 0;
}
#endif
