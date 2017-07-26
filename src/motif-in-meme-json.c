

#include "json-checker.h"
#include "alph-in.h"
#include "linked-list.h"
#include "motif-in-flags.h"
#include "motif-in-meme-json.h"
#include "motif-spec.h"
#include "string-builder.h"

typedef struct {
  char symbol;
  char complement;
  char *comprise;
  char *name;
  int colour;
  char *aliases;
} SYM_T;

typedef struct meme_data DATA_T;
struct meme_data {
  int options;
  int options_found;
  int options_returned;
  int version[3];
  int strands;
  int alph_ncore;
  int alph_flags;
  char *alph_name;
  SYM_T alph_sym;
  char *release;
  ALPH_T* alphabet;
  ARRAY_T *background;
  ARRAYLST_T *scanned_sites;
  int format_match;
  LINKLST_T *warnings;
  LINKLST_T *errors;
  LINKLST_T *motifs;
};

typedef struct mhtml2 MHTML2_T;
struct mhtml2 {
  JSONCHK_T *reader;
  DATA_T *data;
};

/*****************************************************************************
 * An error occurred parsing the MEME HTML+JSON file
 ****************************************************************************/
void mhtml2_verror(void *user_data, char *format, va_list ap) {
  DATA_T *data;
  int len;
  char dummy[1], *msg;
  va_list ap_copy;
  data = (DATA_T*)user_data;
  va_copy(ap_copy, ap);
  len = vsnprintf(dummy, 1, format, ap_copy);
  va_end(ap_copy);
  msg = mm_malloc(sizeof(char) * (len + 1));
  vsnprintf(msg, len + 1, format, ap);
  linklst_add(msg, data->errors);
}

static void error(void *user_data, char *format, ...) {
  va_list ap;
  va_start(ap, format);
  mhtml2_verror(user_data, format, ap);
  va_end(ap);
}

static void warning(void *user_data, char *format, ...) {
  DATA_T *data;
  va_list  argp;
  int len;
  char dummy[1], *msg;
  data = (DATA_T*)user_data;
  va_start(argp, format);
  len = vsnprintf(dummy, 1, format, argp);
  va_end(argp);
  msg = mm_malloc(sizeof(char) * (len + 1));
  va_start(argp, format);
  vsnprintf(msg, len + 1, format, argp);
  va_end(argp);
  linklst_add(msg, data->warnings);
}

/*****************************************************************************
 * callbacks
 ****************************************************************************/
bool mhtml2_program(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  if (strcmp(value, "MEME") == 0 || strcmp(value, "dreme") == 0) {
    // by this point we're pretty sure we've got the right parser
    ((DATA_T*)user_data)->format_match = 3;
    return true;
  } else {
    error(user_data, "Property \"program\" is not set to MEME or dreme.");
    return false;
  }
}

bool mhtml2_version(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T* data;
  int i;
  char *start, *end;
  data = (DATA_T*)user_data;
  for (i = 0, start = (char*)value; i < 3; i++) {
    data->version[i] = strtol(start, &end, 10);
    if (start == end || !(*end == '\0' || *end == '.')) {
      error(data, "The version string is incorrectly formatted.");
      return false;
    }
    if (*end == '\0') break;
    start = end+1;
  }
  for (i++; i < 3; i++) data->version[i] = 0;
  return true;
}

bool mhtml2_release(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  ((DATA_T*)user_data)->release = strdup(value);
  return true;
}

bool mhtml2_alph_name(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  ((DATA_T*)user_data)->alph_name = strdup(value);
  return true;
}

bool mhtml2_alph_extends(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T* data;
  data = (DATA_T*)user_data;
  if (strcmp("rna", value) == 0) {
    data->alph_flags = ALPH_FLAG_EXTENDS_RNA;
  } else if (strcmp("dna", value) == 0) {
    data->alph_flags = ALPH_FLAG_EXTENDS_DNA;
  } else if (strcmp("protein", value) == 0) {
    data->alph_flags = ALPH_FLAG_EXTENDS_PROTEIN;
  } else {
    error(data, "The alphabet extends string is not a standard alphabet.");
    return false;
  }
  return true;
}

bool mhtml2_alph_ncore(void *user_data, void *object, const char *property, long double value) {
  ((DATA_T*)user_data)->alph_ncore = (int)value;
  return true;
}

bool mhtml2_alph_symbols(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T* data;
  data = (DATA_T*)owner;
  data->alphabet = alph_type(value, value_len);
  return data->alphabet != NULL;
}
bool mhtml2_alph_symbols2(void *user_data, void *object, const char *property, void *value) {
  DATA_T* data;
  data = (DATA_T*)user_data;
  data->alphabet = (ALPH_T*)value;
  return data->alphabet != NULL;
}

/*
 * Create the alphabet reader
 */
void* mhtml2_alph_start(void* user_data, void* owner, const int *index) {
  DATA_T *data;
  ALPH_READER_T *reader;
  data = (DATA_T*)user_data;
  if (data->alph_ncore == 0) {
    error(data, "The alphabet property ncore was not specified before the symbols.");
  }
  reader = alph_reader_create();
  alph_reader_header(reader, 1, data->alph_name, data->alph_flags);
  return reader;
}
/*
 * Finish with the alphabet reader and return the ALPH_T
 */
void* mhtml2_alph_end(void* user_data, void* owner, void* self) {
  ALPH_T *alph;
  ALPH_READER_T *reader;
  DATA_T* data;
  PARMSG_T *msg;
  data = (DATA_T*)owner;
  reader = (ALPH_READER_T*)self;
  alph_reader_done(reader);
  while (alph_reader_has_message(reader)) {
    msg = alph_reader_next_message(reader);
    if (msg->severity == SEVERITY_ERROR) {
      error(data, "Alphabet error: %s.", msg->message);
    } else {
      warning(data, "Alphabet warning: %s.", msg->message);
    }
    parmsg_destroy(msg);
  }
  alph = alph_reader_alphabet(reader);
  if (alph == NULL) {
    error(data, "Failed to determine the alphabet.");
  }
  alph_reader_destroy(reader);
  return alph;
}

/*
 * Free the alphabet reader.
 */
void mhtml2_alph_abort(void* self) {
  alph_reader_destroy((ALPH_READER_T*)self);
}

/*
 * Add a symbol entry to the alphabet
 */
bool mhtml2_alph_entry(void *user_data, void *list, const int *index, void *value) {
  DATA_T *data;
  ALPH_READER_T *reader;
  SYM_T *sym;
  char *alias_p, sym_comprise[2];
  data = (DATA_T*)user_data;
  reader = (ALPH_READER_T*)list;
  sym = (SYM_T*)value;
  if (index[0] < data->alph_ncore) {
    // core symbol
    alph_reader_core(reader, sym->symbol, sym->aliases, sym->name, sym->colour, sym->complement);
  } else {
    // ambiguous symbol
    if (sym->comprise != NULL) {
      alph_reader_ambig(reader, sym->symbol, sym->aliases, sym->name, sym->colour, sym->comprise);
    } else {
      error(data, "Alphabet ambiguous symbol %c requires attribute \"equals\".", sym->symbol);
      return false;
    }
  }
  return true;
}

void* mhtml2_alph_sym_reset(void* user_data, void* owner, const int *index) {
  DATA_T *data;
  SYM_T *sym;
  data = (DATA_T*)user_data;
  sym = &(data->alph_sym);
  sym->symbol = '\0';
  sym->complement = '\0';
  if (sym->comprise) free(sym->comprise);
  sym->comprise = NULL;
  if (sym->name) free(sym->name);
  sym->name = NULL;
  sym->colour = -1;
  if (sym->aliases) free(sym->aliases);
  sym->aliases = NULL;
  return sym;
}

bool mhtml2_alph_sym_symbol(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T *data;
  SYM_T *sym;
  data = (DATA_T*)user_data;
  sym = (SYM_T*)owner;
  if (value_len == 1) {
    sym->symbol = value[0];
    return true;
  } else {
    error(data, "The alphabet symbol should be exactly one character.");
    return false;
  }
}
bool mhtml2_alph_sym_name(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  ((SYM_T*)owner)->name = strdup(value);
  return true;
}
bool mhtml2_alph_sym_colour(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T *data;
  SYM_T *sym;
  char *end;
  data = (DATA_T*)user_data;
  sym = (SYM_T*)owner;
  sym->colour = strtol(value, &end, 16);
  if (value_len == 6 && value != end && *end == '\0') {
    return true;
  } else {
    error(data, "The alphabet symbol colour must be 6 hexadecimal digits.");
    return false;
  }
}
bool mhtml2_alph_sym_complement(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T *data;
  SYM_T *sym;
  data = (DATA_T*)user_data;
  sym = (SYM_T*)owner;
  if (value_len == 1) {
    sym->complement = value[0];
    return true;
  } else {
    error(data, "The alphabet symbol complement should be exactly one character.");
    return false;
  }
}
bool mhtml2_alph_sym_equals(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  ((SYM_T*)owner)->comprise = strdup(value);
  return true;
}
bool mhtml2_alph_sym_aliases(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  ((SYM_T*)owner)->aliases = strdup(value);
  return true;
}

bool mhtml2_alph_strands(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  DATA_T* data;
  data = (DATA_T*)owner;
  if (strcmp("both", value) == 0) {
    data->strands = 2;
  } else if (strcmp("forward", value) == 0) {
    data->strands = 1;
  } else if (strcmp("none", value) == 0) {
    data->strands = 0;
  } else {
    error(data, "The strands value should be either \"both\", "
        "\"forward\" or \"none\".");
    return false;
  }
  return true;
}

bool mhtml2_alph_revcomp(void *user_data, void *owner, const char *property, bool value) {
  DATA_T* data;
  data = (DATA_T*)owner;
  data->strands = (value ? 2 : 1);
  return true;
}

bool mhtml2_bgfreqs(void *user_data, void *owner, const char *property, void *value) {
  DATA_T* data;
  data = (DATA_T*)owner;
  data->background = (ARRAY_T*)value;
  return true;
}

void* mhtml2_freqs_create(void *user_data, void *owner, const int *index) {
  DATA_T* data;
  data = (DATA_T*)user_data;
  if (data->alphabet != NULL) {
    return allocate_array(alph_size_core(data->alphabet));
  }
  error(data, "The alphabet should be listed in the file before any frequencies");
  return NULL;
}

void* mhtml2_freqs_finalize(void *user_data, void *owner, void *value) {
  DATA_T* data;
  ARRAY_T* freqs;
  long double sum, delta;
  data = (DATA_T*)user_data;
  freqs = (ARRAY_T*)value;
  sum = array_total(freqs);
  delta = sum - 1.0;
  if (delta < 0) delta = -delta;
  if (delta > 0.01) {
    error(data, "The frequencies do not sum to 1.");
    return NULL;
  }
  return freqs;
}

void mhtml2_freqs_abort(void *value) {
  free_array((ARRAY_T*)value);
}

bool mhtml2_freqs_set(void *user_data, void *list, const int *index, long double value) {
  DATA_T *data;
  ARRAY_T *freqs;
  freqs = (ARRAY_T*)list;
  if (index[0] >= get_array_length(freqs)) {
    data = (DATA_T*)user_data;
    error(data, "Too many frequencies for the alphabet.");
    return false;
  }
  set_array_item(index[0], value, freqs);
  return true;
}


bool mhtml2_motif_set(void *user_data, void *owner, const int *index, void *value) {
  DATA_T* data;
  data = (DATA_T*)owner;
  linklst_add(value, data->motifs);
  return true;
}

void* mhtml2_motif_create(void *user_data, void *owner, const int *index) {
  DATA_T *data;
  MOTIF_T *motif;
  // access the data object to get file scoped variables
  data = (DATA_T*)user_data;
  // setup motif
  motif = mm_malloc(sizeof(MOTIF_T));
  memset(motif, 0, sizeof(MOTIF_T));
  motif->alph = alph_hold(data->alphabet);
  motif->flags = (alph_has_complement(data->alphabet) && data->strands == 2 ? MOTIF_BOTH_STRANDS : 0);
  set_motif_strand('+', motif);
  motif->url = strdup("");
  return motif;
}

void* mhtml2_motif_finalize(void *user_data, void *owner, void *value) {
  //TODO checks
  return value;
}

void mhtml2_motif_abort(void *value) {
  destroy_motif((MOTIF_T*)value);
}

bool mhtml2_motif_id(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  set_motif_id(value, value_len, (MOTIF_T*)owner);
  return true;
}

bool mhtml2_motif_alt(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  set_motif_id2(value, value_len, (MOTIF_T*)owner);
  return true;
}

bool mhtml2_motif_evalue(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  MOTIF_T *motif;
  motif = (MOTIF_T*)owner;
  motif->log_evalue = log10_evalue_from_string(value);
  motif->evalue = pow(10.0, motif->log_evalue);
  return true;
}

bool mhtml2_motif_len(void *user_data, void *owner, const char *property, long double value) {
  DATA_T *data;
  MOTIF_T *motif;
  data = (DATA_T*)user_data;
  motif = (MOTIF_T*)owner;
  // check for negative values
  if (value <= 0) {
    error(data, "Expected positive and non-zero motif length.");
    return false;
  }
  motif->length = (int)value;
  // check for non-integers
  if (value != motif->length) {
    error(data, "Expected integer for motif length.");
    return false;
  }
  return true;
}

bool mhtml2_motif_nsites(void *user_data, void *owner, const char *property, long double value) {
  if (value < 0) {
    error(((DATA_T*)user_data), "Expected positive value for the number of motif sites.");
    return false;
  }
  ((MOTIF_T*)owner)->num_sites = (double)value;
  return true;
}

bool mhtml2_motif_ic(void *user_data, void *owner, const char *property, long double value) {
  ((MOTIF_T*)owner)->complexity = (double)value;
  return true;
}

void* mhtml2_matrix_setup(void *user_data, void *owner, const int *index) {
  DATA_T *data;
  MOTIF_T *motif;
  data = (DATA_T*)user_data;
  motif = (MOTIF_T*)owner;
  if (data->alphabet == NULL) {
    error(data, "The alphabet should be listed in the file before any motifs.");
    return false;
  }
  if (motif->length == 0) {
    error(data, "The motif length should be listed in the file before the "
        "score or probability matricies.");
    return false;
  }
  return allocate_matrix(motif->length,
      alph_size_core(data->alphabet));
}

void* mhtml2_matrix_finalize(void *user_data, void *owner, void *value) {
  MATRIX_T *matrix;
  matrix = (MATRIX_T*)value;
  //TODO validate
  return matrix;
}

void mhtml2_matrix_abort(void *value) {
  free_matrix((MATRIX_T*)value);
}

bool mhtml2_matrix_set(void *user_data, void *owner, const int *index, long double value) {
  DATA_T *data;
  MATRIX_T *matrix;
  data = (DATA_T*)user_data;
  matrix = (MATRIX_T*)owner;
  if (index[0] > get_num_rows(matrix)) {
    error(data, "The matrix has too many rows.");//TODO better message
    return false;
  }
  if (index[1] > get_num_cols(matrix)) {
    error(data, "The matrix has too many columns.");//TODO better message
    return false;
  }
  set_matrix_cell(index[0], index[1], (double)value, matrix);
  return true;
}

bool mhtml2_psm_set(void *user_data, void *owner, const char *property, void *value) {
  ((MOTIF_T*)owner)->scores = (MATRIX_T*)value;
  return true;
}

bool mhtml2_pwm_set(void *user_data, void *owner, const char *property, void *value) {
  ((MOTIF_T*)owner)->freqs = (MATRIX_T*)value;
  return true;
}

void* mhtml2_sseqs_create(void *user_data, void *owner, const int *index) {
  DATA_T *data;
  ARRAYLST_T *sseqs;
  sseqs = arraylst_create();
  data = (DATA_T*)user_data;
  data->scanned_sites = sseqs;
  return sseqs;
}

void* mhtml2_sseq_create(void *user_data, void *owner, const int *index) {
  ARRAYLST_T *sseqs;
  SCANNED_SEQ_T *sseq;
  sseqs = (ARRAYLST_T*)owner;
  sseq = mm_malloc(sizeof(SCANNED_SEQ_T));
  memset(sseq, 0, sizeof(SCANNED_SEQ_T));
  arraylst_add(sseq, sseqs);
  return sseq;
}

bool mhtml2_sseq_name(void *user_data, void *owner, const char *property, const char *value, size_t value_len) {
  SCANNED_SEQ_T *sseq;
  sseq = (SCANNED_SEQ_T*)owner;
  sseq->seq_id = strdup(value);
  return true;
}

bool mhtml2_sseq_length(void *user_data, void *owner, const char *property, long double value) {
  SCANNED_SEQ_T *sseq;
  sseq = (SCANNED_SEQ_T*)owner;
  sseq->length = (long)value;
  return true;
}

void* mhtml2_sseq_select(void *user_data, void *owner, const int *index) {
  DATA_T *data;
  ARRAYLST_T *sseqs;
  data = (DATA_T*)user_data;
  sseqs = data->scanned_sites;
  if (index[0] >= arraylst_size(sseqs)) {
    error(data, "More scanned sequences than were listed in the dataset section!");
    return NULL;
  }
  // lookup the existing entry we already created
  return arraylst_get(index[0], sseqs);
}

bool mhtml2_sseq_pv(void *user_data, void *owner, const char *property, long double value) {
  SCANNED_SEQ_T *sseq;
  sseq = (SCANNED_SEQ_T*)owner;
  sseq->pvalue = (double)value;
  return true;
}

void* mhtml2_ssite_add(void *user_data, void *owner, const int *index) {
  SCANNED_SITE_T *ssite;
  SCANNED_SEQ_T *sseq;
  sseq = (SCANNED_SEQ_T*)owner;
  // as index is generated by json-checker it's only going to increment
  assert(index[0] == sseq->site_count);
  sseq->site_count++;
  sseq->sites = mm_realloc(sseq->sites, sizeof(SCANNED_SITE_T) * sseq->site_count);
  ssite = sseq->sites+(index[0]);
  memset(ssite, 0, sizeof(SCANNED_SITE_T));
  return ssite;
}

bool mhtml2_ssite_motif(void *user_data, void *owner, const char *property, long double value) {
  SCANNED_SITE_T *ssite;
  ssite = (SCANNED_SITE_T*)owner;
  ssite->m_num = (int)value;
  return true;
}

bool mhtml2_ssite_pos(void *user_data, void *owner, const char *property, long double value) {
  SCANNED_SITE_T *ssite;
  ssite = (SCANNED_SITE_T*)owner;
  ssite->position = (long)value;
  return true;
}

bool mhtml2_ssite_rc(void *user_data, void *owner, const char *property, bool value) {
  SCANNED_SITE_T *ssite;
  ssite = (SCANNED_SITE_T*)owner;
  ssite->strand = (value ? '-' : '+');
  return true;
}

bool mhtml2_ssite_pv(void *user_data, void *owner, const char *property, long double value) {
  SCANNED_SITE_T *ssite;
  ssite = (SCANNED_SITE_T*)owner;
  ssite->pvalue = (double)value;
  return true;
}

bool mhtml2_scan_ready(void *user_data, void *owner, const char *propery, void *value) {
  DATA_T *data;
  data = (DATA_T*)user_data;
  data->options_found |= SCANNED_SITES;
  return true;
}

/*****************************************************************************
 * Define the structure of the JSON document
 ****************************************************************************/
JSON_OBJ_DEF_T *json_def(bool scan) {
  JsonSetupFn sseq_select, ssite_add, sseq_create, sseqs_create;
  JsonPropStringFn sseq_name;
  JsonPropNumberFn sseq_pv, ssite_motif, ssite_pos, ssite_pv, sseq_length;
  JsonPropBoolFn ssite_rc;
  JsonPropListFn scan_ready;
  if (scan) {
    sseqs_create = mhtml2_sseqs_create;
    sseq_create = mhtml2_sseq_create;
    sseq_name = mhtml2_sseq_name;
    sseq_length = mhtml2_sseq_length;
    sseq_select = mhtml2_sseq_select;
    sseq_pv = mhtml2_sseq_pv;
    ssite_add = mhtml2_ssite_add;
    ssite_motif = mhtml2_ssite_motif;
    ssite_pos = mhtml2_ssite_pos;
    ssite_rc = mhtml2_ssite_rc;
    ssite_pv = mhtml2_ssite_pv;
    scan_ready = mhtml2_scan_ready;
  } else {
    sseqs_create = NULL;
    sseq_create = NULL; sseq_name = NULL; sseq_length = NULL;
    sseq_select = NULL; sseq_pv = NULL; ssite_add = NULL;
    ssite_motif = NULL; ssite_pos = NULL; ssite_rc = NULL; ssite_pv = NULL;
    scan_ready = NULL;
  }

  return 
    jd_obj(NULL, NULL, NULL, true, 10,
      jd_pstr("program", true, mhtml2_program),
      jd_pstr("version", true, mhtml2_version),
      jd_pstr("release", true, mhtml2_release),
      jd_plst("cmd", true, NULL, 
        jd_lstr(1, NULL, NULL, NULL, NULL)
      ),
      jd_pobj("options", true, NULL, 
        jd_obj(NULL, NULL, NULL, true, 1, // there are lots of other options but they're not important for reading motifs
          jd_pbool("revcomp", false, mhtml2_alph_revcomp)
        )
      ),
      jd_pobj("alphabet", true, NULL,
        jd_obj(NULL, NULL, NULL, false, 6,
          jd_pstr("name", false, mhtml2_alph_name),
          jd_pstr("like", false, mhtml2_alph_extends),
          jd_pnum("ncore", false, mhtml2_alph_ncore),
          // originally the symbols were just a string with all the core symbols
          // of the alphabet but with the addition of custom alphabets it was
          // changed to a list of objects
          jd_pcust("symbols", true, JSON_STRING | JSON_LIST, mhtml2_alph_symbols,
            mhtml2_alph_symbols2, 
            jd_lobj(1, mhtml2_alph_start, mhtml2_alph_end, mhtml2_alph_abort, mhtml2_alph_entry,
              jd_obj(mhtml2_alph_sym_reset, NULL, NULL, false, 6,
                jd_pstr("symbol", true, mhtml2_alph_sym_symbol),
                jd_pstr("name", false, mhtml2_alph_sym_name),
                jd_pstr("colour", false, mhtml2_alph_sym_colour),
                jd_pstr("complement", false, mhtml2_alph_sym_complement),
                jd_pstr("equals", false, mhtml2_alph_sym_equals),
                jd_pstr("aliases", false, mhtml2_alph_sym_aliases)
              )
            )  
          ),
          jd_pstr("strands", false, mhtml2_alph_strands), // removed, now served by combination of alphabet and revcomp
          jd_plst("freqs", false, mhtml2_bgfreqs, // moved outsize alphabet and turned into object
            jd_lnum(1, mhtml2_freqs_create, mhtml2_freqs_finalize, mhtml2_freqs_abort, mhtml2_freqs_set)
          )
        )
      ),
      jd_pobj("background", false, NULL,
        jd_obj(NULL, NULL, NULL, false, 2,
          jd_pstr("source", false, NULL),
          jd_plst("freqs", true, mhtml2_bgfreqs,
            jd_lnum(1, mhtml2_freqs_create, mhtml2_freqs_finalize, mhtml2_freqs_abort, mhtml2_freqs_set)
          )
        )
      ),
      jd_pobj("sequence_db", true, NULL,
        jd_obj(NULL, NULL, NULL, true, 1,
          jd_plst("sequences", false, NULL,
            jd_lobj(1, sseqs_create, NULL, NULL, NULL,
              jd_obj(sseq_create, NULL, NULL, false, 3,
                jd_pstr("name", true, sseq_name),
                jd_pnum("length", true, sseq_length),
                jd_pnum("weight", true, NULL)
              )
            )
          )
        )
      ),  
      jd_plst("motifs", true, NULL,
        jd_lobj(1, NULL, NULL, NULL, mhtml2_motif_set,
          jd_obj(mhtml2_motif_create, mhtml2_motif_finalize, mhtml2_motif_abort, true, 9,
            jd_pnum("db", true, NULL),
            jd_pstr("id", true, mhtml2_motif_id),
            jd_pstr("alt", true, mhtml2_motif_alt),
            jd_pnum("len", true, mhtml2_motif_len),
            jd_pnum("nsites", true, mhtml2_motif_nsites),
            jd_pstr("evalue", true, mhtml2_motif_evalue),
            jd_pnum("ic", false, mhtml2_motif_ic),
            jd_plst("psm", false, mhtml2_psm_set,
              jd_lnum(2, mhtml2_matrix_setup, mhtml2_matrix_finalize, mhtml2_matrix_abort, mhtml2_matrix_set)
            ),
            jd_plst("pwm", true, mhtml2_pwm_set,
              jd_lnum(2, mhtml2_matrix_setup, mhtml2_matrix_finalize, mhtml2_matrix_abort, mhtml2_matrix_set)
            )
          )
        )
      ),
      jd_plst("scan", false, scan_ready,
        jd_lobj(1, NULL, NULL, NULL, NULL,
          jd_obj(sseq_select, NULL, NULL, true, 2,
            jd_pnum("pvalue", true, sseq_pv),
            jd_plst("sites", true, NULL,
              jd_lobj(1, NULL, NULL, NULL, NULL,
                jd_obj(ssite_add, NULL, NULL, true, 4,
                  jd_pnum("motif", true, ssite_motif),
                  jd_pnum("pos", true, ssite_pos),
                  jd_pbool("rc", true, ssite_rc),
                  jd_pnum("pvalue", true, ssite_pv)
                )
              )
            )
          )
        )
      )
    );
}

/*****************************************************************************
 * Create the parser data structure
 ****************************************************************************/
static DATA_T* create_parser_data(int options, const char *optional_filename) {
  DATA_T *data;
  data = mm_malloc(sizeof(DATA_T));
  memset(data, 0, sizeof(DATA_T));
  data->options = options;
  data->format_match = file_name_match("meme", "html", optional_filename);
  data->warnings = linklst_create();
  data->errors = linklst_create();
  data->motifs = linklst_create();
  return data;
}

/*****************************************************************************
 * Destroy the parser data structure
 ****************************************************************************/
static void destroy_parser_data(DATA_T *data) {
  SYM_T *sym;
  if (data->scanned_sites) {
    if (!(data->options_returned & SCANNED_SITES)) {
      arraylst_destroy(sseq_destroy, data->scanned_sites);
    }
  }
  if (data->background) {
    free_array(data->background);
  }
  if (data->release) {
    free(data->release);
  }
  linklst_destroy_all(data->warnings, free);
  linklst_destroy_all(data->errors, free);
  linklst_destroy_all(data->motifs, destroy_motif);
  // free alphabet related allocations
  if (data->alphabet) alph_release(data->alphabet);
  if (data->alph_name) free(data->alph_name);
  sym = &(data->alph_sym);
  if (sym->comprise) free(sym->comprise);
  if (sym->name) free(sym->name);
  if (sym->aliases) free(sym->aliases);
  memset(data, 0, sizeof(DATA_T));
  free(data);
}

/*****************************************************************************
 * Create the parser for MEME HTML with JSON data
 ****************************************************************************/
void* mhtml2_create(const char *optional_filename, int options) {
  MHTML2_T *parser;
  parser = mm_malloc(sizeof(MHTML2_T));
  memset(parser, 0, sizeof(MHTML2_T));
  parser->data = create_parser_data(options, optional_filename);
  parser->reader = jsonchk_create(mhtml2_verror, parser->data,
      json_def((options & SCANNED_SITES) != 0));
  return parser;
}

/*****************************************************************************
 * Destroy the parser
 ****************************************************************************/
void mhtml2_destroy(void *data) {
  MHTML2_T *parser;
  parser = (MHTML2_T*)data;
  jsonchk_destroy(parser->reader);
  destroy_parser_data(parser->data);
  memset(parser, 0, sizeof(MHTML2_T));
  free(parser);
}

/*****************************************************************************
 * Update the parser
 ****************************************************************************/
void mhtml2_update(void *data, const char *chunk, size_t size, short end) {
  jsonchk_update(((MHTML2_T*)data)->reader, chunk, size, end != 0);
}

/*****************************************************************************
 * Return a larger number if we think the format is correct
 ****************************************************************************/
short mhtml2_has_format_match(void *data) {
  return ((MHTML2_T*)data)->data->format_match;
}

/*****************************************************************************
 * 
 ****************************************************************************/
short mhtml2_has_warning(void *data) {
  return linklst_size(((MHTML2_T*)data)->data->warnings) > 0;
}

/*****************************************************************************
 * 
 ****************************************************************************/
char* mhtml2_next_warning(void *data) {
  return linklst_pop(((MHTML2_T*)data)->data->warnings);
}

/*****************************************************************************
 * 
 ****************************************************************************/
short mhtml2_has_error(void *data) {
  return linklst_size(((MHTML2_T*)data)->data->errors) > 0;
}

/*****************************************************************************
 * 
 ****************************************************************************/
char* mhtml2_next_error(void *data) {
  return linklst_pop(((MHTML2_T*)data)->data->errors);
}

/*****************************************************************************
 * 
 ****************************************************************************/
short mhtml2_has_motif(void *data) {
  return linklst_size(((MHTML2_T*)data)->data->motifs) > 0;
}

/*****************************************************************************
 * 
 ****************************************************************************/
MOTIF_T* mhtml2_next_motif(void *data) {
  return (MOTIF_T*)linklst_pop(((MHTML2_T*)data)->data->motifs);
}

/*****************************************************************************
 * 
 * WARNING, does not "hold" alphabet. If users wish to keep the alphabet they
 * must call alph_hold on it themselves.
 ****************************************************************************/
ALPH_T* mhtml2_get_alphabet(void *data) {
  return ((MHTML2_T*)data)->data->alphabet;
}

/*****************************************************************************
 * 
 ****************************************************************************/
int mhtml2_get_strands(void *data) {
  MHTML2_T *parser;
  parser = ((MHTML2_T*)data);
  if (parser->data->alphabet == NULL) return 0;
  return (alph_has_complement(parser->data->alphabet) ? parser->data->strands : 0);
}

/*****************************************************************************
 * 
 ****************************************************************************/
BOOLEAN_T mhtml2_get_bg(void *data, ARRAY_T **bg) {
  MHTML2_T *parser;
  parser = (MHTML2_T*)data;
  if (parser->data->background == NULL) return false;
  *bg = resize_array(*bg, get_array_length(parser->data->background));
  copy_array(parser->data->background, *bg);
  return true;
}

/*****************************************************************************
 * 
 ****************************************************************************/
void* mhtml2_motif_optional(void *data, int option) {
  MHTML2_T *parser;
  parser = (MHTML2_T*)data;
  if (!(parser->data->options & option)) {
    die("Requested value of optional component which was not requested "
        "during construction.\n");
    return NULL;
  }
  // This feature is unused
  return NULL;
}

/*****************************************************************************
 * 
 ****************************************************************************/
void* mhtml2_file_optional(void *data, int option) {
  MHTML2_T *parser;
  parser = (MHTML2_T*)data;
  if (!(parser->data->options & option)) {
    die("Requested value of optional component which was not requested "
        "during construction.\n");
    return NULL;
  }
  if (parser->data->options_found & option) {
    if (parser->data->options_returned & option) {
      die("Sorry, optional values are only returned once. "
          "This is because we can not guarantee that the "
          "previous caller did not deallocate the memory. "
          "Hence this is a feature to avoid memory bugs.\n");
      return NULL;
    }
    parser->data->options_returned |= option;
  } else {
    // Not yet found or unsupported
    return NULL;
  }
  switch (option) {
    case SCANNED_SITES:
      return parser->data->scanned_sites;
    default:
      die("Option code %d does not match any single option. "
          "This means that it must contain multiple options "
          "and only one is allowed at a time\n.", option);
  }
  return NULL; //unreachable! make compiler happy
}
