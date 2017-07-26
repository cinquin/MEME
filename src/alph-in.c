
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <regex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>

#include "alph-in.h"
#include "array-list.h"
#include "linked-list.h"
#include "parser-message.h"
#include "red-black-tree.h"
#include "regex-utils.h"
#include "string-builder.h"

// Stores the compiled regular expression used for parsing.
struct patterns {
  regex_t header;
  regex_t core_single;
  regex_t core_pair;
  regex_t ambig;
};

// The state of the parser.
// In a standalone alphabet it will start in STATE_HEADER and stay
// there until it has seen an alphabet header. For an embeded alphabet
// it will jump straight to STATE_CORE. The parser will stay in
// STATE_CORE until it sees the first ambiguous symbol definition at
// which point it will transition to STATE_AMBIG.
typedef enum {
  STATE_ERROR,
  STATE_HEADER,
  STATE_CORE,
  STATE_AMBIG
} STATE_EN;

// Store a CIELAB colour
typedef struct {
  double l;
  double a;
  double b;
} LAB_COLOUR_T;

// Stores information on a symbol
typedef struct {
  char symbol;
  char *aliases;
  int colour;
  char *name;
  char complement;
  char *comprise;
} ALPH_SYM_T;

// state of the alphabet parser
struct alph_reader {
  // parse input into lines
  bool is_done;
  bool had_error;
  bool had_warning;
  bool in_string;
  bool in_comment;
  bool last_is_nl;
  bool last_is_nl_start;
  char last_value;
  int64_t line_num;
  STR_T* line;
  // patterns to match lines
  struct patterns re;
  // parser
  STATE_EN state;
  int version;
  char* alph_name;
  int flags;
  RBTREE_T* core;
  RBTREE_T* all;
  RBTREE_T* aliases;
  RBTREE_T* merged;
  bool seen_lower_case;
  bool seen_upper_case;
  LINKLST_T *messages;
};

// REGULAR EXPRESSION FRAGMENTS
// match a space (POSIX regular expressions require this monstrosity)
#define S "[[:space:]]"
// match a symbol (only alphabetic)
#define SYM "[A-Za-z0-9\\?\\.\\*\\-]"
// match a 6 digit hexadecimal number
#define COLOUR "[0-9A-Fa-f]{6}"
// match quoted text with backslash escapes (should match a JSON string)
// note that this regular expression uses 2 groups so take that into account
// when counting brackets.
#define NAME "\"(([^\\\"]+|\\[\"\\/bfnrt]|\\u[0-9A-Fa-f]{4})*)\""

// REGULAR EXPRESSIONS
// matches the alphabet header
// match[1] = version number
// match[3] = alphabet name (optional)
static const char* HEADER_RE = "^"S"*ALPHABET("S"+v([1-9][0-9]*))?("S"+"NAME")?("S"+(RNA|DNA|PROTEIN)-LIKE)?";
// matches a core symbol without complement
// match[1] = symbol
// match[3] = symbol name (optional)
// match[6] = symbol colour (optional)
static const char* CORE_SINGLE_RE = "^"S"*("SYM")("S"+"NAME")?("S"+("COLOUR"))?"S"*$";
// matches a core symbol and its complementary core symbol
// match[1]  = symbol 1
// match[3]  = symbol 1 name (optional)
// match[6]  = symbol 1 colour (optional)
// match[7]  = symbol 2
// match[9]  = symbol 2 name (optional)
// match[12] = symbol 2 colour (optional)
static const char* CORE_PAIR_RE = "^"S"*("SYM")("S"+"NAME")?("S"+("COLOUR"))?"S"*~"S"*("SYM")("S"+"NAME")?("S"+("COLOUR"))?"S"*$";
// matches a ambiguous symbol and its list of core symbols
// match[1] = ambiguous symbol
// match[3] = ambiguous symbol name (optional)
// match[6] = ambiguous symbol colour (optional)
// match[7] = core symbols
static const char* AMBIG_RE = "^"S"*("SYM")("S"+"NAME")?("S"+("COLOUR"))?"S"*="S"*("SYM"*)"S"*$";

/*
 * Convert a colour in the HSV colourspace into the RGB colourspace.
 * It is easy to create colours in the HSV colourspace but then we need it in
 * the RGB colourspace to actually use.
 */
static int hsv2rgb(double h, double s, double v) {
  double r, g, b, f, p, q, t;
  int i;
  uint8_t red, blue, green;
  int rgb;
  if (s == 0) {
    // achromatic (grey)
    r = v; g = v; b = v;
  } else {
    h /= 60;
    i = floor(h);
    f = h - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * f));
    t = v * (1.0 - (s * (1.0 - f)));
    switch (i) {
      case 0:  r = v; g = t; b = p; break;
      case 1:  r = q; g = v; b = p; break;
      case 2:  r = p; g = v; b = t; break;
      case 3:  r = p; g = q; b = v; break;
      case 4:  r = t; g = p; b = v; break;
      default: r = v; g = p; b = q; break;
    }
  }
  red = (uint8_t)round(r * 255.0);
  green = (uint8_t)round(g * 255.0);
  blue = (uint8_t)round(b * 255.0);

  rgb = (red << 16) | (green << 8) | blue;
  return rgb;
}

/*
 * Convert a colour in the RGB colourspace into the CIELAB colourspace.
 * The CIELAB colourspace will allow us to determine how similar colours are.
 */
static LAB_COLOUR_T rgb2lab(int rgb) {
  LAB_COLOUR_T out;
  double red, green, blue, c1, c2, c3, x, y, z;
  // Split into RGB components
  red = ((rgb >> 16) & 0xFF);
  green = ((rgb >> 8) & 0xFF);
  blue = (rgb & 0xFF);

  // =========================================================================
  // Convert RGB into XYZ colourspace
  // =========================================================================

  // start by changing the range to 0 to 1
  c1 = red / 255.0;
  c2 = green / 255.0;
  c3 = blue / 255.0;

  // adjust values
  if (c1 > 0.04045) {
    c1 = (c1 + 0.055) / 1.055;
    c1 = pow(c1, 2.4);
  } else {
    c1 = c1 / 12.92;
  }
  if (c2 > 0.04045) {
    c2 = (c2 + 0.055) / 1.055;
    c2 = pow(c2, 2.4);
  } else {
    c2 = c2 / 12.92;
  }
  if (c3 > 0.04045) {
    c3 = (c3 + 0.055) / 1.055;
    c3 = pow(c3, 2.4);
  } else {
    c3 = c3 / 12.92;
  }
  // adjust range again
  c1 *= 100.0;
  c2 *= 100.0;
  c3 *= 100.0;
  // apply matrix transform
  x = (c1 * 0.4124) + (c2 * 0.3576) + (c3 * 0.1805);
  y = (c1 * 0.2126) + (c2 * 0.7152) + (c3 * 0.0722);
  z = (c1 * 0.0193) + (c2 * 0.1192) + (c3 * 0.9505);

  // =========================================================================
  // Convert XYZ into Lab colourspace
  // =========================================================================

  c1 = x / 95.047;
  c2 = y / 100.0;
  c3 = z / 108.883;

  if (c1 > 0.008856) {
    c1 = pow(c1, 1.0 / 3.0);
  } else {
    c1 = (7.787 * c1) + (16.0 / 116.0);
  }
  if (c2 > 0.008856) {
    c2 = pow(c2, 1.0 / 3.0);
  } else {
    c2 = (7.787 * c2) + (16.0 / 116.0);
  }
  if (c3 > 0.008856) {
    c3 = pow(c3, 1.0 / 3.0);
  } else {
    c3 = (7.787 * c3) + (16.0 / 116.0);
  }
  out.l = (116.0 * c2) - 16;
  out.a = 500.0 * (c1 - c2);
  out.b = 200.0 * (c2 - c3);
  return out;
}

/*
 * lab_dist
 * Calculate how different two colours look.
 */
static double lab_dist(LAB_COLOUR_T lab1, LAB_COLOUR_T lab2) {
  double c1, c2, dc, dl, da, db, dh, dh_squared, first, second, third;
  c1 = sqrt((lab1.l * lab1.l) + (lab1.a * lab1.a));
  c2 = sqrt((lab2.l * lab2.l) + (lab2.a * lab2.a));
  dc = c1 - c2;
  dl = lab1.l - lab2.l;
  da = lab1.a - lab2.a;
  db = lab1.b - lab2.b;
  // we don't want NaN due to rounding errors so fudge things a bit...
  dh_squared = (da * da) + (db * db) - (dc * dc);
  dh = (dh_squared > 0 ? sqrt(dh_squared) : 0);
  first = dl;
  second = dc / (1.0 + (0.045 * c1));
  third = dh / (1.0 + (0.015 * c1));
  return sqrt((first * first) + (second * second) + (third * third));
}

/*
 * compile_patterns
 * This just compiles the regular expressions used by the parser.
 */
static void compile_patterns(struct patterns *re) {
  regcomp_or_die("header", &(re->header), HEADER_RE, REG_EXTENDED);
  regcomp_or_die("core single", &(re->core_single), CORE_SINGLE_RE, REG_EXTENDED);
  regcomp_or_die("core pair", &(re->core_pair), CORE_PAIR_RE, REG_EXTENDED);
  regcomp_or_die("ambig", &(re->ambig), AMBIG_RE, REG_EXTENDED);
}

/*
 * free_patterns
 * This frees the previously compiled regular expressions
 */
static void free_patterns(struct patterns *re) {
  regfree(&(re->header));
  regfree(&(re->core_single));
  regfree(&(re->core_pair));
  regfree(&(re->ambig));
}

/*
 * alph_sym_create
 * This wraps the passed parameters in a ALPH_SYM_T structure
 * it does NOT copy the name and comprise fields.
 */
static ALPH_SYM_T* alph_sym_create(char symbol, char *name, int32_t colour, char complement, char* comprise) {
  ALPH_SYM_T* sym;
  sym = mm_malloc(sizeof(ALPH_SYM_T));
  sym->symbol = symbol;
  sym->aliases = NULL;
  sym->colour = colour;
  sym->name = name;
  sym->complement = complement;
  sym->comprise = comprise;
  return sym;
}

/*
 * alph_sym_destroy
 * This frees the ALPH_SYM_T structure and the wrapped parameters
 */
static void alph_sym_destroy(void* a_sym) {
  ALPH_SYM_T* sym;
  sym = (ALPH_SYM_T*)a_sym;
  free(sym->name);
  free(sym->comprise);
  free(sym->aliases);
  memset(sym, 0, sizeof(ALPH_SYM_T));
  free(sym);
}

/*
 * alph_sym_p_compare
 * This compares 2 pointers to ALPH_SYM_T*.
 */
static int alph_sym_p_compare(const void *sym1_p, const void *sym2_p) {
  int i, cmp;
  const ALPH_SYM_T *sym1, *sym2;
  sym1 = *((ALPH_SYM_T**)sym1_p);
  sym2 = *((ALPH_SYM_T**)sym2_p);
  // compare comprise set
  if (sym1->comprise != NULL && sym2->comprise != NULL) {
    int sym1_comprise_len, sym2_comprise_len;
    sym1_comprise_len = strlen(sym1->comprise);
    sym2_comprise_len = strlen(sym2->comprise);
    if (sym1_comprise_len == sym2_comprise_len) {
      for (i = 0; i < sym1_comprise_len; i++) {
        cmp = alph_sym_cmp(sym1->comprise+i, sym2->comprise+i);
        if (cmp != 0) return cmp;
      }
    } else if (sym1_comprise_len < sym2_comprise_len) {
      return 1; // ambiguous with longest comprising set goes first
    } else {
      return -1;
    }
  } else if (sym1->comprise != NULL) {
    return 1; // core syms go first
  } else if (sym2->comprise != NULL) {
    return -1;
  }
  // compare symbol
  cmp = alph_sym_cmp(&(sym1->symbol), &(sym2->symbol));
  if (cmp != 0) return cmp;
  // compare complement
  if (sym1->complement != '\0' && sym2->complement != '\0') {
    cmp = alph_sym_cmp(&(sym1->complement), &(sym2->complement));
    if (cmp != 0) return cmp;
  } else if (sym1->complement != '\0') {
    return 1;
  } else if (sym2->complement != '\0') {
    return -1;
  }
  // compare aliases
  if (sym1->aliases != NULL && sym2->aliases != NULL) {
    cmp = strcmp(sym1->aliases, sym2->aliases);
    if (cmp != 0) return cmp;
  }
  // compare name
  if (sym1->name != NULL && sym2->name != NULL) {
    cmp = strcmp(sym1->name, sym2->name);
    if (cmp != 0) return cmp;
  } else if (sym1->name != NULL) {
    return 1;
  } else if (sym2->name != NULL) {
    return -1;
  }
  // compare colour
  return sym1->colour - sym2->colour;
}

/*
 * is_line_empty
 * Return true if the line of text contains only whitespace
 */
static bool is_line_empty(const char *line) {
  int i;
  const char *c;
  for (c = line; *c != '\0'; c++) {
    if (!isspace(*c)) return false;
  }
  return true;
}

/*
 * add_msg
 * Record an error message for reporting to the calling program
 */
static void add_msg(ALPH_READER_T *reader, PARMSG_T *msg) {
  if (msg->severity == SEVERITY_ERROR) {
    reader->had_error = true;
  } else if (msg->severity == SEVERITY_WARNING) {
    reader->had_warning = true;
  }
  linklst_add(msg, reader->messages);
}

/*
 * decode_string
 * Do an in-place decoding of the backslash escaped, UTF-8 encoded string.
 * The backslash escapes will be converted into UTF-8. Any control characters
 * will be removed. The text will be truncated to 30 codepoints.
 * When an empty string is provided it will be freed and NULL returned.
 */
static char* decode_string(char* str, bool *warn_length, bool *warn_disallowed) {
  int32_t codepoint;
  int name_bytes, char_count, offset, dest_offset, code_unit_length, name_length;
  bool keep_letter;
  char hexnum[5];
  *warn_disallowed = false;
  *warn_length = false;
  if (str != NULL && str[0] != '\0') {
    // need to parse as UTF-8 with backslash escapes to get correct string
    // length and remove disallowed characters
    name_bytes = strlen(str);
    name_length = 0;
    offset = 0;
    dest_offset = 0;
    while (offset < name_bytes) {
      keep_letter = true;
      // first convert to a codepoint
      codepoint = unicode_from_string(str+offset, name_bytes - offset, &char_count);
      if (codepoint < 0) die("Encoding of alphabet is not valid UTF-8");
      if (str[offset] == '\\') {
        // the regular expression should ensure this is a valid escape sequence
        if ((offset + 1) >= name_bytes) die("Not enough characters to handle escape sequence");
        switch(str[offset+1]) {
          // we accept these normal characters
          case '"':
          case '\\':
          case '/':
            codepoint = str[offset+1] - 0x0;
            char_count = 2;
            break;
          // these are all control characters so we reject them
          case 'b':
          case 'f':
          case 'n':
          case 'r':
          case 't':
            char_count = 2;
            keep_letter = false;
            break;
          case 'u':
            // this unicode escape could be anything so we convert it to a codepoint
            if ((offset + 5) >= name_bytes) die("Not enough characters to handle escape sequence");
            char_count = 6;
            strncpy(hexnum, str+(offset+2), 4);
            hexnum[4] = '\0';
            codepoint = strtol(hexnum, NULL, 16);
            break;
        }
      }
      if ((codepoint >= 0x0 && codepoint <= 0x1F) || codepoint == 0x7F || // ASCII range control characters (aka C0)
          (codepoint >= 0x80 && codepoint <= 0x9F) || // C1
          codepoint == 0x2028 || codepoint == 0x2029 || // line separator and paragraph separator
          codepoint == 0x200E || codepoint == 0x200F || (codepoint >= 0x202A && codepoint <= 0x202E) // Bidirectional text control
      ) {
        // these characters are not allowed because they could stuff things up...
        keep_letter = false;
      }
      if (keep_letter) {
        if (name_length < 40) {
          unicode_to_string(codepoint, str+dest_offset, &code_unit_length);
          dest_offset += code_unit_length;
          name_length++;
        } else {
          *warn_length = true;
        }
      } else {
        *warn_disallowed = true;
      }
      offset += char_count;
    }
    str[dest_offset] = '\0';
    if (dest_offset < name_bytes) {
      str = mm_realloc(str, sizeof(char) * (dest_offset + 1));
    }
  }
  // when the str is empty, set it to null
  if (str != NULL && str[0] == '\0') {
    free(str);
    str = NULL;
  }
  return str;
}

static char* decode_name(ALPH_READER_T *reader, const int64_t line_num, char symbol, char* name) {
  bool warn_length, warn_disallowed;
  char *decoded;
  decoded = decode_string(name, &warn_length, &warn_disallowed);
  if (warn_length && warn_disallowed) {
    add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, line_num, -1,
          "name of %c contained disallowed characters and was too long", symbol));
  } else if (warn_length) {
    add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, line_num, -1,
        "name of %c was too long", symbol));
  } else if (warn_disallowed) {
    add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, line_num, -1,
        "name of %c contained disallowed characters", symbol));
  }
  return decoded;
}
/*
 * Group the symbol by the comprising characters
 */
static void track_alias(ALPH_READER_T *reader, ALPH_SYM_T *sym) {
  RBNODE_T *aliases_node;
  bool created;
  ARRAYLST_T *aliases_list;
  char *comprise, buffer[2];
  if (sym->comprise) {
    comprise = sym->comprise;
  } else {
    buffer[0] = sym->symbol;
    buffer[1] = '\0';
    comprise = buffer;
  }
  aliases_node = rbtree_lookup(reader->aliases, comprise, true, &created);
  if (created) {
    aliases_list = arraylst_create();
    rbtree_set(reader->aliases, aliases_node, aliases_list);
  } else {
    aliases_list = (ARRAYLST_T*)rbtree_value(aliases_node);
  }
  arraylst_add(sym, aliases_list);
}

/*
 * process_header
 * Store the version
 */
static void process_header(ALPH_READER_T *reader, const int64_t line_num, int version, char *name, int flags) {
  reader->version = version;
  reader->alph_name = name;
  reader->flags = flags;
}

/*
 * process_core
 * check that a core symbol is new and add it to the set of core symbols along
 * with all its data.
 */
static void process_core(ALPH_READER_T *reader, const int64_t line_num,
    char symbol, char* name, int colour, char complement) {
  ALPH_SYM_T *sym;
  char sym_str[2];
  if (symbol == '?') {
    add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, line_num, -1,
          "symbol '?' is reserved as a wildcard and can not be defined to have any other meaning"));
  }
  sym = alph_sym_create(symbol, name, colour, complement, NULL);
  sym_str[0] = symbol;
  sym_str[1] = '\0';
  if (!rbtree_make(reader->all, sym_str, NULL)) {
    alph_sym_destroy(sym);
    add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, line_num, -1,
          "core symbol %c is already defined", symbol));
    return;
  }
  track_alias(reader, sym);
  rbtree_make(reader->core, sym_str, sym); // track the set of core symbols
  if (symbol >= 'A' && symbol <= 'Z') reader->seen_upper_case = true;
  else if (symbol >= 'a' && symbol <= 'z') reader->seen_lower_case = true;
}

/*
 * sort_and_remove_duplicates
 * Ensure symbols are sorted and there are no duplicates.
 */
static bool sort_and_remove_duplicates(char *syms) {
  int nsyms;
  char *a, *b;
  bool duplicate;
  nsyms = strlen(syms);
  if (nsyms <= 1) return false;
  qsort(syms, nsyms, sizeof(char), alph_sym_cmp);
  duplicate = false;
  a = syms;
  b = syms+1;
  while (*b != '\0') {
    if (*a != *b) {
      a++;
      *a = *b;
      b++;
    } else {
      duplicate = true;
      b++;
    }
  }
  a++;
  *a = '\0';
  return duplicate;
}

/*
 * process_ambig
 * process a ambiguous symbol and the list of equal core symbols
 * check that the referenced core symbols have been defined and
 * check that the ambiguous symbol is completely new. Add the
 * ambiguous symbol to the set of ambiguous symbols.
 */
static void process_ambig(ALPH_READER_T *reader, const int64_t line_num,
    char symbol, char* name, int colour, char *comprise) {
  ALPH_SYM_T *sym;
  char *eq_sym, sym_str[2];
  // Ensure comprise is sorted and there are no duplicates.
  if (sort_and_remove_duplicates(comprise)) {
    add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, line_num, -1,
        "ambiguous symbol %c has a duplication in the comprising symbols", symbol));
  }
  // check that the referenced symbols exist in the core set
  sym_str[1] = '\0';
  for (eq_sym = comprise; *eq_sym != '\0'; eq_sym++) {
    sym_str[0] = *eq_sym;
    if (rbtree_find(reader->core, sym_str) == NULL) {
      add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, line_num, -1,
            "required core symbol %c has not been defined", *eq_sym));
    }
  }
  if (symbol == '?') {
    // check it is a wildcard
    if (strlen(comprise) != rbtree_size(reader->core)) {
      add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, line_num, -1,
            "symbol '?' is reserved as a wildcard and can not be defined to have any other meaning"));
    }
  }
  sym_str[0] = symbol;
  sym = alph_sym_create(symbol, name, colour, '\0', comprise);
  if (!rbtree_make(reader->all, sym_str, NULL)) {
    alph_sym_destroy(sym);
    add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, line_num, -1,
          "ambiguous symbol %c is already defined", symbol));
    return;
  }
  track_alias(reader, sym);
  if (symbol >= 'A' && symbol <= 'Z') reader->seen_upper_case = true;
  else if (symbol >= 'a' && symbol <= 'z') reader->seen_lower_case = true;
}

/*
 * check_complements
 * Ensure that all referenced complement symbols actually exist.
 */
static void check_complements(ALPH_READER_T *reader) {
  RBNODE_T *node;
  ALPH_SYM_T *sym;
  char sym_str[2];
  sym_str[1] = '\0';
  for (node = rbtree_first(reader->core); node != NULL; node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    if (sym->complement == '\0') continue; // no complement
    sym_str[0] = sym->complement;
    if (rbtree_find(reader->core, sym_str) == NULL) {
      add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, reader->line_num, -1, 
            "core symbol %c has complement %c which has not been defined",
            sym->symbol, sym->complement));
    }
  }
}

/*
 * cleanup_alias_list
 * Frees memory used by an alias list. Does not touch the contained symbols as
 * they are assumed to be managed elsewhere.
 */
static void cleanup_alias_list(void* list) {
  // the referenced symbols will be destroyed
  // with the core and ambigs rbtree
  arraylst_destroy(alph_sym_destroy, (ARRAYLST_T*)list);
}

/*
 * alph_reader_create
 * Initialize the ALPH_READER_T
 */
ALPH_READER_T* alph_reader_create() {
  ALPH_READER_T *reader;
  reader = mm_malloc(sizeof(ALPH_READER_T));
  memset(reader, 0, sizeof(ALPH_READER_T));
  compile_patterns(&(reader->re));
  reader->is_done = false;
  reader->had_error = false;
  reader->had_warning = false;
  reader->in_string = false;
  reader->in_comment = false;
  reader->last_is_nl = false;
  reader->last_is_nl_start = false;
  reader->last_value = '\0';
  reader->line_num = 0;
  reader->line = str_create(100);
  reader->state = STATE_HEADER;
  reader->version = 0;
  reader->alph_name = NULL;
  reader->core = rbtree_create(alph_str_cmp, rbtree_strcpy, free, NULL, NULL);
  reader->all = rbtree_create(alph_str_cmp, rbtree_strcpy, free, NULL, NULL);
  reader->aliases = rbtree_create(alph_str_cmp, rbtree_strcpy, free, NULL, cleanup_alias_list);
  reader->merged = rbtree_create(alph_str_cmp, rbtree_strcpy, free, NULL, alph_sym_destroy);
  reader->seen_lower_case = false;
  reader->seen_upper_case = false;
  reader->messages = linklst_create();
  return reader;
}

/*
 * alph_reader_destroy
 * Destroy the ALPH_READER_T
 */
void alph_reader_destroy(ALPH_READER_T *reader) {
  linklst_destroy_all(reader->messages, parmsg_destroy);
  rbtree_destroy(reader->all);
  rbtree_destroy(reader->core);
  rbtree_destroy(reader->aliases);
  rbtree_destroy(reader->merged);
  free(reader->alph_name);
  str_destroy(reader->line, false);
  free_patterns(&(reader->re));
  memset(reader, 0, sizeof(ALPH_READER_T));
  free(reader);
}

/*
 * alph_reader_header
 * Adds the header bypassing the text processing.
 */
void alph_reader_header(ALPH_READER_T *reader, int version, const char *name, int flags) {
  if (reader->state != STATE_HEADER) die("Alphabet header must be specified first!");
  process_header(reader, -1, version, (name == NULL ? NULL : strdup(name)), flags);
  reader->state = STATE_CORE;
}

/*
 * alph_reader_core
 * Adds a core symbol bypassing the text processing.
 */
void alph_reader_core(ALPH_READER_T *reader, char symbol, const char* aliases, const char* name, int colour, char complement) {
  if (reader->state != STATE_CORE) die("Alphabet core symbols must be specified before ambiguous symbols!");
  process_core(reader, -1, symbol, (name == NULL ? NULL : strdup(name)), colour, complement);
  if (aliases != NULL) {
    char comprise[2];
    int i;
    comprise[0] = symbol;
    comprise[1] = '\0';
    for (i = 0; aliases[i] != '\0'; i++) process_ambig(reader, -1, aliases[i], NULL, -1, strdup(comprise));
  }
}

/*
 * alph_reader_ambig
 * Adds a ambiguous symbol bypassing the text processing.
 */
void alph_reader_ambig(ALPH_READER_T *reader, char symbol, const char* aliases, const char* name, int colour, const char *comprise) {
  if (reader->state != STATE_CORE && reader->state != STATE_AMBIG) die("Alphabet header must be specified first!");
  if (reader->state == STATE_CORE) check_complements(reader);
  process_ambig(reader, -1, symbol, (name == NULL ? NULL : strdup(name)), colour, strdup(comprise));
  reader->state = STATE_AMBIG;
  if (aliases != NULL) {
    int i;
    for (i = 0; aliases[i] != '\0'; i++) process_ambig(reader, -1, aliases[i], NULL, -1, strdup(comprise));
  }
}

/*
 * merge_aliases
 * This merges the fields of the aliases and stores it in the merged tree.
 */
static void merge_aliases(ALPH_READER_T *reader) {
  ALPH_SYM_T *sym;
  RBNODE_T *node, *node2;
  ARRAYLST_T *aliases;
  RBTREE_T *alts;
  int i, j, k;
  alts = rbtree_create(rbtree_charcmp, rbtree_charcpy, free, NULL, NULL);
  for (i = 0, node = rbtree_first(reader->aliases); node != NULL; node = rbtree_next(node), i++) {
    ALPH_SYM_T *sym, *alias;
    if (i != 0) rbtree_clear(alts); // reset the alts set
    sym = mm_malloc(sizeof(ALPH_SYM_T));
    aliases = (ARRAYLST_T*)rbtree_value(node);
    arraylst_qsort(alph_sym_p_compare, aliases);
    alias = (ALPH_SYM_T*)arraylst_get(0, aliases);
    sym->symbol = alias->symbol;
    sym->complement = alias->complement;
    sym->comprise = (alias->comprise != NULL ? strdup(alias->comprise) : NULL);
    sym->name = NULL;
    sym->colour = -1;
    sym->aliases = NULL;
    for (j = 0; j < arraylst_size(aliases); j++) {
      alias = (ALPH_SYM_T*)arraylst_get(j, aliases);
      rbtree_make(alts, &(alias->symbol), NULL);
      if (alias->name != NULL && sym->name == NULL) sym->name = strdup(alias->name);
      if (alias->colour != -1 && sym->colour == -1) sym->colour = alias->colour;
      if (alias->aliases != NULL) {
        for (k = 0; alias->aliases[k]; k++) {
          rbtree_make(alts, &(alias->aliases[k]), NULL);
        }
      }
    }
    // remove the main symbol from the alternate symbol set
    rbtree_remove(alts, &(sym->symbol));
    // flatten the tree set into a list of characters.
    sym->aliases = mm_malloc(sizeof(char) * (rbtree_size(alts) + 1));
    for (j = 0, node2 = rbtree_first(alts); node2 != NULL; node2 = rbtree_next(node2), j++) {
      sym->aliases[j] = *((char*)rbtree_key(node2));
    }
    sym->aliases[j] = '\0';
    // sort the alt symbols
    qsort(sym->aliases, rbtree_size(alts), sizeof(char), alph_sym_cmp);
    // finally put the merged symbol into the lookup
    rbtree_put(reader->merged, rbtree_key(node), sym);
  }
  rbtree_destroy(alts);
}

/*
 * set_missing_colours
 * This assigns colours to the symbols that don't have them.
 */
static void set_missing_colours(ALPH_READER_T *reader) {
  RBTREE_T *unique_colours;
  RBNODE_T *node;
  ALPH_SYM_T *sym, *sym2;
  int ncolours, *colours, best, i, j;
  LAB_COLOUR_T *lab_colours, lab_colour;
  double hue, sat, light, step, min, value;
  // 1) Count the number of unique colours already assigned and
  //    the number of symbols with no assigned colour. If a ambiguous symbol
  //    does not have an assigned colour then assume the colour is black and
  //    exclude it from the count. This total will be referred to as X.
  unique_colours = rbtree_create(rbtree_intcmp, rbtree_intcpy, free, NULL, NULL);
  ncolours = 0;
  for (node = rbtree_first(reader->merged); node != NULL; node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    if (sym->colour == -1) {
      if (sym->comprise == NULL) ncolours++; // count core symbols without colour
    } else {
      // track all unique colours
      if (rbtree_make(unique_colours, &(sym->colour), NULL)) ncolours++;
    }
  }
  // 2) Generate X evenly distributed colours in the HSV colourspace and convert to RGB and Lab.
  hue = 0; // red
  sat = 1.0;
  light = 0.4;
  step = 360.0 / ncolours;
  colours = mm_malloc(sizeof(int) * ncolours);
  lab_colours = mm_malloc(sizeof(LAB_COLOUR_T) * ncolours);
  for (i = 0; i < ncolours; i++) {
    colours[i] = hsv2rgb((hue + (step * i)), sat, light);
    lab_colours[i] = rgb2lab(colours[i]);
  }
  // 3) For each of the unique colours find the closest colour in the X
  //    generated colours and remove it.
  for (node = rbtree_first(unique_colours); node != NULL; node = rbtree_next(node)) {
    lab_colour = rgb2lab(*((int*)rbtree_key(node)));
    min = DBL_MAX;
    best = -1;
    for (i = 0; i < ncolours; i++) {
      if (colours[i] == -1) continue;
      value = lab_dist(lab_colour, lab_colours[i]);
      if (value < min) {
        min = value;
        best = i;
      }
    }
    assert(best != -1);
    colours[best] = -1;
  }
  // done with the lab colorspace
  free(lab_colours);
  rbtree_destroy(unique_colours);
  // 4) Assign colours to any symbols without
  for (i = 0, node = rbtree_first(reader->merged); node != NULL; node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    if (sym->colour == -1) {
      // symbol does not have a colour
      if (sym->comprise == NULL) {
        // assign one of the remaining colours to the core symbol
        while (colours[i] == -1) i++;
        sym->colour = colours[i];
        i++;
      } else {
        // assign black to the ambiguous symbol
        sym->colour = 0; // black in RGB
      }
    }
  }
  // finished with generated colours
  free(colours);
}

/*
 * Tries to find a complement for every ambiguous character
 */
static void set_missing_complements(ALPH_READER_T *reader) {
  ALPH_SYM_T *sym, *core, *complement;
  char *opposite, s[2];
  RBNODE_T *node;
  int size, i;
  s[1] = '\0';
  // allocate enough space for opposite to fit all core letters
  size = rbtree_size(reader->core);
  opposite = mm_malloc(sizeof(char) * (size + 1));
  opposite[size] = '\0';
  for (node = rbtree_first(reader->merged); node != NULL; node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    // skip over core symbols
    if (sym->comprise == NULL) continue;
    // skip over symbols with a known complement
    if (sym->complement != '\0') continue;
    // complement the comprising symbols
    for (i = 0; sym->comprise[i]; i++) {
      s[0] = sym->comprise[i];
      core = (ALPH_SYM_T*)rbtree_get(reader->merged, s);
      if (!core->complement) goto next_ambig;
      opposite[i] = core->complement;
    }
    opposite[i] = '\0';
    // sort the comprising symbols of the complement
    sort_and_remove_duplicates(opposite);
    // lookup the complement
    complement = (ALPH_SYM_T*)rbtree_get(reader->merged, opposite);
    if (complement != NULL) sym->complement = complement->symbol;
next_ambig: ;
  }
  free(opposite);
}

/*
 * create_wildcard_if_missing
 * Alphabets must have a wildcard so this creates one using
 * the reserved character '?' if none is present.
 */
static void create_wildcard_if_missing(ALPH_READER_T *reader) {
  RBNODE_T *node;
  ALPH_SYM_T *sym;
  char *core_syms;
  int ncore, i;
  bool created;
  // allocate enough space to fit all the core letters
  ncore = rbtree_size(reader->core);
  core_syms = mm_malloc(sizeof(char) * (ncore + 1));
  core_syms[ncore] = '\0';
  // create string of all core symbols
  for (node = rbtree_first(reader->merged), i = 0; node != NULL; node = rbtree_next(node), i++) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    if (sym->comprise != NULL) break;
    assert(i < ncore);
    core_syms[i] = sym->symbol;
  }
  assert(i == ncore);
  // lookup/create the wildcard
  node = rbtree_lookup(reader->merged, core_syms, true, &created);
  if (created) {
    sym = alph_sym_create('?', NULL, -1, '\0', core_syms);
    rbtree_set(reader->merged, node, sym);
    //add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, -1, -1, "wildcard symbol automatically generated"));
  } else {
    free(core_syms);
  }
}

static void calculate_case_insensitivity(ALPH_READER_T *reader) {
  RBNODE_T *node, *node2;
  RBTREE_T *alias_set;
  ALPH_SYM_T *sym;
  char *alias, symbol;
  bool case_insensitive;
  case_insensitive = true;
  alias_set = rbtree_create(rbtree_charcmp, NULL, NULL, NULL, NULL);
  for (node = rbtree_first(reader->merged); node != NULL; node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    // add all aliases to the set
    rbtree_make(alias_set, &(sym->symbol), NULL);
    for (alias = sym->aliases; *alias != '\0'; alias++) {
      rbtree_make(alias_set, alias, NULL);
    }
    // check if all are case insensitive
    for (node2 = rbtree_first(alias_set); node2 != NULL; node2 = rbtree_next(node2)) {
      symbol = *((char*)rbtree_value(node));
      if (symbol >= 'A' && symbol <= 'Z') {
        symbol += 32;
        if (rbtree_find(alias_set, &symbol) == NULL) {
          case_insensitive = false;
          break;
        }
      } else if (symbol >= 'a' && symbol <= 'z') {
        symbol -= 32;
        if (rbtree_find(alias_set, &symbol) == NULL) {
          case_insensitive = false;
          break;
        }
      }
    }
    if (!case_insensitive) break;
    rbtree_clear(alias_set);
  }
  rbtree_destroy(alias_set);
  if (case_insensitive) reader->flags |= ALPH_CASE_INSENSITIVE;
}

/*
 * check_alphabet_extension
 * If the alphabet is flagged as extending one of the standard alphabets then
 * confirm that it contains all the core prime symbols and their complements.
 */
static void check_alphabet_extension(ALPH_READER_T *reader) {
  char *ext_name, *req_syms, *req_comp;
  char lookup[2];
  ALPH_SYM_T *sym;
  int i;
  lookup[1] = '\0';
  req_syms = NULL;
  req_comp = NULL;
  switch (reader->flags & ALPH_FLAG_EXTENDS) {
    case ALPH_FLAG_EXTENDS_RNA:
      ext_name = "RNA";
      req_syms = "ACGU";
      break;
    case ALPH_FLAG_EXTENDS_DNA:
      ext_name = "DNA";
      req_syms = "ACGT";
      req_comp = "TGCA";
      break;
    case ALPH_FLAG_EXTENDS_PROTEIN:
      ext_name = "protein";
      req_syms = "ACDEFGHIKLMNPQRSTVWY";
      break;
  }
  if (req_syms != NULL) {
    for (i = 0; req_syms[i] != '\0'; i++) {
      lookup[0] = req_syms[i];
      sym = (ALPH_SYM_T*)rbtree_get(reader->merged, lookup);
      if (sym == NULL) {
        add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, -1, -1, 
              "not like %s alphabet as %c is not defined",
              ext_name, req_syms[i]));
      } else {
        ALPH_SYM_T *comp1 = NULL, *comp2 = NULL;
        if (req_comp != NULL) {
          lookup[0] = req_comp[i];
          comp1 = (ALPH_SYM_T*)rbtree_get(reader->merged, lookup);
        }
        if (sym->complement != '\0') {
          lookup[0] = sym->complement;
          comp2 = (ALPH_SYM_T*)rbtree_get(reader->merged, lookup);
        }
        if (comp1 && (comp1 != comp2)) {
          add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, -1, -1, 
                "not like %s alphabet as %c complement rules are incorrect",
                ext_name, req_syms[i]));
        }
      }
    }
  }
}

/*
 * alph_reader_done
 * Finalize the input.
 */
void alph_reader_done(ALPH_READER_T *reader) {
  if (!reader->is_done) {
    if (reader->state != STATE_AMBIG) check_complements(reader);
    merge_aliases(reader);
    create_wildcard_if_missing(reader);
    set_missing_colours(reader);
    set_missing_complements(reader);
    calculate_case_insensitivity(reader);
    check_alphabet_extension(reader);
  }
  reader->is_done = true;
}

/*
 * alph_reader_line
 * Process a line of the file.
 */
void alph_reader_line(ALPH_READER_T *reader, const char *line) {
  regmatch_t matches[13];
  // check the line isn't just empty
  if (is_line_empty(line)) {
    reader->line_num++;
    return;
  }
  switch (reader->state) {
    case STATE_HEADER:
      if (regexec_or_die("header", &(reader->re.header), line, 8, matches, 0)) {
        bool warn_length, warn_disallowed;
        char *name;
        int flags;
        name = decode_string(regex_str(matches+4, line), &warn_length, &warn_disallowed);
        if (warn_length && warn_disallowed) {
          add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, reader->line_num, -1,
                "name of alphabet contained disallowed characters and was too long"));
        } else if (warn_length) {
          add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, reader->line_num, -1,
              "name of alphabet was too long"));
        } else if (warn_disallowed) {
          add_msg(reader, parmsg_create(SEVERITY_WARNING, -1, reader->line_num, -1,
              "name of alphabet contained disallowed characters"));
        }
        flags = 0;
        if (regex_cmp(matches+7, line, "RNA") == 0) {
          flags = ALPH_FLAG_EXTENDS_RNA;
        } else if (regex_cmp(matches+7, line, "DNA") == 0) {
          flags = ALPH_FLAG_EXTENDS_DNA;
        } else if (regex_cmp(matches+7, line, "PROTEIN") == 0) {
          flags = ALPH_FLAG_EXTENDS_PROTEIN;
        }
        process_header(reader, reader->line_num, regex_int(matches+2, line, 1), name, flags);
        reader->state = STATE_CORE;
      } else {
        add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, reader->line_num, -1,
              "no header line"));
        reader->state = STATE_ERROR;
      }
      break;
    case STATE_CORE:
      if (regexec_or_die("core single", &(reader->re.core_single), line, 7, matches, 0)) {
        char symbol;
        symbol = regex_chr(matches+1, line);
        process_core(reader, reader->line_num, symbol,
            decode_name(reader, reader->line_num, symbol, regex_str(matches+3, line)),
            regex_int_with_base(matches+6, line, 16, -1),
            '\0');
        break;
      } else if (regexec_or_die("core pair", &(reader->re.core_pair), line, 13, matches, 0)) {
        char symbol1, symbol2;
        symbol1 = regex_chr(matches+1, line);
        symbol2 = regex_chr(matches+7, line);
        // process symbol 1
        process_core(reader, reader->line_num, symbol1,
            decode_name(reader, reader->line_num, symbol1, regex_str(matches+3, line)),
            regex_int_with_base(matches+6, line, 16, -1),
            symbol2);
        // process symbol 2
        process_core(reader, reader->line_num, symbol2,
            decode_name(reader, reader->line_num, symbol2, regex_str(matches+10, line)),
            regex_int_with_base(matches+12, line, 16, -1),
            symbol1);
        break;
      }
      // fall through
    case STATE_AMBIG:
      if (regexec_or_die("ambig", &(reader->re.ambig), line, 8, matches, 0)) {
        char symbol;
        if (reader->state == STATE_CORE) check_complements(reader);
        symbol = regex_chr(matches+1, line);
        process_ambig(reader, reader->line_num, symbol,
            decode_name(reader, reader->line_num, symbol, regex_str(matches+3, line)),
            regex_int_with_base(matches+6, line, 16, -1),
            regex_str(matches+7, line));
        reader->state = STATE_AMBIG; // don't allow further core entries
      } else {
        add_msg(reader, parmsg_create(SEVERITY_ERROR, -1, reader->line_num, -1,
              "line \"%s\" did not match any expected patterns", line));
      }
      break;
    case STATE_ERROR:
      break;
  }
  reader->line_num++;
}

/*
 * alph_reader_update
 * Update the reader with a chunk of text.
 * This will dispatch the text as lines caching any incomplete lines.
 * Comments will be removed.
 */
void alph_reader_update(ALPH_READER_T *reader, const char *chunk, size_t size, bool end) {
  char value;
  bool is_nl, is_nl_start, is_comment_start;
  int i, offset;
  if (reader->is_done) die("Alphabet reader was already told the input was ended");
  // split into lines
  offset = 0;
  for (i = 0; i < size; i++) {
    is_nl_start = false;
    is_comment_start = false;
    value = chunk[i];
    is_nl = (value == '\n' || value == '\r');
    // check if we've already seen a newline
    if (reader->last_is_nl) {
      // see what the current char is
      if (!is_nl) {
        // not a newline, so line is incremented by previous newline
        reader->in_comment = false;
        reader->in_string = false; // strings are not multi-line
      } else {
        // current char is a newline, so line is only incremented if 
        // the previous is identical char or it is not a newline beginning
        if (reader->last_value == value || !reader->last_is_nl_start) {
          // previous char was a differnt newline so increment line and reset
          reader->in_comment = false;
          reader->in_string = false; // strings are not multi-line
          // so by deduction the current char is a newline beginning
          is_nl_start = true;
        }
      }
    } else {
      // previous char is not a newline
      if (is_nl) {
        // so by deduction the current char is a newline beginning
        is_nl_start = true;
      }
    }
    // check to see if we're in a string or comment
    if (!reader->in_comment) {
      if (reader->in_string) {
        if (value == '"' && reader->last_value != '\\') {
          reader->in_string = false;
        }
      } else {
        if (value == '"') {
          reader->in_string = true;
        } else if (value == '#') {
          reader->in_comment = true;
          // send the line up to this point
          is_comment_start = true;
        }
      }
    }
    // update the line parsing variables
    reader->last_value = value;
    reader->last_is_nl = is_nl;
    reader->last_is_nl_start = is_nl_start;

    // now process the lines (with comments removed)
    if ((is_nl_start && !reader->in_comment) || is_comment_start) {
      str_append(reader->line, chunk+offset, (i - offset));
      alph_reader_line(reader, str_internal(reader->line));
      str_clear(reader->line);
    }
    // handle multiple lines in one chunk
    if (is_nl) {
      offset = i + 1;
    }
  }
  // store any left over text
  if (!reader->in_comment) {
    if (end) {
      str_append(reader->line, chunk+offset, (size - offset));
      alph_reader_line(reader, str_internal(reader->line));
      str_clear(reader->line);
    } else {
      str_append(reader->line, chunk+offset, (size - offset));
    }
  }
  if (end) alph_reader_done(reader);
}

/*
 * add_encoding
 * Adds a letter to an encoding array. If the alphabet is is specified in
 * all upper case or all lower case then we accept the opposite case letters,
 * this does not effect numbers or punctuation symbols.
 */
static inline void add_encoding(uint8_t* encode, bool both_case, char symbol, uint8_t index) {
  encode[(uint8_t)symbol] = index;
  if (both_case) {
    if (symbol >= 'A' && symbol <= 'Z') { // convert to lower case
      encode[(uint8_t)symbol + 32] = index;
    } else if (symbol >= 'a' && symbol <= 'z') { // convert to upper case
      encode[(uint8_t)symbol - 32] = index;
    } // else not alphabetic and hence no other case option
  }
}

/*
 * add_encoding
 * Adds a letter and the alternates to an encoding array. 
 * See add_encoding for more details.
 */
static inline void add_encodings(uint8_t* encode, bool both_case, char symbol, char *alts, uint8_t index) {
  int i;
  add_encoding(encode, both_case, symbol, index);
  if (alts != NULL) {
    for (i = 0; alts[i] != '\0'; i++) {
      add_encoding(encode, both_case, alts[i], index);
    }
  }
}

/*
 * alph_reader_alphabet
 * Get the loaded alphabet.
 * This will return NULL if the parser is not done or there are errors.
 */
ALPH_T* alph_reader_alphabet(ALPH_READER_T *reader) {
  ALPH_T *alphabet, *reference;
  RBNODE_T *node;
  ALPH_SYM_T *sym;
  char symbol[2];
  char *temp_complements;
  int i, j, nfull, ncomprise;
  bool both_case;
  // check that we have all the information to do this
  if (!reader->is_done || reader->had_error) return NULL;
  // see if we should treat uppercase and lowercase the same
  both_case = reader->seen_lower_case != reader->seen_upper_case;
  // count the total number of symbols in the alphabet
  nfull = rbtree_size(reader->merged);
  // initilise tempory storage for complement symbols
  temp_complements = mm_malloc(sizeof(char) * (nfull + 1));
  // allocate memory for a new alphabet
  alphabet = mm_malloc(sizeof(ALPH_T));
  memset(alphabet, 0, sizeof(ALPH_T));
  // initially there is one owner
  alphabet->references = 1;
  // set flags
  alphabet->flags = reader->flags;
  // set the number of core symbols
  alphabet->ncore = rbtree_size(reader->core);
  // set the total number of symbols
  alphabet->nfull = nfull;
  // Now allocate some memory for fields. Note that we always allocate an 
  // extra position for the error state at index 0.
  // Allocate symbols (allocate extra for null byte as well as error state)
  alphabet->symbols = mm_malloc(sizeof(char) * (nfull + 2));
  // Allocate aliases
  alphabet->aliases = mm_malloc(sizeof(char*) * (nfull + 1));
  // Allocate names.
  alphabet->names = mm_malloc(sizeof(char*) * (nfull + 1));
  // Allocate colours.
  alphabet->colours = mm_malloc(sizeof(uint32_t) * (nfull + 1));
  // Allocate comprising symbol counts.
  alphabet->ncomprise = mm_malloc(sizeof(uint8_t) * (nfull + 1));
  // allocate comprising symbol lists
  alphabet->comprise = mm_malloc(sizeof(uint8_t*) * (nfull + 1));
  // Allocate complements
  alphabet->complement = mm_malloc(sizeof(uint8_t) * (nfull + 1));
  // set the values for the ERROR state
  alphabet->symbols[0] = '?';
  alphabet->aliases[0] = strdup("");
  alphabet->names[0] = strdup("Unknown");
  alphabet->colours[0] = 0;
  alphabet->ncomprise[0] = 0;
  alphabet->comprise[0] = mm_malloc(sizeof(uint8_t));
  alphabet->comprise[0][0] = 0;
  alphabet->complement[0] = 0;
  // initialize the mapping
  for (i = 0; i <= UCHAR_MAX; i++) {
    alphabet->encode[i] = 0;
    alphabet->encode2core[i] = 0;
  }
  // now iterate over the core symbols
  symbol[1] = '\0';
  for (i = 1, node = rbtree_first(reader->merged); node != NULL; i++, node = rbtree_next(node)) {
    sym = (ALPH_SYM_T*)rbtree_value(node);
    alphabet->symbols[i] = sym->symbol;
    alphabet->aliases[i] = strdup(sym->aliases != NULL ? sym->aliases : "");
    // name
    symbol[0] = sym->symbol;
    alphabet->names[i] = strdup(sym->name != NULL ? sym->name : symbol);
    // colour
    alphabet->colours[i] = sym->colour;
    // assign the mapping of symbol to index
    add_encodings(alphabet->encode, both_case, sym->symbol, sym->aliases, i);
    if (sym->comprise == NULL) {
      // assign the core only mapping of symbol to index
      add_encodings(alphabet->encode2core, both_case, sym->symbol, sym->aliases, i);
      // core symbols comprise themselves only
      alphabet->ncomprise[i] = 1;
      alphabet->comprise[i] = mm_malloc(sizeof(uint8_t) * 2);
      alphabet->comprise[i][0] = i;
      alphabet->comprise[i][1] = 0; // zero terminate comprising symbols
    } else {
      // count the number of comprising symbols
      ncomprise = strlen(sym->comprise);
      alphabet->ncomprise[i] = ncomprise;
      alphabet->comprise[i] = mm_malloc(sizeof(uint8_t) * (ncomprise + 1));
      for (j = 0; j < ncomprise; j++) {
        alphabet->comprise[i][j] = alphabet->encode2core[(uint8_t)sym->comprise[j]];
      }
      alphabet->comprise[i][j] = 0; // zero terminate comprising symbols
    }
    // can't do the complement in this loop because we're still creating the encode array
    temp_complements[i] = sym->complement;
  }
  // null terminate symbols
  alphabet->symbols[i] = '\0';
  // set the complements
  for (i = 1; i <= nfull; i++) {
    alphabet->complement[i] = alphabet->encode[(uint8_t)temp_complements[i]];
  }
  free(temp_complements);
  // set alphabet name
  if (reader->alph_name == NULL || strlen(reader->alph_name) == 0) {
    // when no name set, copy the core symbols as the alphabet name
    alphabet->alph_name = mm_malloc(sizeof(char) * (alphabet->ncore + 1));
    strncpy(alphabet->alph_name, alphabet->symbols+1, alphabet->ncore);
    alphabet->alph_name[alphabet->ncore] = '\0';
  } else {
    alphabet->alph_name = strdup(reader->alph_name);
  }
  // create safe encoding arrays
  for (i = 0; i <= UCHAR_MAX; i++) {
    alphabet->encodesafe[i] = (alphabet->encode[i] ? alphabet->encode[i] : alphabet->ncore + 1) - 1;
    alphabet->encodesafe2core[i] = (alphabet->encode2core[i] ? alphabet->encode2core[i] : alphabet->ncore + 1) - 1;
  }
  // finally if the built-in flag is not set but the alphabet is
  // extending one of the built-in alphabets check if it actually
  // is one of the built in alphabets.
  if (!alph_is_builtin(alphabet) && alph_extends(alphabet)) {
    if (alph_extends_dna(alphabet)) {
      reference = alph_dna();
    } else if (alph_extends_rna(alphabet)) {
      reference = alph_rna();
    } else {
      reference = alph_protein();
    }
    if (alph_equal(reference, alphabet)) {
      alphabet->flags |= ALPH_FLAG_BUILTIN;
    }
    alph_release(reference);
  }

  return alphabet;
}

/*
 * alph_reader_had_warning
 * Return true if an input caused a recoverable warning.
 */
bool alph_reader_had_warning(ALPH_READER_T *reader) {
  return reader->had_warning;
}

/*
 * alph_reader_had_error
 * Return true if an input caused an unrecoverable error.
 */
bool alph_reader_had_error(ALPH_READER_T *reader) {
  return reader->had_error;
}

/*
 * alph_reader_has_message
 * Check if there are any pending warning or error mesages.
 */
bool alph_reader_has_message(ALPH_READER_T *reader) {
  return (linklst_size(reader->messages) > 0);
}

/*
 * alph_reader_next_message
 * Return the next pending warning or error message.
 * The caller is responsible for freeing memory.
 */
PARMSG_T* alph_reader_next_message(ALPH_READER_T *reader) {
  return (PARMSG_T*)linklst_pop(reader->messages);
}


#ifdef MAIN
// try loading an alphabet
int main(int argc, char **argv) {
  ALPH_T *alphabet, *dna, *rna, *protein;
  if (argc >= 2) {
    dna = alph_dna();
    rna = alph_rna();
    protein = alph_protein();
    alphabet = alph_load(argv[1], true);
    if (alphabet != NULL) {
      alph_print(alphabet, true, stdout);
      if (alph_equal(dna, alphabet)) {
        fprintf(stdout, "Alphabet equals builtin DNA!\n");
      }
      if (alph_equal(rna, alphabet)) {
        fprintf(stdout, "Alphabet equals builtin RNA!\n");
      }
      if (alph_equal(protein, alphabet)) {
        fprintf(stdout, "Alphabet equals builtin Protein!\n");
      }
      alph_release(alphabet);
    }
    alph_release(dna);
    alph_release(rna);
    alph_release(protein);
  }
  return 0;
}
#endif
