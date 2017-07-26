
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/xmlstring.h>

#include "binary-search.h"
#include "sax-parser-utils.h"
#include "regex-utils.h"

/*****************************************************************************
 * Compare two pointers to strings. Used to binary search the attributes.
 ****************************************************************************/
int compare_pstrings(const void *v1, const void *v2) {
  char *str1, *str2;
  str1 = *((char**)v1);
  str2 = *((char**)v2);
  return strcmp(str1, str2);
}

/*****************************************************************************
 * Initilize the ATTRBUF_T
 ****************************************************************************/
void attrbuf_init(ATTRBUF_T *buf) {
  memset(buf, 0, sizeof(ATTRBUF_T));
  regcomp_or_die("XML entitity pattern", &(buf->pattern), "&(quot|amp|apos|lt|gt|#[0-9]+|#x[A-Za-z0-9]+);", REG_EXTENDED);
  buf->attrs_size = 4;
  buf->attrs = mm_malloc(sizeof(xmlChar*) * buf->attrs_size);
  buf->buffer_size = 128;
  buf->buffer = mm_malloc(sizeof(char) * buf->buffer_size);
}

/*****************************************************************************
 * Free the ATTRBUF_T
 ****************************************************************************/
void attrbuf_free(ATTRBUF_T *buf) {
  if (buf->buffer != NULL) free(buf->buffer);
  if (buf->attrs != NULL) free(buf->attrs);
  regfree(&(buf->pattern));
  memset(buf, 0, sizeof(ATTRBUF_T));
}

/*****************************************************************************
 * Translates any &X; into the represented symbol, stores the values temporally
 * in the ATTRBUF_T.
 ****************************************************************************/
const xmlChar** translate_attributes(ATTRBUF_T *buf, const xmlChar **attrs) {
  int i, length, entityLen, codepointLen;
  long codepoint;
  regmatch_t matches[2];
  char *buffer;
  // check that we actually have attributes
  if (attrs == NULL) return NULL;
  // measure the strings and count the attributes
  for (i = 0, length = 0; attrs[i] != NULL && attrs[i+1] != NULL; i += 2) {
    if (regexec_or_die("XML entity check", &(buf->pattern), (const char *)attrs[i+1], 2, matches, 0)) {
      length += strlen((const char *)attrs[i+1]) + 1;
    }
  }
  // if we don't need to modify anything then go ahead and use the original attributes
  if (length == 0) return attrs;
  // allocate enough space for a null terminated list of attributes
  if (buf->attrs_size < i) {
    buf->attrs_size = i + 1;
    buf->attrs = mm_realloc(buf->attrs, sizeof(xmlChar*) * buf->attrs_size);
  }
  // allocate enough space for the strings
  if (buf->buffer_size < length) {
    buf->buffer_size = length;
    buf->buffer = mm_realloc(buf->buffer, sizeof(char) * buf->buffer_size);
  }
  // rewrite the attributes that have character entities
  buffer = buf->buffer;
  for (i = 0; attrs[i] != NULL && attrs[i+1] != NULL; i += 2) {
    char *value;
    // set the key attribute
    buf->attrs[i] = attrs[i];
    // set the value attribute (but we will change this if it contains character entities)
    buf->attrs[i+1] = attrs[i+1];
    // determine if we need to modify the value attribute
    value = (char*)attrs[i+1];
    if (regexec_or_die("XML entity check", &(buf->pattern), value, 2, matches, 0)) {
      // point the value attribute at the buffer
      buf->attrs[i+1] = (xmlChar*)buffer;
      // copy the string and replace character entities
      do {
        // copy the bit before the &
        strncpy(buffer, value, matches[0].rm_so);
        buffer += matches[0].rm_so;
        // now translate the character entity
        if (value[matches[1].rm_so] == '#') { // number
          entityLen = matches[1].rm_eo - matches[1].rm_so;
          if (value[matches[1].rm_so + 1] == 'x') { // base 16
            strncpy(buffer, value+(matches[1].rm_so + 2), entityLen - 2);
            buffer[entityLen - 2] = '\0';
            codepoint = strtol(buffer, NULL, 16);
          } else { // base 10
            strncpy(buffer, value+(matches[1].rm_so + 1), entityLen - 1);
            buffer[entityLen - 1] = '\0';
            codepoint = strtol(buffer, NULL, 10);
          }
          unicode_to_string(codepoint, buffer, &codepointLen);
          buffer += codepointLen;
        } else if (regex_cmp(matches+1, value, "lt") == 0) {
          buffer[0] = '<';
          buffer++;
        } else if (regex_cmp(matches+1, value, "gt") == 0) {
          buffer[0] = '>';
          buffer++;
        } else if (regex_cmp(matches+1, value, "amp") == 0) {
          buffer[0] = '&';
          buffer++;
        } else if (regex_cmp(matches+1, value, "quot") == 0) {
          buffer[0] = '"';
          buffer++;
        } else if (regex_cmp(matches+1, value, "apos") == 0) {
          buffer[0] = '\'';
          buffer++;
        } else {
          die("This should be unreachable because the pattern should not match any other alternatives");
        }
        // move the value pointer to after the character entity
        value += matches[0].rm_eo;
      } while (regexec_or_die("XML entity check", &(buf->pattern), value, 2, matches, 0));
      // copy the remaining characters
      for (; value[0] != '\0'; value++, buffer++) buffer[0] = value[0];
      buffer[0] = '\0';
      buffer++;
    }
  }
  // NULL terminate the buffer
  buf->attrs[i] = NULL;
  return buf->attrs;
}

/*****************************************************************************
 * Compare a program version and a constant version
 ****************************************************************************/
int prog_version_cmp(const struct prog_version *version, const int major, const int minor, const int patch) {
  if (version->major > major) {
    return 1;
  } else if (version->major == major) {
    if (version->minor > minor) {
      return 1;
    } else if (version->minor == minor) {
      if (version->patch > patch) {
        return 1;
      } else if (version->patch == patch) {
        return 0;
      }
    }
  }
  return -1;
}

/*****************************************************************************
 * Load a string attribute when used with parse attributes. 
 * Expects that the data attribute should be a pointer to a string.
 * Note that users will need to duplicate the string should they wish to 
 * keep it.
 ****************************************************************************/
int ld_str(char *value, void *data) {
  *((char**)data) = value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a character attribute.
 * Expects that the data attribute should be a pointer to a character.
 ****************************************************************************/
int ld_char(char *value, void *data) {
  int i;
  char c;
  for (i = 0; value[i] != '\0'; ++i) {
    if (!isspace(value[i])) break;
  }
  //check that we found something
  if (value[i] == '\0') {
    return -1;
  }
  c = value[i];
  //check that it was the only thing to be found
  for (++i; value[i] != '\0'; ++i) {
    if (!isspace(value[i])) return -2;
  }
  //set the value
  *((char*)data) = c;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a hexadecimal integer attribute.
 * Expects that the data attribute should be a pointer to an integer.
 ****************************************************************************/
int ld_hex(char *value, void *data) {
  long parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtol(value, &end_ptr, 16);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  if (parsed_value > INT_MAX || parsed_value < INT_MIN) return ERANGE;
  *((int*)data) = (int)parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a integer attribute.
 * Expects that the data attribute should be a pointer to an integer.
 ****************************************************************************/
int ld_int(char *value, void *data) {
  long parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtol(value, &end_ptr, 10);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  if (parsed_value > INT_MAX || parsed_value < INT_MIN) return ERANGE;
  *((int*)data) = (int)parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a long integer attribute.
 * Expects that the data attribute should be a pointer to a long integer.
 ****************************************************************************/
int ld_long(char *value, void *data) {
  long parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtol(value, &end_ptr, 10);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  *((long*)data) = parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a double precision floating point 
 * number attribute.
 * Expects that the data attribute should be a pointer to a double.
 ****************************************************************************/
int ld_double(char *value, void *data) {
  double parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtod(value, &end_ptr);
  if (end_ptr == value) return -1;
  // allow out of range values, mainly because evalues can be
  // stupendously small or large
  if (errno && errno != ERANGE) return errno;
  *((double*)data) = parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a double that is in the range 0 to 1
 * inclusive and hence could be used to represent a probability.
 * Expects that the data attribute should be a pointer to a double.
 ****************************************************************************/
int ld_pvalue(char *value, void *data) {
  double parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtod(value, &end_ptr);
  if (end_ptr == value) return -1;
  // allow out of range values because pvalues can be really small
  if (errno && errno != ERANGE) return errno;
  if (parsed_value < 0 || parsed_value > 1) return -1;
  *((double*)data) = parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load an E-value which may be too small to
 * represent with a double.
 * Expects that the data attribute should be a pointer to a double.
 ****************************************************************************/
int ld_log10_ev(char *value, void *data) {
  double parsed_value;
  errno = 0;
  parsed_value = log10_evalue_from_string(value);
  if (errno) return errno;
  *((double*)data) = parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load an p-value which may be too small to
 * represent with a double.
 * Expects that the data attribute should be a pointer to a double.
 ****************************************************************************/
int ld_log10_pv(char *value, void *data) {
  double parsed_value;
  errno = 0;
  parsed_value = log10_evalue_from_string(value);
  if (errno) return errno;
  // p-value, so must be smaller than log(1)
  if (parsed_value > 0) return 1;
  *((double*)data) = parsed_value;
  return 0;
}

/*****************************************************************************
 * Use with parse attributes to load a multi-valued attribute.
 * Expects that the data attribute is a pointer to a structure of type 
 * MULTI_T which contains a list of possible strings for the multi-valued 
 * attribute and a list of integer output codes to represent them as well
 * as the pointer to the integer that is to store the output. Note
 * that the list of strings must be alphabetically ordered.
 ****************************************************************************/
int ld_multi(char *value, void *data) {
  int found;
  MULTI_T *attr;
  attr = (MULTI_T*)data;
  found = binary_search(&value, attr->options, attr->count, sizeof(char*), compare_pstrings);
  if (found >= 0) {
    *(attr->target) = attr->outputs[found];
    return 0;
  } else {
    return -1;
  }
}

/*****************************************************************************
 * A function to load the program version from a string and put it into the
 * prog_version structure passed as data.
 ****************************************************************************/
int ld_version(char *value, void *data) {
  struct prog_version *ver;
  char *start, *end;
  ver = (struct prog_version *)data;
  ver->major = 0;
  ver->minor = 0;
  ver->patch = 0;
  if (*value == '\0') return -1;
  start = value;
  ver->major = strtol(start, &end, 10);
  if (start == end) return -1;
  if (*end == '\0') return 0;
  if (*end != '.') return -1;
  start = end+1;
  ver->minor = strtol(start, &end, 10);
  if (start == end) return -1;
  if (*end == '\0') return 0;
  if (*end != '.') return -1;
  start = end+1;
  ver->patch = strtol(start, &end, 10);
  if (start == end) return -1;
  if (*end != '\0') return -1;
  return 0;
}

/*****************************************************************************
 * Parses the list of attributes given to the libxml startElement callback.
 * To do the parsing it takes a list of alphabetically ordered attribute names
 * and a list of associated parsers, in-out data to be passed to the parsers,
 * booleans to indicate if the attribute is required and booleans to indicate
 * if the attribute was actually read.
 ****************************************************************************/
void parse_attributes(void (*error)(void *, int, const char*, const char*, const char*), 
    void *state, char *tag, const xmlChar **attrs, 
    int count, char **names, int (**parsers)(char*, void*), void **parser_data, 
    BOOLEAN_T *required, BOOLEAN_T *done) {
  int i, found;
  BOOLEAN_T done_data[count];
  if (done == NULL) done = done_data;
  for (i = 0; i < count; ++i) done[i] = FALSE;
  for (i = 0; attrs[i] != NULL; i += 2) {
    found = binary_search(attrs+i, names, count, sizeof(char*), compare_pstrings);
    if (found >= 0) {
      int (*parser)(char*, void*);
      void* data;
      if (done[found]) {
        error(state, PARSE_ATTR_DUPLICATE, tag, names[found], NULL);
        continue;
      }
      done[found] = TRUE;
      parser = parsers[found];
      data = parser_data[found];
      if (parser((char*)attrs[i+1], data)) {
        error(state, PARSE_ATTR_BAD_VALUE, tag, names[found], (const char*)attrs[i+1]);
        continue;
      }
    }
  }
  for (i = 0; i < count; ++i) {
    if (required[i] && !done[i]) {
      error(state, PARSE_ATTR_MISSING, tag, names[i], NULL);
    }
  }
}

/*****************************************************************************
 * Character filter that accepts everything
 ****************************************************************************/
int accept_all_chars(char ch) {
  return 1;
}

/*****************************************************************************
 * Character filter that accepts everything that isn't a space
 ****************************************************************************/
int accept_all_chars_but_isspace(char ch) {
  return !isspace(ch);
}

/*****************************************************************************
 * Stores characters which pass the acceptance test into the buffer 
 * automatically expanding the buffer if required.
 ****************************************************************************/
void store_xml_characters(CHARBUF_T *buf, const char *ch, int len) {
  int (*accept) (char);
  int i, start, end, accepted;

  accept = buf->accept;
  if (accept) {
    i = 0;
    while (i < len) {
      //find how many to skip
      for (; i < len; ++i) if (accept(ch[i])) break;
      if (i >= len) return; //nothing left
      start = i;
      //find how many to accept
      for (++i; i < len; ++i) if (!(accept(ch[i]))) break;
      end = i;
      accepted = end - start;
      //increase buffer size if needed
      if ((buf->pos + accepted) >= buf->size) {
        buf->size = buf->pos + accepted + 1;
        buf->buffer = mm_realloc(buf->buffer, buf->size * sizeof(char));
      }
      //copy over the characters to the buffer
      for (i = start; i < end; ++i) {
        buf->buffer[(buf->pos)++] = ch[i];
      }
      buf->buffer[buf->pos] = '\0';
    }
  }
}
