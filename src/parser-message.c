
#include <stdarg.h>

#include "parser-message.h"
#include "utils.h"

/*
 * parmsg_vset
 * Copies the parser message into a preallocated structure. 
 */
static void parmsg_vset(PARMSG_T* dest, SEVERITY_EN severity, int64_t offset, int64_t line, int64_t column, const char* message_format, va_list ap) {
  va_list ap_copy;
  char temp[1];
  size_t length;
  dest->severity = severity;
  dest->offset = offset;
  dest->line = line;
  dest->column = column;
  // measure the message length
  va_copy(ap_copy, ap);
  length = vsnprintf(temp, 1, message_format, ap_copy);
  va_end(ap_copy);
  // write the message
  dest->message = mm_malloc(sizeof(char) * (length + 1));
  vsnprintf(dest->message, length + 1, message_format, ap);
}

/*
 * parmsg_set
 * Copies the parser message into a preallocated structure. 
 */
void parmsg_set(PARMSG_T* dest, SEVERITY_EN severity, int64_t offset, int64_t line, int64_t column, const char* message_format, ...) {
  va_list ap;
  va_start(ap, message_format);
  parmsg_vset(dest, severity, offset, line, column, message_format, ap);
  va_end(ap);
}

/*
 * parmsg_free
 * Frees the parser message leaving the preallocated structure. 
 */
void parmsg_free(PARMSG_T* msg) {
  free(msg->message);
  memset(msg, 0, sizeof(PARMSG_T));
}

/*
 * parmsg_create
 * Creates the parser message.
 */
PARMSG_T* parmsg_create(SEVERITY_EN severity, int64_t offset, int64_t line, int64_t column, const char* message_format, ...) {
  va_list ap;
  PARMSG_T *msg;
  msg = mm_malloc(sizeof(PARMSG_T));
  va_start(ap, message_format);
  parmsg_vset(msg, severity, offset, line, column, message_format, ap);
  va_end(ap);
  return msg;
}

/*
 * parmsg_destroy
 * Destroys the parser message.
 */
void parmsg_destroy(void* msg) {
  PARMSG_T* msg2;
  msg2 = (PARMSG_T*)msg;
  parmsg_free(msg2);
  free(msg2);
}

/*
 * parmsg_print
 * Print out the parser message.
 */
void parmsg_print(PARMSG_T* msg, FILE* out) {
  //TODO implement better
  if (msg->severity == SEVERITY_WARNING) {
    fputs("WARNING: ", out);
  } else if (msg->severity == SEVERITY_ERROR) {
    fputs("ERROR: ", out);
  }
  fprintf(out, "%s\n", msg->message);
}
