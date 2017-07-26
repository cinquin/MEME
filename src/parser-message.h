/**************************************************************************
 * FILE: parser-message.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 28 May 2014
 * PROJECT: shared
 * COPYRIGHT: UQ, 2014
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Enables storing messages and file coordinates from parsers
 **************************************************************************/

#ifndef PARSER_MESSAGE_H
#define PARSER_MESSAGE_H

#include <stdio.h>
#include <sys/types.h>

/*
 * An enumeration of error/warning severity
 */
typedef enum severity {
  SEVERITY_ERROR,
  SEVERITY_WARNING
} SEVERITY_EN;

/*
 * Structure for holding a parser message
 */
typedef struct parmsg {
  SEVERITY_EN severity;
  int64_t offset; // set to -1 when unused
  int64_t line; // set to -1 when unused
  int64_t column; // set to -1 when unused
  char* message; // short description; should not contain newline characters
} PARMSG_T;

/*
 * parmsg_set
 * Copies the parser message into a preallocated structure. 
 */
void parmsg_set(PARMSG_T* dest, SEVERITY_EN severity, int64_t offset, int64_t line, int64_t column, const char* message_format, ...);

/*
 * parmsg_free
 * Frees the parser message leaving the preallocated structure. 
 */
void parmsg_free(PARMSG_T* msg);

/*
 * parmsg_create
 * Creates the parser message.
 */
PARMSG_T* parmsg_create(SEVERITY_EN severity, int64_t offset, int64_t line, int64_t column, const char* message_format, ...);

/*
 * parmsg_destroy
 * Destroys the parser message.
 */
void parmsg_destroy(void* msg);

/*
 * parmsg_print
 * Print out the parser message.
 */
void parmsg_print(PARMSG_T* msg, FILE* out);

#endif
