/********************************************************************
 * FILE: xlate-in.h
 * AUTHOR: James Johnson
 * CREATE DATE: 
 * DESCRIPTION: Create the structure for translating one alphabet into another.
 ********************************************************************/
#ifndef XLATE_IN_H
#define XLATE_IN_H

#include "alphabet.h"

typedef struct xlate_reader XLATE_READER_T;

XLATE_READER_T* xlate_reader_create(ALPH_T *src, ALPH_T *dest);

void xlate_reader_destroy(XLATE_READER_T* reader);

void xlate_reader_add(XLATE_READER_T* reader, const char *source, const char *dest);

void xlate_reader_done(XLATE_READER_T* reader);

/*
 * xlate_reader_had_warning
 * Return true if an input caused a recoverable warning.
 */
bool xlate_reader_had_warning(XLATE_READER_T *reader);

/*
 * xlate_reader_had_error
 * Return true if an input caused an unrecoverable error.
 */
bool xlate_reader_had_error(XLATE_READER_T *reader);

/*
 * xlate_reader_has_message
 * Check if there are any pending warning or error mesages.
 */
bool xlate_reader_has_message(XLATE_READER_T *reader);

/*
 * xlate_reader_next_message
 * Return the next pending warning or error message.
 * The caller is responsible for freeing memory.
 */
PARMSG_T* xlate_reader_next_message(XLATE_READER_T *reader);

/*
 * Adds in all missing translations using the wildcard when
 * there is no better option.
 */
XLATE_T* xlate_reader_translator(XLATE_READER_T* reader);

#endif
