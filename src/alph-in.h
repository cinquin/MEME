#ifndef ALPH_IN_H
#define ALPH_IN_H

#include <stdbool.h>
#include <stdint.h>

#include "alphabet.h"
#include "parser-message.h"

typedef struct alph_reader ALPH_READER_T;

/*
 * alph_reader_create
 * Initialize the ALPH_READER_T
 */
ALPH_READER_T* alph_reader_create();

/*
 * alph_reader_destroy
 * Destroy the ALPH_READER_T
 */
void alph_reader_destroy(ALPH_READER_T *reader);

/*
 * alph_reader_header
 * Adds the header bypassing the text processing.
 */
void alph_reader_header(ALPH_READER_T *reader, int version, const char *name, int flags);

/*
 * alph_reader_core
 * Adds a core symbol bypassing the text processing.
 */
void alph_reader_core(ALPH_READER_T *reader, char symbol, const char* aliases, const char* name, int colour, char complement);

/*
 * alph_reader_ambig
 * Adds a ambiguous symbol bypassing the text processing.
 */
void alph_reader_ambig(ALPH_READER_T *reader, char symbol, const char* aliases, const char* name, int colour, const char *comprise);

/*
 * alph_reader_done
 * Finalize the input.
 */
void alph_reader_done(ALPH_READER_T *reader);

/*
 * alph_reader_line
 * Process a line of the file.
 */
void alph_reader_line(ALPH_READER_T *reader, const char *line);

/*
 * alph_reader_update
 * Update the reader with a chunk of text.
 */
void alph_reader_update(ALPH_READER_T *reader, const char *chunk, size_t size, bool end);

/*
 * alph_reader_had_warning
 * Return true if an input caused a recoverable warning.
 */
bool alph_reader_had_warning(ALPH_READER_T *reader);

/*
 * alph_reader_had_error
 * Return true if an input caused an unrecoverable error.
 */
bool alph_reader_had_error(ALPH_READER_T *reader);

/*
 * alph_reader_has_message
 * Check if there are any pending warning or error mesages.
 */
bool alph_reader_has_message(ALPH_READER_T *reader);

/*
 * alph_reader_next_message
 * Return the next pending warning or error message.
 * The caller is responsible for freeing memory.
 */
PARMSG_T* alph_reader_next_message(ALPH_READER_T *reader);

/*
 * alph_reader_alphabet
 * Get the loaded alphabet.
 * This will return NULL if the parser is not done or there are errors.
 */
ALPH_T* alph_reader_alphabet(ALPH_READER_T *reader);

#endif
