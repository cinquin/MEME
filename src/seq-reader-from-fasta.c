/******************************************************************************
 * FILE: seq-reader-from-fasta.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-17
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the concrete implementation for a 
 * data-block-reader UDT for reading sequence segments from a FASTA file.
 * CONSIDER: This and the PSP_DATA_BLOCK_READER_T data structure are almost 
 * identical. Can they be refactored to share more of the implementation?
 *****************************************************************************/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <regex.h>
#include "alphabet.h"
#include "seq-reader-from-fasta.h"

const size_t BUFFER_SIZE = 5000;

typedef struct seq_reader_from_fasta {
  BOOLEAN_T at_end_of_file;
  BOOLEAN_T at_start_of_line;
  BOOLEAN_T at_end_of_seq;
  BOOLEAN_T parse_genomic_coord;
  int current_position;
  char *seq_buffer;
  size_t seq_buffer_index;
  size_t seq_buffer_last;
  char* filename;
  size_t filename_len; // Includes trailing '\0'
  char* seq_header;
  size_t seq_header_len; // Includes trailing '\0'
  size_t header_buffer_len;
  char* seq_name;
  size_t seq_name_len;
  FILE *fasta_file;
  long prev_block_position; 
  ALPH_T* alphabet;
} SEQ_READER_FROM_FASTA_T;

// Forward declarations

void free_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T seq_reader_from_fasta_is_eof(DATA_BLOCK_READER_T *reader);
BOOLEAN_T reset_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T close_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T read_seq_header_from_seq_reader_from_fasta(
  SEQ_READER_FROM_FASTA_T *fasta_reader
);
BOOLEAN_T get_seq_name_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  char **name // OUT
);
BOOLEAN_T go_to_next_sequence_in_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader
);
BOOLEAN_T get_next_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  DATA_BLOCK_T *data_block
);
BOOLEAN_T unget_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader
);

/******************************************************************************
 * This function creates an instance of a data block reader UDT for reading
 * sequence segments from a FASTA file.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_seq_reader_from_fasta(
  BOOLEAN_T parse_genomic_coord, 
  ALPH_T *alph, 
  const char *filename
) {
  SEQ_READER_FROM_FASTA_T *fasta_reader = mm_malloc(sizeof(SEQ_READER_FROM_FASTA_T) * 1);
  fasta_reader->at_end_of_file = FALSE;
  fasta_reader->at_start_of_line = TRUE;
  fasta_reader->at_end_of_seq = FALSE;
  fasta_reader->parse_genomic_coord = parse_genomic_coord;
  int filename_len = strlen(filename) + 1;
  fasta_reader->filename = mm_malloc(sizeof(char)* filename_len);
  fasta_reader->filename_len = filename_len;
  strncpy(fasta_reader->filename, filename, filename_len);
  fasta_reader->prev_block_position = 0;
  fasta_reader->current_position = 0;
  fasta_reader->seq_header = NULL;
  fasta_reader->seq_header_len = 0;
  fasta_reader->seq_name = NULL;
  fasta_reader->seq_name_len = 0;
  fasta_reader->alphabet = alph_hold(alph);
  fasta_reader->seq_buffer = mm_calloc(sizeof(char), BUFFER_SIZE);
  fasta_reader->seq_buffer_index = 0;
  fasta_reader->seq_buffer_last = 0;
  if (
    open_file(
      filename, 
      "r", 
      TRUE, 
      "FASTA", 
      "sequences", 
      &(fasta_reader->fasta_file
    )) == FALSE) {
    die("Couldn't open the file %s.\n", filename);
  }

  DATA_BLOCK_READER_T *reader = new_data_block_reader(
    (void *) fasta_reader,
    free_seq_reader_from_fasta,
    close_seq_reader_from_fasta,
    reset_seq_reader_from_fasta,
    seq_reader_from_fasta_is_eof,
    get_next_data_block_from_seq_reader_from_fasta,
    unget_data_block_from_seq_reader_from_fasta,
    go_to_next_sequence_in_seq_reader_from_fasta,
    get_seq_name_from_seq_reader_from_fasta
  );
  return reader;
}

/******************************************************************************
 * This function frees an instance of the sequence FASTA reader UDT.
 *****************************************************************************/
void free_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  myfree(fasta_reader->seq_buffer);
  fasta_reader->seq_buffer_index = 0;
  fasta_reader->seq_buffer_last = 0;
  myfree(fasta_reader->filename);
  fasta_reader->filename_len = 0;
  myfree(fasta_reader->seq_header);
  fasta_reader->seq_header_len = 0;
  fasta_reader->seq_header_len = 0;
  myfree(fasta_reader->seq_name);
  fasta_reader->seq_name_len = 0;
  fasta_reader->seq_name_len = 0;
  alph_release(fasta_reader->alphabet);
  myfree(fasta_reader);
}

/******************************************************************************
 * This function closes a sequence FASTA reader UDT.
 *****************************************************************************/
BOOLEAN_T close_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  fasta_reader->current_position = 0;
  if (fasta_reader->fasta_file != NULL) {
    if (fclose(fasta_reader->fasta_file) == EOF) {
      die(
        "Error closing file: %s.\nError message: %s\n", 
        fasta_reader->filename, 
        strerror(errno)
      );
    }
    else {
      result = TRUE;
    }
  }
  return result;
}

/******************************************************************************
 * This function resets a sequence FASTA reader UDT.
 *****************************************************************************/
BOOLEAN_T reset_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  if (fasta_reader->fasta_file == stdin) {
    // TODO it would be nice if this could be made to work using temporary files
    die("Unable to rewind when reading sequence from standard input\n");
  }
  else {
    rewind(fasta_reader->fasta_file);
  }
  memset(fasta_reader->seq_buffer, 0, BUFFER_SIZE);
  fasta_reader->seq_buffer_index = 0;
  fasta_reader->seq_buffer_last = 0;
  fasta_reader->current_position = -1;
  fasta_reader->at_end_of_file = FALSE;
  fasta_reader->at_start_of_line = TRUE;
  fasta_reader->at_end_of_seq = FALSE;
  return TRUE;
}

/******************************************************************************
 * This function reports on whether a prior reader has reached EOF
 * Returns TRUE if the reader is at EOF
 *****************************************************************************/
BOOLEAN_T seq_reader_from_fasta_is_eof(DATA_BLOCK_READER_T *reader) {
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  return feof(fasta_reader->fasta_file) ? TRUE : FALSE;
}

/******************************************************************************
 * This function attempts to parse genomic coordinates from the
 * current sequence header. If successful it will set the sequence
 * name to the chromosome name, and set the starting sequence position.
 *
 * Returns TRUE if it was able to find genomic coordinates, FALSE otherwise.
 *****************************************************************************/
static BOOLEAN_T parse_genomic_coordinates(
  SEQ_READER_FROM_FASTA_T *fasta_reader
) {
  // Copy chromosome name and position to reader
  int end_pos;
  return (
    parse_genomic_coordinates_helper(
      fasta_reader->seq_header,
      &(fasta_reader->seq_name),
      &(fasta_reader->seq_name_len),
      &(fasta_reader->current_position),
      &end_pos
    )
  );
}

/******************************************************************************
 * This function attempts to parse genomic coordinates from the
 * current sequence header. If successful it will return the 
 * chromosome name, name length and starting position.
 *
 * There are two supported formats.
 * The first is the UCSC/fastaFromBed style: name:start-stop[(+|-)][_id]
 * e.g., ">chr1:1000-1010(-)_xyz"
 * The second is the Galaxy "Fetch sequence" style: genome_name_start_stop_strand
 * where strand is +|-.
 * e.g., ">mm9_chr18_75759530_7575972_-"
 *
 * Returns TRUE if it was able to find genomic coordinates, FALSE otherwise.
 *****************************************************************************/
BOOLEAN_T parse_genomic_coordinates_helper(
  char*  header,	  // sequence name 
  char** chr_name_ptr,    // chromosome name ((chr[^:]))
  size_t * chr_name_len_ptr,// number of characters in chromosome name
  int *  start_ptr,	  // start position of sequence (chr:(\d+)-)
  int *  end_ptr	  // end position of sequence (chr:\d+-(\d+))
) {

  #define NUM_SUBMATCHES 6
  #define ERROR_MESSAGE_SIZE 100
  #define BUFFER_SIZE 512

  static BOOLEAN_T first_time = TRUE;
  static regex_t ucsc_header_regex;
  static regex_t galaxy_header_regex;
  static regmatch_t matches[NUM_SUBMATCHES];

  char error_message[ERROR_MESSAGE_SIZE];
  int status = 0;

  if (first_time == TRUE) {

    // Initialize regular express for extracting chromsome coordinates;

    status = regcomp(
      &ucsc_header_regex, 
      "([^[:space:]:]+):([[:digit:]]+)-([[:digit:]]+)(\\([+-]\\))?(_[^[:space:]]+)?", 
      REG_EXTENDED | REG_ICASE | REG_NEWLINE
    );

    if (status != 0) {
      regerror(status, &ucsc_header_regex, error_message, ERROR_MESSAGE_SIZE);
      die(
        "Error while intitializing regular expression\n"
        "for parsing UCSC style genome coordinates from FASTA header: %s\n", 
        error_message
      );
    }

    status = regcomp(
      &galaxy_header_regex, 
      "([^[:space:]_]+_[^[:space:]_]+)_([[:digit:]]+)_([[:digit:]]+)(_[+-])?", 
      REG_EXTENDED | REG_ICASE | REG_NEWLINE
    );

    if (status != 0) {
      regerror(status, &galaxy_header_regex, error_message, ERROR_MESSAGE_SIZE);
      die(
        "Error while intitializing regular expression\n"
        "for parsing Galaxy style genome coordinates from FASTA header: %s\n", 
        error_message
      );
    }

    first_time = FALSE;
  }

  BOOLEAN_T found_coordinates = FALSE;

  // Try UCSC style first
  status = regexec(&ucsc_header_regex, header, NUM_SUBMATCHES, matches, 0);
  if(status && status != REG_NOMATCH ){
    regerror(status, &ucsc_header_regex, error_message, 100);
    die("Error trying to parse UCSC style genome coordinates from sequence header: %s\n"
        "error message is %s\n",
        header, 
        error_message
    );
  }

  if (status == REG_NOMATCH) {
    // UCSC didn't work, try GALAXY style
    status = regexec(&galaxy_header_regex, header, NUM_SUBMATCHES, matches, 0);
    if(status && status != REG_NOMATCH ){
      regerror(status, &ucsc_header_regex, error_message, 100);
      die("Error trying to parse GALAXY style genome coordinates from sequence header: %s\n"
          "error message is %s\n",
          header, 
          error_message
      );
    }
  }

  if(!status) {

    // The sequence header contains genomic coordinates
    found_coordinates = TRUE;

    // Get chromosome name (required)
    char buffer[BUFFER_SIZE];
    int chr_name_len = matches[1].rm_eo - matches[1].rm_so;
    strncpy(buffer, header + matches[1].rm_so, chr_name_len);
    buffer[chr_name_len] = 0;
    char *chr_name = strdup(buffer);
    if (chr_name == NULL) {
      die("Unable to allocate memory while parsing sequence header.\n");
    }
    *chr_name_ptr = chr_name;
    *chr_name_len_ptr = chr_name_len;

    // Get chromosome start position (required)
    size_t pos_len = matches[2].rm_eo - matches[2].rm_so;
    strncpy(buffer, header + matches[2].rm_so, pos_len);
    buffer[pos_len] = 0;
    *start_ptr = atol(buffer) - 1;

    // Get chromosome end position (required)
    pos_len = matches[3].rm_eo - matches[3].rm_so;
    strncpy(buffer, header + matches[3].rm_so, pos_len);
    buffer[pos_len] = 0;
    *end_ptr = atol(buffer) - 1;

    // Get the strand (optional)
    if (matches[4].rm_so >= 0) {
      pos_len = matches[4].rm_eo - matches[4].rm_so;
      char strand = *(header + matches[4].rm_so + 1);
    }

    // Get the id tag (optional)
    if (matches[5].rm_so >= 0) {
      pos_len = matches[5].rm_eo - matches[5].rm_so;
      strncpy(buffer, header + matches[5].rm_so + 1, pos_len);
      buffer[pos_len] = 0;
    }

  }

  return found_coordinates;

  #undef NUM_SUBMATCHES
  #undef ERROR_MESSAGE_SIZE
  #undef BUFFER_SIZE
}

/******************************************************************************
 * This function attempts to parse the sequence name from the
 * current sequence header. The sequence name is the string up to the first
 * white space in the sequence header.
 *
 * Returns TRUE if it was able to genomic coordinates, FALSE otherwise.
 *****************************************************************************/
static BOOLEAN_T parse_seq_name(
  SEQ_READER_FROM_FASTA_T *fasta_reader
) {

  #define NUM_SUBMATCHES 2
  #define ERROR_MESSAGE_SIZE 100
  #define BUFFER_SIZE 512

  static BOOLEAN_T first_time = TRUE;
  static regex_t ucsc_header_regex;
  static regmatch_t matches[NUM_SUBMATCHES];

  char error_message[ERROR_MESSAGE_SIZE];
  int status = 0;

  if (first_time == TRUE) {
    // Initialize regular express for extracting chromsome coordinates;
    // Expected format is based on fastaFromBed name:start-stop[(+|-)][_id]
    // e.g., ">chr1:1000-1010(-)_xyz"
    char *regexp_str = fasta_reader->parse_genomic_coord ?
      "([^[:space:]:]+)[[:space:]]*" : "([^[:space:]]+)[[:space:]]*";
    status = regcomp(
      &ucsc_header_regex,
      regexp_str,
      REG_EXTENDED | REG_ICASE | REG_NEWLINE
    );

    if (status != 0) {
      regerror(status, &ucsc_header_regex, error_message, ERROR_MESSAGE_SIZE);
      die(
        "Error while intitializing regular expression\n"
        "for parsing sequence name from FASTA header: %s\n", 
        error_message
      );
    }

    first_time = FALSE;
  }

  BOOLEAN_T found_name = FALSE;
  char* header = fasta_reader->seq_header;
  status = regexec(
    &ucsc_header_regex, 
    header,
    NUM_SUBMATCHES, 
    matches, 
    0
  );
  if(!status) {
    // The sequence header contains genomic coordinates
    found_name = TRUE;

    // Copy sequence name to reader
    char buffer[BUFFER_SIZE];
    int name_len = matches[1].rm_eo - matches[1].rm_so;
    strncpy(buffer, header + matches[1].rm_so, name_len);
    buffer[name_len] = 0;
    if (fasta_reader->seq_name) {
      myfree(fasta_reader->seq_name);
      fasta_reader->seq_name = NULL;
    }
    fasta_reader->seq_name = strdup(buffer);
    if (fasta_reader->seq_name == NULL) {
      die("Unable to allocate memory while parsing sequence header.\n");
    }
    fasta_reader->seq_name_len = name_len;
  }
  else if(status != REG_NOMATCH ){
    regerror(status, &ucsc_header_regex, error_message, 100);
    die("Error trying to parse name from sequence header: %s\n"
        "error message is %s\n",
        header, 
        error_message
    );
  }
  else {
    found_name = FALSE;
  }

  return found_name;

  #undef NUM_SUBMATCHES
  #undef ERROR_MESSAGE_SIZE
  #undef BUFFER_SIZE
}


/******************************************************************************
 * This function reads the entire sequence header at the start of a new sequence.
 * The current position is assumed to be start of a new sequence.
 * Read from the current position to the end of the current line.
 *
 * Returns TRUE if it was able to read the sequence text, FALSE if 
 * EOF reached before the terminal newline was found. Dies if other errors
 * are encountered.
 *****************************************************************************/
BOOLEAN_T read_seq_header_from_seq_reader_from_fasta(
  SEQ_READER_FROM_FASTA_T *fasta_reader
) {

  int result = FALSE;

  // Initial allocation of sequence buffer
  const size_t initial_buffer_len = 100;
  if (fasta_reader->seq_header == NULL) {
     fasta_reader->seq_header = mm_malloc(sizeof(char) * initial_buffer_len);
     fasta_reader->header_buffer_len = initial_buffer_len;
  }

  // Look for EOL
  int c = 0;
  int seq_index = 0;
  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {

    if (seq_index >= fasta_reader->header_buffer_len) {
      // Need to grow buffer
      fasta_reader->seq_header
        = mm_realloc(
            fasta_reader->seq_header, 
            2 * fasta_reader->header_buffer_len
          );
      fasta_reader->header_buffer_len = 2 * fasta_reader->header_buffer_len;
    }

    if (c == '\n') {
      // Found EOL
      fasta_reader->seq_header[seq_index] = '\0';
      fasta_reader->seq_header_len = seq_index + 1;
      fasta_reader->at_start_of_line = TRUE;
      result = TRUE;
      break;
    }
    else {
      // Keep looking for EOL
      fasta_reader->seq_header[seq_index] = c;
      ++seq_index;
    }

  }

  // At this point c is EOL or EOF
  if (c == EOF) {
    if (ferror(fasta_reader->fasta_file)) {
      // EOF could actually indicate an error.
      die(
        "Error reading file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }
    else if (feof(fasta_reader->fasta_file)) {
        // Reached EOF before reaching EOL for the sequence.
        fasta_reader->seq_header[0] = '\0';
        fasta_reader->seq_header_len = 0;
    }
  }

  return result;

}

/******************************************************************************
 * This function gets the name of the current sequence from the data block
 * reader. The name of the sequence is passed using the name parameter.
 * The caller is responsible for freeing the memory for the sequence name.
 *
 * Returns TRUE if successful, FALSE if there is no current sequence, as 
 * at the start of the file.
 *****************************************************************************/
BOOLEAN_T get_seq_name_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  char **name // OUT
) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  if (fasta_reader->seq_name == NULL || fasta_reader->seq_name_len <= 0) {
    result = FALSE;
  }
  else {
    *name = strdup(fasta_reader->seq_name);
    result = TRUE;
  }

  return result;
}

/******************************************************************************
 * Read from the current position in the file to the first symbol after the
 * start of the next sequence. Set the value of the current sequence.
 *
 * Returns TRUE if it was able to advance to the next sequence, FALSE if 
 * EOF reached before the next sequence was found. Dies if other errors
 * encountered.
 *****************************************************************************/
BOOLEAN_T go_to_next_sequence_in_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  fasta_reader->current_position = 0;
  int c = 0;
  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {
    if (fasta_reader->at_start_of_line == TRUE && c == '>') {
      break;
    }
    else if (c == '\n') {
      fasta_reader->at_start_of_line = TRUE;
    }
    else {
      fasta_reader->at_start_of_line = FALSE;
    }
  }
  // At this point c is '>' or EOF
  if (c == '>') {
    BOOLEAN_T found_genomic_coordinates = FALSE;
    fasta_reader->at_end_of_seq = FALSE;
    fasta_reader->seq_buffer_index = 0;
    fasta_reader->seq_buffer_last = 0;
    memset(fasta_reader->seq_buffer, 0, BUFFER_SIZE);
    result = read_seq_header_from_seq_reader_from_fasta(fasta_reader);
    if (result == TRUE && fasta_reader->parse_genomic_coord == TRUE) {
      // Look for genomic coordinates in header
      found_genomic_coordinates = parse_genomic_coordinates(fasta_reader);
    }
    if (found_genomic_coordinates == FALSE) {
      //  Look for whitespace in header
      //  The sequence name is the string before the white space.
      BOOLEAN_T found_name = FALSE;
      found_name = parse_seq_name(fasta_reader);
      if (found_name == FALSE) {
        die(
            "Unable to find sequence name in header %s.\n",
            fasta_reader->seq_header
        );
      }
    }
  }
  else {
    if (ferror(fasta_reader->fasta_file)) {
      die(
        "Error reading file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }
    else if (feof(fasta_reader->fasta_file)) {
        // Reached EOF before reaching the start of the sequence
        result = FALSE;
    }
  }
  return result;
}

size_t get_current_pos_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader
) {

  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  return fasta_reader->current_position;
}


/****************************************************************************
 * Read up to max_chars letters of one sequence from a DATA_BLOCK_T readder
 * and copy them in to the raw sequence in the SEQ_T object starting at the
 * given buffer offset. 
 ****************************************************************************/
void read_one_fasta_segment_from_reader(
   DATA_BLOCK_READER_T *fasta_reader,
   size_t max_size,
   size_t offset,
   SEQ_T *sequence
) {


  assert(sequence != NULL);
  assert(offset < max_size);

  // Get the raw sequence buffer from the SEQ_T
  char *raw_sequence = get_raw_sequence(sequence);
  if (raw_sequence == NULL) {
    // Allocate space for raw sequence if not done yet.
    raw_sequence = mm_malloc(sizeof(char) * max_size + 1);
    raw_sequence[0] = 0;
  }

  // Read a block of sequence charaters into the
  // raw sequence buffer for the SEQ_T, starting at offset.
  BOOLEAN_T is_complete = read_raw_sequence_from_reader(
    fasta_reader,
    max_size - offset,
    raw_sequence + offset
  );
  set_raw_sequence(raw_sequence, is_complete, sequence);
}

/****************************************************************************
 * Read raw sequence until a new sequence is encountered or too many letters
 * are read.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
BOOLEAN_T read_raw_sequence_from_reader(
   DATA_BLOCK_READER_T *reader, // Sequence source
   unsigned int max_chars, // Maximum chars in raw_sequence.
   char* raw_sequence // Pre-allocated sequence buffer.
) {

  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  // Read sequence into temp. buffer from the sequence file.
  char buffer[max_chars];
  long start_file_pos = ftell(fasta_reader->fasta_file);
  size_t seq_index = 0;
  size_t total_read = 0;
  while (seq_index < max_chars) {

    size_t num_char_read = fread(
      buffer,
      sizeof(char), 
      max_chars - seq_index,
      fasta_reader->fasta_file
    );
    fasta_reader->current_position += num_char_read;
    total_read += num_char_read;

    if (feof(fasta_reader->fasta_file)) {
       fasta_reader->at_end_of_file = TRUE;
    }
    else if (num_char_read < (max_chars - seq_index)) {
      die(
        "Error while reading sequence from file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }

    size_t i;
    for(i = 0; i < num_char_read; ++i) {
      char c = buffer[i];
      assert(c != 0);
      if (isspace(c)) {
        // Skip over white space
        fasta_reader->at_start_of_line = (c == '\n');
      }
      else if (c == '>' && fasta_reader->at_start_of_line == TRUE) {
        // We found the start of a new sequence while trying
        // to fill the buffer. Leave the buffer incomplete.
        // and wind back the file
        fseek(fasta_reader->fasta_file, start_file_pos + i - 1, SEEK_SET);
        fasta_reader->current_position = start_file_pos + i - 1;
        fasta_reader->at_end_of_seq = TRUE;
        fasta_reader->at_start_of_line = FALSE;
        fasta_reader->at_end_of_file = FALSE;
        break;
      }
      else {
        fasta_reader->at_start_of_line = FALSE;
        // Check that character is legal in alphabet. 
        // If not, replace with wild card character.
        if (alph_is_known(fasta_reader->alphabet, c)) {
          raw_sequence[seq_index] = c;
        }
        else {
          raw_sequence[seq_index] = alph_wildcard(fasta_reader->alphabet);
          fprintf(
            stderr, 
            "Warning: %c is not a valid character in %s alphabet.\n"
            "         Converting %c to %c.\n",
            c,
            alph_name(fasta_reader->alphabet),
            c,
            raw_sequence[i]
          );
        }
        ++seq_index;
      }
    }
    if (fasta_reader->at_end_of_seq | fasta_reader->at_end_of_file) {
      break;
    }
  }

  raw_sequence[seq_index] = '\0';
  return(fasta_reader->at_end_of_seq | fasta_reader->at_end_of_file);
}

/******************************************************************************
 * Populates the data block for the with the next block of sequence. 
 *
 * During the first call for the sequence it fills in a buffer from a file,
 * The sequence pointer in the data block is set to point at the start of the buffer.
 * On successive calls, the sequence pointer in the block is shifted down one position
 * in the buffer. When the end of the buffer is reached, it is filled again from the file.
 * 
 * Returns TRUE if it was able to completely fill the block, FALSE if 
 * the next sequence or EOF was reached before the block was filled.
 * Dies if other errors encountered.
 *****************************************************************************/
BOOLEAN_T get_next_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  DATA_BLOCK_T *data_block
) {

  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  int w = get_block_size_from_data_block(data_block);

  if ((fasta_reader->seq_buffer_last - fasta_reader->seq_buffer_index) > w) {
    // Enough data in buffer to provide next block
    // Update data block with updated pointer into seq buffer
    // and we're done.
    set_sequence_in_data_block(
      data_block,
      fasta_reader->seq_buffer + fasta_reader->seq_buffer_index
    );
    set_num_read_into_data_block(data_block, w);
    ++fasta_reader->seq_buffer_index;
    ++fasta_reader->current_position;
    set_start_pos_for_data_block(data_block, fasta_reader->current_position);
    return TRUE;
  }

  int num_char_to_read;
  if (fasta_reader->seq_buffer_last == 0) {
    // Seq buffer is completely empty
    fasta_reader->seq_buffer_index = 0;
    num_char_to_read = BUFFER_SIZE;
  }
  else if (!fasta_reader->at_end_of_seq && !fasta_reader->at_end_of_file) {
    // At the end of the seq buffer. Copy the last data block's worth
    // of sequence from the end of the buffer to the start of the buffer
    memcpy(
      fasta_reader->seq_buffer,
      fasta_reader->seq_buffer + fasta_reader->seq_buffer_index,
      w
    );
    fasta_reader->seq_buffer_index = 0;
    fasta_reader->seq_buffer_last = w;
    num_char_to_read = BUFFER_SIZE - w;
  }

  int num_copied = 0;
  if (!fasta_reader->at_end_of_seq && !fasta_reader->at_end_of_file) {

    // Read more sequence into the tmp buffer from the sequence file.
    char raw_buffer[BUFFER_SIZE];
    long start_file_pos = ftell(fasta_reader->fasta_file);
    size_t num_char_read = fread(
      raw_buffer,
      sizeof(char), 
      num_char_to_read,
      fasta_reader->fasta_file
    );

    if (feof(fasta_reader->fasta_file)) {
       fasta_reader->at_end_of_file = TRUE;
    }
    else if (num_char_read < num_char_to_read) {
      die(
        "Error while reading sequence from file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }
  
    // Copy from tmp buffer to seq buffer
    int i;
    for(i = 0; i < num_char_read; ++i) {
      char c = raw_buffer[i];
      assert(c != 0);
      if (isspace(c)) {
        // Skip over white space
        if (c == '\n') {
          fasta_reader->at_start_of_line = TRUE;
        }
        else {
          fasta_reader->at_start_of_line = FALSE;
        }
      }
      else if (c == '>' && fasta_reader->at_start_of_line == TRUE) {
        // We found the start of a new sequence while trying
        // to fill the buffer. Leave the buffer incomplete.
        // and wind back the file
        fseek(fasta_reader->fasta_file, start_file_pos + i - 1, SEEK_SET);
        fasta_reader->at_end_of_seq = TRUE;
        fasta_reader->at_start_of_line = FALSE;
        fasta_reader->at_end_of_file = FALSE;
        break;
      }
      else {
        fasta_reader->at_start_of_line = FALSE;
        // Check that character is legal in alphabet. 
        // If not, replace with wild card character.
        if (alph_is_known(fasta_reader->alphabet, c)) {
          fasta_reader->seq_buffer[fasta_reader->seq_buffer_last] = c;
        }
        else{
          fasta_reader->seq_buffer[fasta_reader->seq_buffer_last] 
            = alph_wildcard(fasta_reader->alphabet);
          fprintf(
            stderr, 
            "Warning: %c is not a valid character in %s alphabet.\n"
            "         Converting %c to %c.\n",
            c,
            alph_name(fasta_reader->alphabet),
            c,
            fasta_reader->seq_buffer[fasta_reader->seq_buffer_last]
          );
        }
        ++fasta_reader->seq_buffer_last;
        ++num_copied;
      }
    }
  }

  // Update data block with updated pointer into seq buffer
  set_sequence_in_data_block(
    data_block,
    fasta_reader->seq_buffer + fasta_reader->seq_buffer_index
  );
  // Did we find enough characters to fill the data block?
  int num_char_remaining =
    fasta_reader->seq_buffer_last - fasta_reader->seq_buffer_index;
  if (num_char_remaining >= w) {
    set_num_read_into_data_block(data_block, w);
  }
  else {
    set_num_read_into_data_block(data_block, num_char_remaining);
  }
  ++fasta_reader->seq_buffer_index;
  ++fasta_reader->current_position;
  set_start_pos_for_data_block(data_block, fasta_reader->current_position);

  return num_char_remaining >= w;
}

/*
 if buffer is empty
   read from disk into tmp buffer
   set current position index to 0
   copy number of bytes read from tmp buffer into buffer starting at origin
     if start of new sequence found rewind file to point just before new sequence
     set set top position index
 else if current position in buffer >= (top of buffer  - width of block)
   copy from (top of buffer - width of block + 1) to origin of buffer
   read (BUFFER_SIZE - width of block) characters from disk into tmp buffer
   copy number of bytes read from tmp buffer into buffer starting at (origin + width of block)
     if encounter start of new sequence rewind file to point just before new sequence
     don't copy white space
   update top of buffer
 else 
   // current position in buffer is < (top of buffer - width of block)
   increase current position by 1
 update data block data pointer with pointer to current position in buffer
 if there is a full block's worth of data available
   return true
 else
   return fail
 */

/******************************************************************************
 * Fills in the next data block for the sequence. 
 * During the first call for the sequence it fills in the full data block.
 * On successive calls, shifts the sequence in the block down one position
 * and reads one more character.
 * 
 * Returns TRUE if it was able to completely fill the block, FALSE if 
 * the next sequence or EOF was reached before the block was filled.
 * Dies if other errors encountered.
 *****************************************************************************/
BOOLEAN_T get_next_data_block_from_seq_reader_from_fasta_old(
  DATA_BLOCK_READER_T *reader, 
  DATA_BLOCK_T *data_block
) {

  BOOLEAN_T result = FALSE;
  char *raw_seq = get_sequence_from_data_block(data_block);
  int block_size = get_block_size_from_data_block(data_block);
  int num_read = get_num_read_into_data_block(data_block);

  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  if (num_read == block_size) {
    // Block is alread full, shift all elements in the block down by one position
    // FIXME CEG: Inefficient, replace with circular buffer.
    memmove(raw_seq, raw_seq + 1, block_size - 1);
    num_read = block_size - 1;
    raw_seq[num_read] = 0;
  }

  long save_pos = ftell(fasta_reader->fasta_file);
  int c = 0;

  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {
    if (isspace(c)) {
      // Skip over white space
      if (c == '\n') {
        fasta_reader->at_start_of_line = TRUE;
      }
      else {
        fasta_reader->at_start_of_line = FALSE;
      }
      continue;
    }
    else if (c == '>' && fasta_reader->at_start_of_line == TRUE) {
      // We found the start of a new sequence while trying
      // to fill the block. Leave the block incomplete.
      c = ungetc(c, fasta_reader->fasta_file);
      if (ferror(fasta_reader->fasta_file)) {
        die(
          "Error while reading file:%s.\nError message: %s\n", 
          fasta_reader->filename,
          strerror(ferror(fasta_reader->fasta_file))
        );
      }
      raw_seq[num_read] = 0;
      break;
    }
    else {
      // Fill in another character in the block
      raw_seq[num_read] = c;
      // Check that character is legal in alphabet. 
      // If not, replace with wild card character.
      if (!alph_is_known(fasta_reader->alphabet, c)) {
        raw_seq[num_read] = alph_wildcard(fasta_reader->alphabet);
        fprintf(
          stderr, 
          "Warning: %c is not a valid character in %s alphabet.\n"
          "         Converting %c to %c.\n",
          c,
          alph_name(fasta_reader->alphabet),
          c,
          raw_seq[num_read]
        );
      }
      ++num_read;
      if (num_read == block_size) {
        // block is full
        result = TRUE;
        break;
      }
    }
  }

  if (c == EOF && ferror(fasta_reader->fasta_file)) {
    die(
      "Error while reading file:%s.\nError message: %s\n", 
      fasta_reader->filename,
      strerror(ferror(fasta_reader->fasta_file))
    );
  }

  ++fasta_reader->current_position;
  fasta_reader->prev_block_position = save_pos;
  set_start_pos_for_data_block(data_block, fasta_reader->current_position);
  set_num_read_into_data_block(data_block, num_read);
  return result;

}

/******************************************************************************
 * Sets the state of the reader to the positon before the last data block was
 * read
 * 
 * Returns TRUE if it was able to rewind the reader, FALSE otherwise.
 * Dies if other errors encountered.
 *****************************************************************************/
BOOLEAN_T unget_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader
) {

  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  if (fasta_reader->prev_block_position) {
    errno = 0;
    int e = fseek(fasta_reader->fasta_file, fasta_reader->prev_block_position, SEEK_SET);
    if (e) {
      die(
        "Error when trying to unget data block in prior PSP file: "
        "%s.\nError message: %s.\n", 
        fasta_reader->filename, 
        strerror(errno)
      );
    }
    else {
      --fasta_reader->current_position;
      fasta_reader->prev_block_position = 0;
      result = TRUE;
    }
  }
  return result;

}
