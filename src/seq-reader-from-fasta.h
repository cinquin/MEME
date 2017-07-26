/******************************************************************************
 * FILE: seq-reader-from-fasta.h
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-17
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the public declarations for the concrete implementation of a 
 * data-block-reader UDT for reading sequence segments from a FASTA file.
 *****************************************************************************/

#ifndef SEQ_READER_FROM_FASTA_H
#define SEQ_READER_FROM_FASTA_H

#include "alphabet.h"
#include "data-block-reader.h"

/******************************************************************************
 * This function creates an instance of a data block reader UDT for reading
 * sequence segments from a FASTA file.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_seq_reader_from_fasta(
  BOOLEAN_T parse_genomic_coord, 
  ALPH_T* alph, 
  const char *filename
);

/******************************************************************************
 * This function parses a FASTA sequence header and returns the
 * the chromosome name and starting position if the header has the format
 *	"chrXX:start-end", where start and end are positive integers.
 *****************************************************************************/
BOOLEAN_T parse_genomic_coordinates_helper(
  char* header,           // FASTA sequence header
  char** chr_name,  	  // chromosome name ((chr[^:]))
  size_t * chr_name_len_ptr, // number of characters in chromosome name
  int * start_ptr,        // start position of sequence (chr:(\d+)-)
  int * end_ptr           // end position of sequence (chr:\d+-(\d+))
);

/****************************************************************************
 * Read raw sequence until a new sequence is encountered or too many letters
 * are read.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
BOOLEAN_T read_raw_sequence_from_reader(
   DATA_BLOCK_READER_T *fasta_reader, // Sequence source
   unsigned int max_chars, // Maximum chars in raw_sequence.
   char* raw_sequence // Pre-allocated sequence.
);

/****************************************************************************
 * Read up to max_chars letters of one sequence from a DATA_BLOCK_T readder
 * and copy them in to the raw sequence in the SEQ_T object starting at the
 * given buffer offset. 
 ****************************************************************************/
void read_one_fasta_segment_from_reader(
   DATA_BLOCK_READER_T *reader,
   size_t max_size,
   size_t offset,
   SEQ_T *sequence
);


size_t get_current_pos_from_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
#endif
