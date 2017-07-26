/*
 * motif_regexp.c
 *
 *  Created on: 22/09/2008
 *      Author: rob
 */
#include "motif_regexp.h"
#include "motif-spec.h"
#include <math.h>

#define MAX_MOTIF_WIDTH 100

#define xstr(s) str(s)
#define str(s) #s

#define PSEUDO 0.01
#define DEFAULT_SITES 100
ARRAYLST_T* read_regexp_file(
   ALPH_T *alph, // IN Alphabet
   char *filename, // IN Name of file containing regular expressions
   ARRAYLST_T *motifs // IN/OUT The retrieved motifs, allocated and returned if NULL
) 
{
  FILE*      motif_file; // file containing the motifs.
  const char * pattern = "%" xstr(MAX_MOTIF_ID_LENGTH) "s\t%" xstr(MAX_MOTIF_WIDTH) "s";
  char motif_name[MAX_MOTIF_ID_LENGTH + 1];
  char motif_regexp[MAX_MOTIF_WIDTH + 1];
  char *symbol;
  ARRAY_T* these_freqs;
  MOTIF_T* m;
  int i, j, sym_index, core_count;
  MTYPE freq;
  ARRAY_T *bg;

  //Set things to the defaults.
  if (motifs == NULL) {
    motifs = arraylst_create();
  }
  bg = get_uniform_frequencies(alph, NULL);

  // Open the given MEME file.
  if (open_file(filename, "r", TRUE, "motif", "motifs", &motif_file) == 0)
    exit(1);

  while (fscanf(motif_file, pattern, motif_name, motif_regexp) == 2) {
    // Now we:
    // 1. Allocate new motif
    // 2. Assign name
    // 3. Convert regexp into frequency table.
    m = mm_malloc(sizeof(MOTIF_T*));
    memset(m, 0, sizeof(MOTIF_T));
    set_motif_id(motif_name, strlen(motif_name), m);
    set_motif_strand('+', m);
    m->length = strlen(motif_regexp);
    // Store the evalue of the motif.
    m->evalue = 0;
    m->log_evalue = -HUGE_VAL;
    /* Store the alphabet size in the motif. */
    m->alph = alph_hold(alph);
    m->flags = MOTIF_HAS_AMBIGS | (alph_has_complement(alph) ? MOTIF_BOTH_STRANDS : 0);
    /* Allocate memory for the matrix. */
    m->freqs = allocate_matrix(m->length, alph_size_full(alph));
    init_matrix(0, m->freqs);
    m->url = strdup("");
    // extra values that we use the defaults
    m->num_sites = 0;
    m->trim_left = 0;
    m->trim_right = 0;
    //Set motif frequencies here.
    for (i = 0; i < m->length; i++) {
      sym_index = alph_index(alph, motif_regexp[i]);
      if (sym_index != -1) {
        core_count = alph_ncomprise(alph, sym_index);
        freq = 1.0 / core_count;
        for (j = 0; j < core_count; j++) {
          set_matrix_cell(i, alph_comprise(alph, sym_index, j), freq, m->freqs);
        }
      }
      calc_ambigs(alph, FALSE, get_matrix_row(i, m->freqs));
    }
    // Compute and store the motif complexity.
    m->complexity = compute_motif_complexity(m);
    // Compute the motif scores
    m->scores = convert_freqs_into_scores(m->alph, m->freqs, bg, DEFAULT_SITES, PSEUDO);
    // add the new motif to the list
    arraylst_add(m, motifs);
  }
  fclose(motif_file);
  free_array(bg);
  return motifs;
}

