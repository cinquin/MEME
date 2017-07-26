/**************************************************************************
 * FILE: tomtom.c
 * AUTHOR: Shobhit Gupta, Timothy Bailey, William Stafford Noble
 * CREATE DATE: 02-21-06
 * PROJECT: TOMTOM
 * DESCRIPTION: Motif-Motif comparison using the score statistic.
 **************************************************************************/
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#ifndef WIFEXITED
#include <wait.h>
#endif

#include "array-list.h"
#include "dir.h"
#include "macros.h"
#include "motif.h"       // Needed for complement_dna_freqs.
#include "matrix.h"      // Routines for floating point matrices.
#include "array.h"       // Routines for floating point arrays.
#include "alphabet.h"    // The alphabet.
#include "motif-in.h"
#include "simple-getopt.h"
#include "projrel.h"
#include "pssm.h"
#include "fitevd.h"
#include "config.h"
#include "ceqlogo.h"
#include "io.h"
#include "qvalue.h"
#include "string-list.h"
#include "object-list.h"
#include "red-black-tree.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

#define BINS 100
#define LOGOHEIGHT 10 // height of an alignment of logos in cm.
#define MIN_TARGET_DATABASE_SIZE 50
#define SHIFT_QUANTILE 0.5  // UK: quantile for shifting scores 
// (this could move to an input variable)

#define STRINGIZER(arg) #arg
#define STR_VALUE(arg) STRINGIZER(arg)

static const char *XML_FILENAME = "tomtom.xml";
static const char *HTML_STYLESHEET = "tomtom-to-html.xsl";
static const char *HTML_FILENAME = "tomtom.html";
static const char *OLD_HTML_FILENAME = "old_tomtom.html";
static const char *TXT_FILENAME = "tomtom.txt";
const char* tomtom_dtd = 
"<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n";

// structure for storing command line options of Tomtom.
typedef struct {
  ARRAYLST_T *q_files;
  ARRAYLST_T *t_files;
  char *bg_file;
  char *cs_file; // columnwise scores file (DEV only option)
  RBTREE_T *ids;
  RBTREE_T *idxs;
  bool clobber;
  char *outdir;
  bool text_only;
  bool html;
  double sig_thresh;
  bool sig_type_q;
  char *dist_name;
  void (*dist_func)(MOTIF_T *, MOTIF_T *, ARRAY_T *, int, MATRIX_T **);
  bool dist_dna_only;
  bool dist_allow_zeros;
  bool internal;
  int min_overlap;
  double pseudo;
  bool rc;
  bool xalph;
  bool png;
  bool eps;
  bool ssc;
  bool complete_scores;
  time_t now;
} TOMTOM_T;

/**************************************************************************
 * An object for storing information about a file containing a set of motifs
 **************************************************************************/
typedef struct {
  int id; //an id assigned to the database based on the order it was listed
  char *source;
  char *name;
  time_t last_mod;
  int loaded;
  int excluded;
  int list_index;//the index that in that list of motifs that this db starts at
  int list_entries;//the number of entries in the list of motifs from this db
  RBTREE_T *matched_motifs;
} MOTIF_DB_T;

/**************************************************************************
 * An object for storing the list of unique motif lengths in the target motifs
 **************************************************************************/
typedef struct {
  int min;
  int max;
  int total;
  int unique;
  int *lookup;
  int *lengths;
} MLEN_T;

/**************************************************************************
 * An object for storing information about one query-target alignment.
 **************************************************************************/
typedef struct {
  MOTIF_DB_T *db; // target motif database
  MOTIF_T *target; // the target motif, may be the reverse complement
  MOTIF_T *original; // the non reverse complemented motif
  int offset;
  int overlap;
  double pvalue;
  double evalue;
  double qvalue;
} TOMTOM_MATCH_T;

/**************************************************************************
 * An object for identifying a target motif uniquely
 **************************************************************************/
typedef struct {
  int db_id;
  char *motif_id;
} TOMTOM_ID_T;

/**************************************************************************
 * Create a new query-target alignment object.
 * Initially, there is no q-value.
 **************************************************************************/
TOMTOM_MATCH_T* new_tomtom_match (MOTIF_DB_T *target_db, MOTIF_T *target,
    MOTIF_T *original, int offset, int overlap, double pvalue, double evalue) {
  TOMTOM_MATCH_T *match;
  match = mm_malloc(sizeof(TOMTOM_MATCH_T));
  match->db = target_db;
  match->target = target;
  match->original = original;
  match->offset = offset;
  match->overlap = overlap;
  match->pvalue = pvalue;
  match->evalue = evalue;
  match->qvalue = 1.0;
  return match;
}

/**************************************************************************
 * Compare 2 matches and put the best ones first.
 **************************************************************************/
static int tomtom_match_compare(const TOMTOM_MATCH_T **p1, const TOMTOM_MATCH_T **p2) {
  // 1st sort on pvalue (ascending)
  if ((*p1)->pvalue == (*p2)->pvalue) {
    // 2nd sort on overlap (descending)
    if ((*p1)->overlap == (*p2)->overlap) {
      // 3rd sort on target db (ascending)
      if ((*p1)->db->id == (*p2)->db->id) {
        // Finally sort on target motif id with strand which should be unique (ascending)
        return strcmp(get_motif_st_id((*p1)->target), get_motif_st_id((*p2)->target));
      } else if ((*p1)->db->id < (*p2)->db->id) {
        return -1;
      } else {
        return 1;
      }
    } else if ((*p1)->overlap > (*p2)->overlap) {
      return -1;
    } else {
      return 1;
    }
  } else if ((*p1)->pvalue < (*p2)->pvalue) {
    return -1;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Copy a target motif identifier
 **************************************************************************/
static void* copy_tomtom_id(void *p) {
  TOMTOM_ID_T *id, *id_copy;
  id = (TOMTOM_ID_T*)p;
  id_copy = mm_malloc(sizeof(TOMTOM_ID_T));
  id_copy->db_id = id->db_id;
  id_copy->motif_id = id->motif_id; // assumes won't be deallocated
  return id_copy;
}

/**************************************************************************
 * Compare 2 target motif identifiers
 **************************************************************************/
static int compare_tomtom_id(const void *p1, const void *p2) {
  TOMTOM_ID_T *id1, *id2;
  id1 = (TOMTOM_ID_T*)p1;
  id2 = (TOMTOM_ID_T*)p2;
  if (id1->db_id == id2->db_id) {
    return strcmp(id1->motif_id, id2->motif_id);
  } else if (id1->db_id < id2->db_id) {
    return -1;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Create a new motif database entry.
 * Makes a copy of the strings so there is no confusion on what can be
 * freed later.
 **************************************************************************/
static MOTIF_DB_T* create_motif_db(int id, char *source, char *name, time_t last_mod) {
  MOTIF_DB_T *db;

  assert(source != NULL);
  assert(name != NULL);

  db = mm_malloc(sizeof(MOTIF_DB_T));
  memset(db, 0, sizeof(MOTIF_DB_T));
  db->id = id;
  copy_string(&(db->source), source);
  copy_string(&(db->name), name);
  db->last_mod = last_mod;
  // the strings that I will be using are fixed size buffers of the MOTIF_T object
  // and they won't be deallocated in the life of the database entry so I don't need
  // a key copy or key free function.
  db->matched_motifs = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);

  return db;
}

/**************************************************************************
 * Destroys the motif database entry.
 **************************************************************************/
void destroy_motif_db(void *ptr) {
  MOTIF_DB_T *db = (MOTIF_DB_T*)ptr;
  rbtree_destroy(db->matched_motifs);
  free(db->source);
  free(db->name);
  memset(db, 0, sizeof(MOTIF_DB_T));
  free(db);
}

/*****************************************************************************
 * outputs the file name when given a path to that file.
 *
 * searches the path string for '/' chars and returns the pointer to the position
 * after the last '/'. Note that the path is not copied or modified in any way.
 * This is not intended to be perfect as I'm sure there's some valid way to have 
 * a '/' in a file name but I expect it will work most of the time.
 *****************************************************************************/
static char* file_name_from_path(char *path) {
  char *where, *start, *end, *name, *c;
  int i, len;
  // scan forwards and find the last '/'
  start = path;
  for (where = path; *where != '\0'; where++) {
    if (*where == '/') start = where+1;
  }
  // now work backwards and find the last '.' after the last '/'
  end = where;
  for (; where > start; where--) {
    if (*where == '.') {
      end = where;
      break;
    }
  }
  len = end - start;
  name = mm_malloc(sizeof(char) * (len + 1));
  for (i = 0, c = start; c < end; c++, i++) {
    name[i] = (*c == '_' ? ' ' : *c);
  }
  name[len] = '\0';
  return name;
}

/**************************************************************************
 * Creates motif database entries from the file names.
 **************************************************************************/
static ARRAYLST_T* create_motif_dbs(ARRAYLST_T *files, char *type, time_t now, bool *stdin_used) {
  ARRAYLST_T *dbs;
  char *path, *name;
  time_t last_mod;
  struct stat stbuf;
  int i, result;

  dbs = arraylst_create_sized(arraylst_size(files));
  for (i = 0; i < arraylst_size(files); i++) {
    path = (char*)arraylst_get(i, files);
    if (strcmp(path, "-") == 0) {
      if (*stdin_used) {
        die("Database file %s %d made use of standard input when it is "
            "already in use for another file.\n", type, i+1);
      }
      *stdin_used = true;
      name = strdup("standard input");
      last_mod = now;
    } else {
      // check if the path exists and determine last modified time
      errno = 0;
      if (stat(path, &stbuf) == -1) {
        if (errno == ENOENT) { // File not found
          die("Database file %s %d points to the path \"%s\" which does "
              "not exist.\n", type, i+1, path);
        } else { // stat failed for some other reason
          die("Unable to check for status of file '%s'.\nError: %s.\n",
            path, strerror(errno));
        }
      }
      last_mod = stbuf.st_mtime;
      // get the file name from the path
      name = file_name_from_path(path);
    }
    arraylst_add(create_motif_db(i, path, name, last_mod), dbs);
    free(name);
  }
  return dbs;
}

/**************************************************************************
 * Creates a logo of the aligned target and query motif
 **************************************************************************/
static void create_logo(MOTIF_DB_T *target_db, MOTIF_T *target,
    MOTIF_DB_T *query_db, MOTIF_T *query, int offset, 
    bool png, bool eps, bool ssc, char *output_dirname) {
  STR_T *logo_filename;
  char *logo_path;
  // create filename
  logo_filename = str_create(50);
  str_setf(
      logo_filename, 
      "align_%s_%d_%s", 
      get_motif_id(query), 
      target_db->id, 
      get_motif_st_id(target)
    );
  // create output file path
  logo_path = make_path_to_file(output_dirname, str_internal(logo_filename));
  // create logo
  CL_create2(target, get_motif_id(target), query, get_motif_id(query), true,
      ssc, LOGOHEIGHT, 0, -offset, "Tomtom", logo_path, eps, png);
  // clean up
  myfree(logo_path);
  str_destroy(logo_filename, false);
}


/**************************************************************************
 * Create a lookup table for the pv.
 *
 * FIXME: document what this really does!
 *
 **************************************************************************/
void get_pv_lookup_new(MATRIX_T* scores, int query_len, int targets_len, int range,
    MATRIX_T* reference_matrix, MATRIX_T* pv_lookup_matrix, MATRIX_T* pmf_matrix) {
  int i, j, k, s, ii, size, start, address_index, max;
  double uniform_bg, old, new, p;
  ARRAY_T *pdf_old, *pdf_new, *pv;

  address_index = 0; // next free position in matrix
  size = (query_len * range) + 1;
  pdf_old = allocate_array(size);
  pdf_new = allocate_array(size);
  // Uniform background for the infinite alphabet
  uniform_bg = 1.0 / targets_len;

  // Compute pv tables for groups of motif columns [start...]
  // UK: "start" is the first query column in the overlap
  for (start = query_len - 1; start >= 0; start--) {

    // Compute the pdf recursively.
    init_array(0, pdf_new);
    set_array_item(0, 1, pdf_new);    // Prob(0)
    // UK: "start + i" is the last query column in the overlap 
    // (so "i+1" is the overlap size)
    for (i = 0; (start + i) < query_len; i++) {
      max = i * range;
      SWAP(ARRAY_T*, pdf_new, pdf_old);
      init_array(0, pdf_new); // UK: this should be safer
      // UK: "j" runs over the target column indices
      for (j = 0; j < targets_len; j++) { 
        s = (int) get_matrix_cell(start + i, j, scores);
        for (k = 0; k <= max; k++) {
          old = get_array_item(k, pdf_old);
          if (old != 0) {
            new = get_array_item(k+s, pdf_new) + (old * uniform_bg);
            set_array_item(k+s, new, pdf_new);
          } // old
        } // k
      } // j

      // Compute 1-cdf for motif consisting of columns [start, start+i]
      // This is the p-value lookup table for those columns.
      pv = allocate_array(get_array_length(pdf_new));
      copy_array(pdf_new, pv);
      for (ii = size - 2; ii >= 0; ii--) {
        p = get_array_item(ii, pv) + get_array_item(ii + 1, pv);
        set_array_item(ii, MIN(1.0, p), pv);
      }

      // copy the pv lookup table into the lookup matrix
      set_matrix_row(address_index, pv, pv_lookup_matrix);
      free_array(pv);
      // UK: copy the pmf into the pmf matrix
      set_matrix_row(address_index, pdf_new, pmf_matrix);

      // Store the location of this pv lookup table in the reference array
      set_matrix_cell(start, start + i, (double) address_index, reference_matrix);
      address_index++;
    } // i
  } // start position

  // Free space.
  free_array(pdf_new);
  free_array(pdf_old);
}

/**************************************************************************
 * Parse a list of Tomtom output strings, converting p-values to q-values.
 **************************************************************************/
static void convert_tomtom_p_to_q(int match_count, TOMTOM_MATCH_T **match_list) {
  ARRAY_T *pvalues;
  int i;

  // Extract all of the p-values.
  pvalues = allocate_array(match_count);
  for (i = 0; i < match_count; i++) set_array_item(i, match_list[i]->pvalue, pvalues);

  // Convert p-values to q-values.
  compute_qvalues(
    false, // Don't stop with FDR.
    true, // Estimate pi-zero.
    NULL, // Don't store pi-zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    match_count, 
    pvalues,
    NULL // No sampled p-values provided
  );

  // Set q-values
  for (i = 0; i < match_count; i++) match_list[i]->qvalue = get_array_item(i, pvalues);
  // clean up
  free_array(pvalues);
}

/**************************************************************************
 * Gets the consensus sequence from a motif.
 * Caller is responsible for deallocating memory.
 **************************************************************************/
static char* get_cons(MOTIF_T *motif) {
  ALPH_T *alph;
  MATRIX_T *freqs;
  char *cons;
  int i, j, len;
  // use float to maintain old rounding behaviour
  float max;
  alph = get_motif_alph(motif);
  freqs = get_motif_freqs(motif);
  len = get_motif_length(motif);
  cons = (char*)mm_malloc(sizeof(char) * (len + 1));
  for (i = 0; i < len; i++) {
    cons[i] = alph_char(alph, 0);
    max = get_matrix_cell(i, 0, freqs);
    //find the column with the best score
    for (j = 1; j < alph_size_core(alph); j++) {
      if (get_matrix_cell(i, j, freqs) > max) {
        max = get_matrix_cell(i, j, freqs);
        cons[i] = alph_char(alph, j);
      }
    }
  }
  cons[len] = '\0';
  return cons;
}


/**************************************************************************
 * Computes the score for a particular configuration
 **************************************************************************/
static double compute_overlap_score(
    MATRIX_T* scores, int query_len, int target_len, int target_offset,
    int configuration, int *index_start, int *index_stop) {
  int target_col, query_col;
  double score;
  // Initialize the columns from where the comparison needs to be started
  if (configuration <= query_len) {
    target_col = 0;
    query_col = query_len - configuration;
  } else {
    query_col = 0;
    target_col = configuration - query_len;
  }
  // Sum the pairwise column scores.
  score = 0;
  *index_start = query_col;
  do {
    score += get_matrix_cell(query_col, target_offset + target_col, scores);
    query_col++; target_col++;
  } while ((query_col < query_len) && (target_col < target_len));
  *index_stop = query_col - 1; // -1 for the last increment
  return score;
}

/**************************************************************************
 * Computes the optimal offset and score corresponding to the smallest
 * pvalue.
 **************************************************************************/
static void compare_motifs (
    MATRIX_T *scores, int query_len, int target_len, int target_offset,
    MATRIX_T *reference_matrix, MATRIX_T *pv_lookup_matrix,
    bool internal, int min_overlap, bool complete_scores,
    double scores_scale, double scores_offset, int *configurations,
    int *optimal_offset, double *optimal_pvalue, int *optimal_overlap) {
  int start, end, mo, configuration;
  // determine acceptable overlaps
  start = 1;
  end = query_len + target_len - 1;
  if (internal) {
    // motifs must be fully overlapping
    if (query_len < target_len) {
      start = query_len; end = target_len;
    } else {
      start = target_len; end = query_len;
    }
  } else if (min_overlap > 1) {
    // Override min-overlap if query/target is smaller than min-overlap
    mo = min_overlap;
    if (query_len < mo) mo = query_len;
    if (target_len < mo) mo = target_len;
    start = mo;
    end = query_len + target_len - mo;
  }
  // Slide one profile over another
  *configurations = end - start + 1;
  *optimal_pvalue = BIG;
  *optimal_overlap = -1;
  *optimal_offset = 0; // stop compiler complaining
  for (configuration = start; configuration <= end; configuration++) {
    double pvalue, temp_score;
    int overlap, i_start, i_stop, pv_index;
    // Compute the score for the current configuration
    temp_score = compute_overlap_score(scores, query_len, target_len,
        target_offset, configuration, &i_start, &i_stop);
    // Compute the pvalue of the score
    overlap = i_stop - i_start + 1;
    assert(overlap > 0);
    if (complete_scores) {
      // UK: undo the offset when (quantile-) shifted scores are used, 
      // also pvalue is not really probability in this case
      // JIJ: I assume the negative is to make it so the smaller scores are the
      // better ones. This seems to be the case as it is reversed upon return...
      // I had originally assumed that the
      // "temp_score + overlap * scores_offset * scores_scale" was meant to
      // reverse the effect of scaling the matrix to go between 0 and bins but
      // the formula doesn't seem right for that. I would have expected:
      // scores_scale * temp_score + overlap * scores_offset
      // I wonder if it was a mistake?
      // Another posibility is that this is part of the method used to take into
      // account non-overlapping columns, as before calculating p-values
      // in calc_shifted_scores_pvalues_per_query it does:
      // "temp_score - query_len * scores_offset * scores_scale"
      // as query_len is always going to be larger than overlap this results
      // in subracting away "scores_offset * scores_scale" for every
      // non-overlapping query column. However the odd thing about that theory
      // is that the full correction isn't applied here.
      pvalue = - rint(temp_score  + overlap * scores_offset * scores_scale);
    } else {
      pv_index = (int)get_matrix_cell(i_start, i_stop, reference_matrix);
      pvalue = get_matrix_cell(pv_index, (int)temp_score, pv_lookup_matrix);
    }
    // Keep the configuration with the lowest pvalue (or -score).
    // If pvalue is the same, keep the largest overlap 
    if (pvalue == *optimal_pvalue) {
      if (*optimal_overlap < overlap) {
        *optimal_offset = configuration - query_len;
        *optimal_pvalue = pvalue;
        *optimal_overlap = overlap;
      }
    } else {
      if (pvalue < *optimal_pvalue) {
        *optimal_offset = configuration - query_len;
        *optimal_pvalue = pvalue;
        *optimal_overlap = overlap;
      }
    }
  }
  // undo the negation trick used to make the best scores the smallest
  if (complete_scores) *optimal_pvalue = -(*optimal_pvalue);
}

/**************************************************************************
 * Computes ALLR scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_allr.
 **************************************************************************/
void allr_scores(
    MOTIF_T* query_motif,
    MOTIF_T* target_motif,
    ARRAY_T* background,
    int offset,
    MATRIX_T** score
    ) {
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_DENOMINATOR = 0;
      int index; /* Index corresponds to the alphabets. */

      for (index = 0; index < get_motif_alph_size(query_motif); index++) {
        /* The numerator of the ALLR */
        double nr_part1, nr_part2;
        /*  Part 1 of the numerator  */
        nr_part1
          = (
              get_motif_nsites(target_motif)
              * (
                get_matrix_cell(index_target, index, get_motif_freqs(target_motif))
                /* nb for target */
                )
              * log(
                get_matrix_cell(index_query, index, get_motif_freqs(query_motif))
                / get_array_item(index, background)
                )
            ); /* Likelihood for query */
        /*  Part 2 of the numerator */
        nr_part2 
          = (
              get_motif_nsites(query_motif)
              * (
                get_matrix_cell(index_query, index, get_motif_freqs(query_motif))
                /* nb for query */
                )
              * log(
                get_matrix_cell(index_target, index, get_motif_freqs(target_motif))
                / get_array_item(index, background)
                )
            ); /* Likelihood for target */

        /* Sum of the two parts. */
        SCORE_NUMERATOR += nr_part1 + nr_part2;

        /*  The denominoator of the ALLR */
        SCORE_DENOMINATOR
          += (
              (
               get_motif_nsites(query_motif)
               *  get_matrix_cell(
                 index_query, 
                 index, 
                 get_motif_freqs(query_motif)
                 )
               // nb for query
              )
              + (
                get_motif_nsites(target_motif)
                * get_matrix_cell(index_target, index, get_motif_freqs(target_motif))
                ) // nb for target
             ); /* Sum of nb for query and nb for target */

      } /* Nr and Dr of the score summed over the entire alphabet. */
      SCORE_F = SCORE_NUMERATOR / SCORE_DENOMINATOR;
      set_matrix_cell(index_query, offset + index_target, SCORE_F, *score);

    } /* Target motif */
  } /* Query motif */
} /* Function ALLR scores */

/**************************************************************************
 * Computes euclidian distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void ed_scores(
    MOTIF_T* query_motif,
    MOTIF_T* target_motif,
    ARRAY_T* background,
    int offset,
    MATRIX_T** score
    ) {

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < get_motif_alph_size(query_motif); index++) {

        double query_freq = get_matrix_cell(index_query, index, get_motif_freqs(query_motif));
        double target_freq = get_matrix_cell(index_target, index, get_motif_freqs(target_motif));
        /* Euclidean distance */
        double sq;
        sq = pow((query_freq - target_freq), 2);
        SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = - sqrt(SCORE_T);
      set_matrix_cell(index_query, offset + index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function ed_scores */

/**************************************************************************
 * Computes sandelin distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void sandelin_scores(
    MOTIF_T* query_motif,
    MOTIF_T* target_motif,
    ARRAY_T* background,
    int offset,
    MATRIX_T** score
    ) {

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < get_motif_alph_size(query_motif); index++) {

        double query_freq
          = get_matrix_cell(index_query, index, get_motif_freqs(query_motif));
        double target_freq 
          = get_matrix_cell(index_target, index, get_motif_freqs(target_motif));
        /* Euclidean distance */
        double sq;
        sq = pow((query_freq - target_freq), 2);
        SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = 2 - SCORE_T;
      set_matrix_cell(index_query, offset + index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function sandelin*/

/**************************************************************************
 * Computes Kullback-Leiber for all columns of the query matrix against target
 * motif. Note we use
 * a modified Kullback-leibler for computing single column scores. We do not
 * divide the score by the width. This division is replaced by our DP based
 * p-value computation.
 **************************************************************************/
void kullback_scores(
    MOTIF_T* query_motif,
    MOTIF_T* target_motif,
    ARRAY_T* background,
    int offset,
    MATRIX_T** score
    ) {
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {
      double SCORE_F = 0;
      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < get_motif_alph_size(query_motif); index++) {
        double query_freq
          = get_matrix_cell(index_query, index, get_motif_freqs(query_motif));
        double target_freq
          = get_matrix_cell(index_target, index, get_motif_freqs(target_motif));
        /* Average Kullback */
        double avg_kull;
        avg_kull = ((query_freq * log10(query_freq/target_freq))
            + (target_freq * log10(target_freq/query_freq))) / 2;
        SCORE_F = SCORE_F + avg_kull;
      }
      SCORE_F = SCORE_F * -1;
      set_matrix_cell(index_query, offset + index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function Kullback*/


/**************************************************************************
 * Computes Pearson scores for all columns of the query matrix against target
 * motif. 
 **************************************************************************/
void pearson_scores(
    MOTIF_T* query_motif,
    MOTIF_T* target_motif,
    ARRAY_T* background,
    int offset,
    MATRIX_T** score
    ) {

  /*   Define the average values */
  double x_bar = 1.0 / get_motif_alph_size(query_motif);
  double y_bar = 1.0 / get_motif_alph_size(target_motif);

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_SQ_DENOMINATOR_1 = 0;
      double SCORE_SQ_DENOMINATOR_2 = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < get_motif_alph_size(query_motif); index++) {
        /* The numerator of the Pearson */
        SCORE_NUMERATOR
          += (get_matrix_cell(index_query, index, get_motif_freqs(query_motif)) - x_bar)
          * (get_matrix_cell(index_target, index, get_motif_freqs(target_motif)) - y_bar);

        /*  The denominoator of the Pearson */
        SCORE_SQ_DENOMINATOR_1
          += pow(
              (get_matrix_cell(index_query, index, get_motif_freqs(query_motif))
               - x_bar),
              2
              );
        SCORE_SQ_DENOMINATOR_2
          += pow(
              (get_matrix_cell(index_target, index, get_motif_freqs(target_motif))
               - y_bar),
              2
              );
      } /* Nr and Dr components summed over the entire alphabet. */

      /*       Compute the Denominator */
      double SCORE_SQ_DENOMINATOR 
        = SCORE_SQ_DENOMINATOR_1 * SCORE_SQ_DENOMINATOR_2;
      double SCORE_DENOMINATOR = sqrt(SCORE_SQ_DENOMINATOR);

      SCORE_F = SCORE_DENOMINATOR != 0 ? (SCORE_NUMERATOR / SCORE_DENOMINATOR) : 0;
      set_matrix_cell(index_query, offset + index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function pearson */

/**************************************************************************
 * Estimate one Dirichlet component.
 *
 * {\hat p}_i
 *     = \frac{n_i + \alpha_i}{\sum_{j \in \{A,C,G,T\}} (n_j + \alpha_j)}
 **************************************************************************/
static ARRAY_T* one_dirichlet_component(
    int      hyperparameters[],
    ARRAY_T* counts
    ) {
  int alph_size = get_array_length(counts);
  ARRAY_T* return_value = allocate_array(alph_size);

  // Calculate the denominator.
  double denominator = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    denominator += get_array_item(i_alph, counts)
      + (int)(hyperparameters[i_alph]);
  }

  // Estimate the source distribution.
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    double numerator = get_array_item(i_alph, counts) 
      + hyperparameters[i_alph];
    set_array_item(i_alph, numerator / denominator, return_value);
  }
  return(return_value);
}

/**************************************************************************
 * Likelihood of the data, given a Dirichlet component.
 *
 * This is the final equation on p. 13 of Habib et al.
 **************************************************************************/
static double compute_log_likelihood (int hypers[], ARRAY_T* counts) {

  // Sum the counts and the hyperparameters.
  double hyper_sum = 0.0;
  double count_sum = 0.0;
  int i_alph;
  int alphabet_size = get_array_length(counts);
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    hyper_sum += (float)(hypers[i_alph]);
    count_sum += get_array_item(i_alph, counts);
  }

  // Compute the log of likelihood.
  double log_likelihood = lgamma(hyper_sum)- lgamma(count_sum+hyper_sum);
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    double hyper = (float)(hypers[i_alph]);
    double count = get_array_item(i_alph, counts);
    log_likelihood += lgamma(hyper + count) - lgamma(hyper);
  }
  return(log_likelihood);
}  

/********************************************
 * log(a+b)=loga+log(1+exp(logb-loga))
 *********************************************/
static double logsum(double log_a, double log_b) {
  return (log_a < log_b ?
      log_b + log(1.0 + exp(log_a - log_b)) : log_a + log(1.0 + exp(log_b - log_a)));
}

/**************************************************************************
 * Estimate the source distribution from counts, using a 1-component
 * or 5-component Dirichlet prior.
 *
 * The 5-component prior only works for DNA.
 **************************************************************************/
  static ARRAY_T* estimate_source
(int      num_dirichlets,
 ARRAY_T* counts)
{
  int dna_hypers[6][4] = {
    {1, 1, 1, 1}, // Single component
    {5, 1, 1, 1}, // Five components ...
    {1, 5, 1, 1},
    {1, 1, 5, 1},
    {1, 1, 1, 5},
    {2, 2, 2, 2}
  };
  int protein_hypers[20] 
    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  // Single component or 5-component mixture?
  int alph_size = get_array_length(counts);
  if (num_dirichlets == 1) {
    if (alph_size == 4) {
      return(one_dirichlet_component(dna_hypers[0], counts));
    } else if (alph_size == 20) {
      return(one_dirichlet_component(protein_hypers, counts));
    } else {
      die("Invalid alphabet size (%d) for BLiC score.\n", alph_size);
    }
  }

  // Die if trying to do proteins.
  if (alph_size != 4) {
    die("Sorry, BLiC5 is only implemented for DNA motifs.\n");
  }

  // Compute the denominator.

  int i_hyper;
  double loga=log(1.0/(float) num_dirichlets)+compute_log_likelihood(dna_hypers[1], counts);
  double logb = 0.0;
  double log_sum = 0.0;
  for (i_hyper = 2; i_hyper <= num_dirichlets; i_hyper++) {
    logb=log(1.0/(float) num_dirichlets)+compute_log_likelihood(dna_hypers[i_hyper], counts);
    log_sum = logsum(loga,logb);
    loga = log_sum; 
  }

  double log_denominator=log_sum;
  ARRAY_T* return_value = allocate_array(alph_size);
  for (i_hyper = 1; i_hyper <= num_dirichlets; i_hyper++) {

    /* 
       Compute the posterior.
       \Pr(k|n) = \frac{\Pr(k) \Pr(n|k)}{\sum_j \Pr(j) \Pr(n|j)}
       */
    double log_posterior = log(1.0 / (float) num_dirichlets) 
      + compute_log_likelihood(dna_hypers[i_hyper], counts) - log_denominator;

    // Compute the Dirichlet component.
    ARRAY_T* one_component
      = one_dirichlet_component(dna_hypers[i_hyper], counts);

    // Multiply them together and add to the return value.
    scalar_mult(exp(log_posterior), one_component);
    sum_array(one_component, return_value);

    free_array(one_component);
  }

  return(return_value);
}

/**************************************************************************
 * Computes the BLiC (Bayesian Likelihood 2-Component) score described
 * by Equation (2) in Habib et al., PLoS CB 4(2):e1000010.
 *
 * Thus far, BLiC is only implemented for DNA.  It would be relatively
 * straightforward to extend to amino acids.
 **************************************************************************/
static double one_blic_score(
    int      num_dirichlets,
    ARRAY_T* query_counts, 
    ARRAY_T* target_counts,
    ARRAY_T* background
    ) {

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += 2 * get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }

  // Fourth term.
  double fourth_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    fourth_term += get_array_item(i_alph, common_counts)
      * log(get_array_item(i_alph, background));
  }

  return(first_term - (second_term + third_term + fourth_term));

}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 **************************************************************************/
static void blic_scores(
    int num_dirichlets,     // Number of Dirichlets [IN]
    MOTIF_T* query_motif,   // Single query motif [IN].
    MOTIF_T* target_motif,  // Single target motif [IN].
    ARRAY_T* background,    // Background distribution [IN].
    int offset,
    MATRIX_T** scores)      // BLiC score matrix, indexed by start. [OUT]
{
  // Traverse the query and target motifs
  int query_i;
  for (query_i = 0; query_i < get_motif_length(query_motif); query_i++) {
    ARRAY_T* query_counts = get_motif_counts(query_i, query_motif);

    int target_i;
    for (target_i = 0; target_i < get_motif_length(target_motif); target_i++) {
      ARRAY_T* target_counts = get_motif_counts(target_i, target_motif);

      // Compute the BLiC score according to Equation (2).
      double my_score 
        = one_blic_score(
            num_dirichlets,
            query_counts,
            target_counts,
            background
            );
      free_array(target_counts);

      // Store the computed score.
      set_matrix_cell(query_i, offset + target_i, my_score, *scores);
    }
    free_array(query_counts);
  }
}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 * The number of dirichlets is set to 1.
 **************************************************************************/
static void blic1_scores(
    MOTIF_T* query_motif,   // Single query motif [IN].
    MOTIF_T* target_motif,  // Single target motif [IN].
    ARRAY_T* background,    // Background distribution [IN].
    int offset,
    MATRIX_T** scores       // BLiC score matrix, indexed by start. [OUT]
    ) {
  blic_scores(1, query_motif, target_motif, background, offset, scores);
}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 * The number of dirichlets is set to 5.
 **************************************************************************/
static void blic5_scores(
    MOTIF_T* query_motif,  // Single query motif [IN].
    MOTIF_T* target_motif, // Single target motif [IN].
    ARRAY_T* background,   // Background distribution [IN].
    int offset,
    MATRIX_T** scores      // BLiC score matrix, indexed by start. [OUT]
    ) {
  blic_scores(5, query_motif, target_motif, background, offset, scores);
}

/**************************************************************************
 * Computes the llr score 
 **************************************************************************/
static double one_llr_score(
    int      num_dirichlets,
    ARRAY_T* query_counts, 
    ARRAY_T* target_counts,
    ARRAY_T* background
    ) {

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }

  return(first_term - (second_term + third_term));

}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 **************************************************************************/
static void llr_scores(
    int num_dirichlets,     // Number of Dirichlets [IN]
    MOTIF_T* query_motif,   // Single query motif [IN].
    MOTIF_T* target_motif,  // Single target motif [IN].
    ARRAY_T* background,    // Background distribution [IN].
    int offset,
    MATRIX_T** scores       // BLiC score matrix, indexed by start. [OUT]
    ) {
  // Traverse the query and target motifs
  int index_query;
  for (index_query = 0; index_query < get_motif_length(query_motif); index_query++) {
    ARRAY_T* query_counts = get_motif_counts(index_query, query_motif);

    int index_target;
    for (index_target = 0; index_target < get_motif_length(target_motif); index_target++) {
      ARRAY_T* target_counts = get_motif_counts(index_target, target_motif);


      double my_score 
        = one_llr_score(
            num_dirichlets,
            query_counts,
            target_counts,
            background
            );
      free_array(target_counts);

      // Store the computed score.
      set_matrix_cell(index_query, offset + index_target, my_score, *scores);
    }
    free_array(query_counts);
  }
}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 * The number of Dirichlets is set to 1.
 **************************************************************************/
  static void llr1_scores
(MOTIF_T* query_motif,  // Single query motif [IN].
 MOTIF_T* target_motif, // Single target motif [IN].
 ARRAY_T* background,   // Background distribution [IN].
 int offset,
 MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  llr_scores(1, query_motif, target_motif, background, offset, scores);
}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 * The number of Dirichlets is set to 5.
 **************************************************************************/
  static void llr5_scores
(MOTIF_T* query_motif,  // Single query motif [IN].
 MOTIF_T* target_motif, // Single target motif [IN].
 ARRAY_T* background,   // Background distribution [IN].
 int offset,
 MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  llr_scores(5, query_motif, target_motif, background, offset, scores);
}

/**************************************************************************
 * Compare scores and sort ascending.
 **************************************************************************/
static int score_compare(const double *score1, const double *score2) {
  return ((*score1 < *score2) ? -1 : ((*score1 == *score2) ? 0 : 1));
}

/**************************************************************************
 * UK: Scores are shifted by shift quantile
 *
 * JIJ:
 * This method shifts the scores by the -median score so the median score
 * becomes zero.
 * The idea here is that we don't want Tomtom to align on uninformative
 * columns when there are better options.
 * To do that we want to give a score to unaligned columns so that failing
 * to align an uninformative column gets a higher score than failing to align
 * an informative one.
 * The suggestion from the paper is that the median of the set of scores from
 * aligning a particular query column to every possible target column will
 * forfill this critera.
 * Paraphrased from the paper...
 * Assuming a target column 'T' from the set 'T_all' of all possible target
 * columns in a database, a query column 'A' which represents strong preference
 * for a single letter, a query column 'U' which represents a uniform
 * preference and a query column 'X' which represents a preference for
 * the letters not preferred by 'A'.
 * Also assuming a scoring function S(C_1,C_2) which compares a
 * column C_1 with a column C_2.
 * Then:
 *  max{S(U, T) : T in T_all} = max{S(A, T) : T in T_all}
 * And:
 *  S(A, X) < S(U, X)
 * So it follows that:
 *  min{S(U, T) : T in T_all} > min{S(A, T) : T in T_all}
 * This explains why the median of the set
 *  {S(Q, T) : T in T_all}
 * will be higher when Q = U (or a uninformative column) then when Q = A (or an
 * informative column).
 * To make things even simple the paper proposes that we can get around having
 * to add scores for each unaligned column by making the score to be added
 * (aka the median score for the column) equal to zero.
 **************************************************************************/
static void shift_all_pairwise_scores(MATRIX_T* scores,
    int query_len, int targets_len, double shift_quant) {
  int t, q, index_1, index_2;
  double *one_col, quantile;
  // allocate space
  one_col = malloc(targets_len * sizeof(double));
  // determine the indexes we'll use for calculating the quantile
  // Note that by default shift_quant is set to 0.5 so we're selecting
  // the median here (well not quite - it's actually off by 1 position,
  // but close enough that I don't want to change the behaviour).
  index_1 = (int) floor(shift_quant * targets_len) - 1;
  index_2 = (int) ceil(shift_quant * targets_len) - 1;
  // handle edge case of very small target length
  if (index_1 < 0) index_1 = 0;
  if (index_2 < index_1) index_2 = index_1;
  // loop over each query column adjusting scores
  for (q = 0; q < query_len; ++q) {
    // copy the scores into a sortable form
    for (t = 0; t < targets_len; ++t) {
      one_col[t] = get_matrix_cell(q, t, scores);
    }
    // There are probably faster ways to find the scores at the indexes, see:
    // http://en.wikipedia.org/wiki/Selection_algorithm
    qsort(one_col, targets_len, sizeof(double), (void *)score_compare);
    quantile = (one_col[index_1] + one_col[index_2]) / 2.0;
    // shift the column
    for (t = 0; t < targets_len; ++t) {
      incr_matrix_cell(q, t, -quantile, scores);
    }
  }
  // cleanup
  free(one_col);
}

/**************************************************************************
 * Print the score matrix to the specified file.
 **************************************************************************/
static void print_columnwise_scores(char *cs_file, MATRIX_T *cs) {
  FILE* out;
  fprintf(stderr, "Storing %d by %d matrix of column scores in %s.\n",
      get_num_rows(cs), get_num_cols(cs), cs_file);
  if (!open_file(cs_file, "w", false, "column scores", "column scores", &out)) {
    exit(1);
  }
  print_matrix(cs, 10, 7, false, out);
  fclose(out);
}

/*****************************************************************************
 * outputs the xml for a motif.
 *
 * Rather than using the motif id embeded in the motif object a passed string is
 * printed. This is because the motif may have been reverse complemented and the
 * strand printed on the front of the id (and the caller may not want that in 
 * the file). As query and target motifs may have conflicting ids a prefix
 * is appended to ensure uniqueness.
 *
 * The indent and tab strings specifiy respectivly, the text string appended to the
 * front of each line, and the text string appended in front of nested tags.
 * 
 * the buffer must be expandable using realloc and the buffer_len must be set
 * to reflect the currently allocated size or zero if the buffer is unallocated
 * and set to null.
 *****************************************************************************/
static void print_xml_motif(FILE *file, int db, MOTIF_T *motif, char *indent, char *tab) {
  int i, j, len;
  char *id, *alt, *url;
  STR_T *b;
  ALPH_T *alph;
  MATRIX_T *freqs;
  double A, C, G, T, log_ev;

  b = str_create(10);
  id = get_motif_id(motif);
  alph = get_motif_alph(motif);
  alt = get_motif_id2(motif);
  len = get_motif_length(motif);
  url = get_motif_url(motif);
  freqs = get_motif_freqs(motif);
  log_ev = get_motif_log_evalue(motif);

  fprintf(file, "%s<motif ", indent);
  fprintf(file, "db=\"%d\" id=\"%s\" ", db, xmlify(id, b, true));
  if (alt && alt[0] != '\0') fprintf(file, "alt=\"%s\" ", xmlify(alt, b, true));
  fprintf(file, "length=\"%d\"", len);
  if (get_motif_nsites(motif) > 0)
    fprintf(file, " nsites=\"%g\"", get_motif_nsites(motif));
  if (log_ev > -HUGE_VAL && log_ev < HUGE_VAL)
    fprintf(file, " evalue=\"%s\"", str_evalue(b, log_ev, 1));
  if (url && url[0] != '\0') {
    fprintf(file, "\n%s%s%surl=\"%s\"", indent, tab, tab, xmlify(url, b, true));
  }
  fprintf(file, ">\n");
  for (i = 0; i < len; ++i) {
    fprintf(file, "%s%s<pos", indent, tab);
    for (j = 0; j < alph_size_core(alph); j++) {
      fprintf(file, " %s=\"%g\"", alph_xml_id(alph, j, b),
          get_matrix_cell(i, j, freqs));
    }
    fprintf(file, "/>\n");
  }
  fprintf(file, "%s</motif>\n", indent);
  str_destroy(b, false);
}

/*****************************************************************************
 * outputs the motifs
 *****************************************************************************/
static RBTREE_T* print_xml_motifs(FILE *out, ARRAYLST_T *dbs, const char* name) {
  RBTREE_T *index_lookup;
  int i, index;
  index_lookup = rbtree_create(compare_tomtom_id, copy_tomtom_id, free, rbtree_intcpy, free);
  fprintf(out, "  <%s>\n", name);
  for (i = 0, index = 0; i < arraylst_size(dbs); i++) {
    MOTIF_DB_T *db;
    RBNODE_T *node;
    db = arraylst_get(i, dbs);
    for (node = rbtree_first(db->matched_motifs); node != NULL; node = rbtree_next(node), index++) {
      MOTIF_T *motif;
      TOMTOM_ID_T id;
      motif = (MOTIF_T*)rbtree_value(node);
      print_xml_motif(out, i, motif, "    ", "  ");
      id.db_id = db->id;
      id.motif_id = get_motif_id(motif);
      if (!rbtree_make(index_lookup, &id, &index)) abort();
    }
  }
  fprintf(out, "  </%s>\n", name);
  return index_lookup;
}

/*****************************************************************************
 * outputs the model section of the xml results
 *
 * the buffer must be expandable using realloc and the buffer_len must be set
 * to reflect the currently allocated size or zero if the buffer is unallocated
 * and set to null.
 *****************************************************************************/
static void print_xml_model(FILE *xml_output,
    int argc, char**argv, 
    ALPH_T *alph, bool rc,
    char *distance_measure, 
    bool sig_type_q, double sig_thresh, 
    char *bg_file, ARRAY_T *bg_freqs_target,
    time_t now) {
  int i;
  bool bg_from_file;
  STR_T *b;
  b = str_create(10);
  fprintf(xml_output, "  <model>\n");
  fprintf(xml_output, "    <command_line>tomtom");
  for (i = 1; i < argc; ++i) {
    char *arg;
    //note that this doesn't correctly handle the case of a " character in a filename
    //in which case the correct output would be an escaped quote or \"
    //but I don't think that really matters since you could guess the original command
    arg = xmlify(argv[i], b,  true);
    if (strchr(argv[i], ' ')) {
      fprintf(xml_output, " &quot;%s&quot;", arg);
    } else {
      fprintf(xml_output, " %s", arg);
    }
  }
  fprintf(xml_output, "</command_line>\n");
  fprintf(xml_output, "    <distance_measure value=\"%s\"/>\n", distance_measure);
  fprintf(xml_output, "    <threshold type=\"%s\">%g</threshold>\n", (sig_type_q ? "qvalue" : "evalue"), sig_thresh);

  alph_print_xml(alph, "alphabet", "    ", "  ", xml_output);
  fprintf(xml_output, "    <strands>%s</strands>\n",  
      (alph_has_complement(alph) ? (rc ? "both" : "forward") : "none"));
  bg_from_file = (bg_file != NULL && strcmp(bg_file, "motif-file") != 0);
  fprintf(xml_output, "    <background from=\"%s\"",
      (bg_from_file ? "file" : "first_target"));
  for (i = 0; i < alph_size_core(alph); i++) {
    fprintf(xml_output, " %s=\"%g\"", alph_xml_id(alph, i, b),
        get_array_item(i, bg_freqs_target));
  }

  if (bg_from_file) {
    fprintf(xml_output, " file=\"%s\"/>\n", xmlify(bg_file, b, true));
  } else {
    fprintf(xml_output, "/>\n");
  }

  fprintf(xml_output, "    <host>%s</host>\n", hostname());
  fprintf(xml_output, "    <when>%s</when>\n", strtok(ctime(&now),"\n"));
  fprintf(xml_output, "  </model>\n");
  str_destroy(b, false);
}

/*****************************************************************************
 * outputs the motif dbs
 *****************************************************************************/
static void print_xml_dbs(FILE *out, ARRAYLST_T *dbs, const char* name) {
  MOTIF_DB_T *db;
  STR_T *b;
  int i;
  b = str_create(10);
  fprintf(out, "  <%s>\n", name);
  for (i = 0; i < arraylst_size(dbs); ++i) {
    db = arraylst_get(i, dbs);
    fprintf(out, "    <db ");
    fprintf(out, "source=\"%s\" ", xmlify(db->source, b, true)); 
    fprintf(out, "name=\"%s\" ", xmlify(db->name, b, true)); 
    fprintf(out, "loaded=\"%d\" ", db->loaded);
    fprintf(out, "excluded=\"%d\" ", db->excluded);
    fprintf(out, "last_mod_date=\"%s\"/>\n", strtok(ctime(&(db->last_mod)),"\n"));
  }
  fprintf(out, "  </%s>\n", name);
  str_destroy(b, false);
}

/*****************************************************************************
 * outputs the xml results
 *****************************************************************************/
static void print_xml_results(FILE *xml_output, 
    int argc, char**argv, 
    ALPH_T *alph,
    bool rc,
    char *distance_measure, 
    bool sig_type_q, double sig_thresh, 
    char *bg_file, ARRAY_T *bg_freqs_target,
    time_t now,
    ARRAYLST_T *query_dbs,
    ARRAYLST_T *target_dbs,
    FILE *matches) {
  RBTREE_T *query_index, *target_index;
  int i, index, last_query;
  double cycles;
  //print dtd
  fprintf(xml_output, "%s", tomtom_dtd);
  //print xml body
  fprintf(xml_output, "<tomtom version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n");
  // print model
  print_xml_model(xml_output, argc, argv, alph, rc, distance_measure, sig_type_q,
      sig_thresh, bg_file, bg_freqs_target, now);
  // output query databases
  print_xml_dbs(xml_output, query_dbs, "query_dbs");
  // output target databases
  print_xml_dbs(xml_output, target_dbs, "target_dbs");
  // output query motifs and track the index we use for each
  query_index = print_xml_motifs(xml_output, query_dbs, "queries");
  // output target motifs and track the index we used for each
  target_index = print_xml_motifs(xml_output, target_dbs, "targets");
  // output matches
  last_query = -1;
  fprintf(xml_output, "  <matches>\n");
  if (fseek(matches, 0L, SEEK_SET) == -1) die("Failed to seek on matches temporary file");
  while (!feof(matches) && !ferror(matches)) {
    int q_db_id, t_db_id, offset, read;
    char q_id[MAX_MOTIF_ID_LENGTH+1];
    char t_id[MAX_MOTIF_ID_LENGTH+1];
    char orient;
    double pvalue, evalue, qvalue;
    read = fscanf(matches, "%d %" STR_VALUE(MAX_MOTIF_ID_LENGTH) "s "
        "%d %" STR_VALUE(MAX_MOTIF_ID_LENGTH) "s %c %d %lf %lf %lf\n",
        &q_db_id, q_id, &t_db_id, t_id, &orient, &offset, &pvalue,
        &evalue, &qvalue);
    if (read == 9) {
      TOMTOM_ID_T id;
      int q_index, t_index;
      // lookup the query index
      id.db_id = q_db_id;
      id.motif_id = q_id;
      q_index = *((int*)rbtree_get(query_index, &id));
      if (q_index != last_query) {
        if (last_query != -1) fprintf(xml_output, "    </query>\n");
        fprintf(xml_output, "    <query idx=\"%d\">\n", q_index);
        last_query = q_index;
      }
      // lookup the target index
      id.db_id = t_db_id;
      id.motif_id = t_id;
      t_index = *((int*)rbtree_get(target_index, &id));
      // write out the match
      fprintf(xml_output, "      <target idx=\"%d\"", t_index);
      if (rc) fprintf(xml_output, " rc=\"%s\"", (orient == '-' ? "y" : "n"));
      fprintf(xml_output, " off=\"%d\" pv=\"%.2e\" ev=\"%.2e\" qv=\"%.2e\"/>\n",
          offset, pvalue, evalue, qvalue);
    }
  }
  if (last_query != -1) fprintf(xml_output, "    </query>\n");
  fprintf(xml_output, "  </matches>\n");

  cycles = myclock(); //measures number of cycles since first call to myclock (unfortunately prone to overflow for large run times).

  fprintf(xml_output, "  <runtime cycles=\"%.0f\" seconds=\"%.3f\"/>\n", cycles, cycles / CLOCKS_PER_SEC);
  fprintf(xml_output, "</tomtom>\n");
  // cleanup
  rbtree_destroy(query_index);
  rbtree_destroy(target_index);
}

void pmf_to_cdf (ARRAY_T* pmf, ARRAY_T* cdf) {
  // cum sum of pmf is cdf (no checking of dimensions!)
  ATYPE last=0;
  int i;
  for (i = 0; i < get_array_length(pmf); i++) {
    last += get_array_item(i, pmf);
    set_array_item(i, last, cdf);
  }
}

void pmf_to_pv (ARRAY_T* pmf, ARRAY_T* pv) {
  // like pmf_to_cdf for pv (so cumsum of reversed pmf)
  ATYPE last=0;
  int i;
  for (i = get_array_length(pmf)-1; i >= 0; i--) {
    last += get_array_item(i, pmf);
    set_array_item(i, last, pv);
  }
}

/*****************************************************************************
 * Calculates the p-values of the shifted scores by finding the distribution of
 * the maximal shifted alignments score for each possible target motif length UK
 *****************************************************************************/
void calc_shifted_scores_pvalues_per_query(bool rc, int query_len, MLEN_T *tlen,
    int target_count, MATRIX_T* reference_matrix, MATRIX_T* pmf_matrix,
    double scores_offset, double scores_scale,
    TOMTOM_MATCH_T **match_list) {
  int len, start_col, end_col, j, i, k, pmf_index, offset_correction, iorient, norient;
  int pmf_len, ext_pmf_len, num_pv_lookup_array;
  MATRIX_T *pv_of_max_algn_matrix, *extended_pmf_matrix;
  ARRAY_T *pmf_new, *pmf_old, *pmf_current, *cdf_old, *cdf_current;

  norient = (rc ? 2 : 1);

  pmf_len = get_num_cols(pmf_matrix);
  // pmf size if offsets are peeled off, note that for an positive quantile 
  // shifted scores offset is always >= 0
  ext_pmf_len = (int) ceil(pmf_len - (query_len-1) * scores_offset * scores_scale);

  // This can be made smaller if the attained scores are scanned first
  pv_of_max_algn_matrix = allocate_matrix(tlen->unique, ext_pmf_len);

  num_pv_lookup_array = get_num_rows(pmf_matrix);

  // for embeding the pmf with offset peeled
  extended_pmf_matrix = allocate_matrix(num_pv_lookup_array, ext_pmf_len);
  init_matrix(0, extended_pmf_matrix);

  for (start_col = 0; start_col < query_len; start_col++) {
    // properly embed each alignment pmf in the extended array
    for (end_col = start_col; end_col < query_len; end_col++) {
      pmf_index = get_matrix_cell(start_col, end_col, reference_matrix);
      offset_correction = (int) rint(scores_offset * (query_len - 1 - (end_col-start_col)) * scores_scale);

      for (j = 0; j < pmf_len; j++) {
        set_matrix_cell(
            pmf_index, 
            j-offset_correction, 
            get_matrix_cell(pmf_index, j, pmf_matrix), 
            extended_pmf_matrix
            );
      }
    }
  }

  pmf_new = allocate_array(ext_pmf_len);
  pmf_old = allocate_array(ext_pmf_len);
  cdf_old = allocate_array(ext_pmf_len);
  cdf_current = allocate_array(ext_pmf_len);

  for (j = 0; j < tlen->unique; j++) {
    int nalgns;
    // compute the probability mass function (pmf) of the best alignment
    // shifted score for each target motif width
    init_array(0, pmf_new);
    set_array_item(0, 1.0, pmf_new);
    len = tlen->lengths[j];
    nalgns = query_len + len - 1;

    for (i = 0; i < nalgns; i++) {
      if (i < query_len) {
        // find the start and end columns of this alignment
        start_col = query_len - i - 1;
        end_col = MIN(query_len, start_col+len) - 1;
      } else {
        start_col = 0;
        end_col = MIN(query_len, len + query_len - i - 1) - 1;
      }

      pmf_index = get_matrix_cell(start_col, end_col, reference_matrix);
      pmf_current = get_matrix_row(pmf_index, extended_pmf_matrix);
      pmf_to_cdf(pmf_current, cdf_current);

      for (iorient = 0; iorient < norient; iorient++) {
        // remember that the target is tested in both orientations
        SWAP(ARRAY_T*, pmf_new, pmf_old)
          pmf_to_cdf(pmf_old, cdf_old);

        // next is the heart of the calculation of the pmf 
        // (can be made faster if observed scores are scanned first)
        set_array_item(0, get_array_item(0, pmf_current) * get_array_item(0, pmf_old), pmf_new);
        for (k = 1; k < ext_pmf_len; k++) {
          ATYPE new;
          new = get_array_item(k, pmf_current) * get_array_item(k-1, cdf_old) 
            + get_array_item(k, pmf_old) * get_array_item(k-1, cdf_current)
            + get_array_item(k, pmf_current) * get_array_item(k, pmf_old);
          set_array_item(k, new, pmf_new);
        }  // loop on possible values (k)
      }  // loop on two orientations (iorient)
    }  // loop on alignments (i)

    pmf_to_pv(pmf_new, pmf_old);  // pmf_old is just for storage, it's really pv

    set_matrix_row(j, pmf_old, pv_of_max_algn_matrix);

  }  // loop on target lengths (j)

  // Traverse the list of matches and replace the pvalue which is really 
  // the complete score with its pvalue
  TOMTOM_MATCH_T *match;
  int score;
  for (i = 0; i < target_count; i++) {
    match = match_list[i];
    score = rint(match->pvalue - query_len * scores_offset * scores_scale);
    len = get_motif_length(match->target);
    // copy the double into a float to replicate the old rounding behaviour
    float truncated_pvalue = get_matrix_cell(tlen->lookup[len], score, pv_of_max_algn_matrix);
    match->pvalue =  truncated_pvalue;
    match->evalue = truncated_pvalue * target_count / norient;
  }

  free_matrix(pv_of_max_algn_matrix);
  free_matrix(extended_pmf_matrix);
  free_array(pmf_new);
  free_array(pmf_old);
  free_array(cdf_old);
  free_array(cdf_current);

}

/*****************************************************************************
 * Calculate a list of unique motif lengths as well as min, max and total.
 * These help save time later.
 *****************************************************************************/
static MLEN_T* create_mlen(ARRAYLST_T *motifs) {
  int i, nmotifs, len;
  MLEN_T *mlen;
  mlen = mm_malloc(sizeof(MLEN_T));
  mlen->min = 0;
  mlen->max = 0;
  mlen->total = 0;
  nmotifs = arraylst_size(motifs);
  for (i = 0; i < nmotifs; i++) {
    len = get_motif_length((MOTIF_T*)arraylst_get(i, motifs));
    // determine minimum and maximum length
    if (mlen->min == 0 || len < mlen->min) {
      mlen->min = len;
    }
    if (mlen->max == 0 || len > mlen->max) {
      mlen->max = len;
    }
    // determine total length
    mlen->total += len;
  }
  assert(mlen->min <= mlen->max);
  // count number of unique lengths
  mlen->unique = 0;
  mlen->lookup = mm_calloc(mlen->max + 1, sizeof(int));
  for (i = 0; i < nmotifs; i++) {
    len = get_motif_length((MOTIF_T*)arraylst_get(i, motifs));
    assert(len >= mlen->min);
    assert(len <= mlen->max);
    if (mlen->lookup[len] == 0) {
      mlen->unique++;
      mlen->lookup[len] = 1;
    }
  }
  // create list of possible lengths
  mlen->lengths = mm_malloc(mlen->unique * sizeof(int));
  for (i = 0, len = mlen->min; len <= mlen->max; len++) {
    if (mlen->lookup[len] != 0) {
      mlen->lengths[i] = len;
      mlen->lookup[len] = i;
      i++;
    }
  }
  assert(i == mlen->unique);
  return mlen;
}

/*****************************************************************************
 * Destroy the list of motif lengths.
 *****************************************************************************/
static void destroy_mlen(MLEN_T *mlen) {
  free(mlen->lookup);
  free(mlen->lengths);
  free(mlen);
}

/*****************************************************************************
 * Load motifs from a list of databases.
 * Optionally filter if ids or idxs are specified.
 * Store/check alphabet with alph.
 * Optionally store first motif db background in bg.
 *****************************************************************************/
static ARRAYLST_T* load_motifs(ARRAYLST_T *dbs, RBTREE_T *ids, RBTREE_T *idxs, bool allow_zeros, double pseudocount, bool xalph, ALPH_T **alph, ARRAY_T **bg) {
  ARRAYLST_T *motifs;
  MOTIF_DB_T *db;
  MREAD_T *mread;
  MOTIF_T *motif;
  RBTREE_T *dups;
  ALPH_T *db_alph;
  int i, idx;
  bool printed_id_warn, printed_zeros_warn;
  printed_id_warn = false;
  printed_zeros_warn = false;
  // create motif array
  motifs = arraylst_create();
  for (i = 0; i < arraylst_size(dbs); i++) {
    // get the db
    db = (MOTIF_DB_T*)arraylst_get(i, dbs);
    // set the start index of the list
    db->list_index = arraylst_size(motifs);
    // now start reading
    mread = mread_create(db->source, OPEN_MFILE);
    mread_set_pseudocount(mread, pseudocount);
    if (*bg != NULL) mread_set_background(mread, *bg);
    // check the alphabet
    db_alph = mread_get_alphabet(mread);
    if (*alph != NULL) {
      // check if the target motif has an alphabet that matches the query motif
      if (mread_has_motif(mread) && !alph_equal(*alph, db_alph)) {
        if (xalph) {
          // the alphabet does not match but we have been told to convert it so
          // we must check that such a thing is possible
          switch(alph_core_subset(db_alph, *alph)) {
            case 0:
              die("The motifs in \"%s\" have the %s alphabet which is not "
                  "a subset of the %s alphabet.", db->name, 
                  alph_name(db_alph), alph_name(*alph));
              break;
            case -1:
              fprintf(stderr, "Warning: the alphabet expansion from %s to %s"
                  " requires changing complementing rules.\n", 
                  alph_name(db_alph), alph_name(*alph));
              break;
          }
          mread_set_conversion(mread, *alph, *bg);
        } else {
          die("The motifs in \"%s\" have the %s alphabet which is not "
              "the same as the expected %s alphabet.", db->name, 
              alph_name(db_alph), alph_name(*alph));
        }
      }
    } else {
      // store the alphabet of the query motif
      *alph = alph_hold(db_alph);
    }
    // store the background (should happen on query motif)
    if (*bg == NULL) *bg = mread_get_background(mread);
    // create lookup for duplicate checking within the database
    dups = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
    printed_id_warn = false;
    // now read each motif
    while (mread_has_motif(mread)) {
      // read one motif
      motif = mread_next_motif(mread);
      db->loaded += 1;
      // check that this motif has a unique identifier
      if (!rbtree_make(dups, get_motif_id(motif), NULL)) {
        // duplicate motif ID!
        if (printed_zeros_warn) {
          fprintf(stderr, "\n");
          printed_zeros_warn = false;
        }
        if (!printed_id_warn) {
          fprintf(stderr, "Warning, the following motifs in \"%s\" were "
              "removed because they had a duplicate name:", db->name);
          printed_id_warn = true;
        }
        fprintf(stderr, " %s", get_motif_id(motif));
        destroy_motif(motif);
        db->excluded += 1;
        continue;
      }
      // check that this motif has not been filtered
      if ((ids && rbtree_size(ids) > 0) || (idxs && rbtree_size(idxs) > 0)) {
        idx = get_motif_idx(motif);
        if (ids && rbtree_find(ids, get_motif_id(motif))) {
          // keep motif as the ID is selected
        } else if (idxs && rbtree_find(idxs, &idx)) {
          // keep motif as the index is selected
        } else {
          // not in the selected motifs
          destroy_motif(motif);
          db->excluded += 1;
          continue;
        }
      }
      if (!allow_zeros && has_motif_zeros(motif)) {
        if (printed_id_warn) {
          fprintf(stderr, "\n");
          printed_id_warn = false;
        }
        if (!printed_zeros_warn) {
          fprintf(stderr, "Warning, the following motifs in \"%s\" were "
              "removed because they contained probabilities of zero in some "
              "columns (add pseudocounts and rerun!):", db->name);
          printed_zeros_warn = true;
        }
        fprintf(stderr, " %s", get_motif_id(motif));
        destroy_motif(motif);
        db->excluded += 1;
        continue;
      }
      // add the motif to the list
      arraylst_add(motif, motifs);
    }
    if (printed_zeros_warn || printed_id_warn) {
      fprintf(stderr, "\n");
    }
    db->list_entries = db->loaded - db->excluded;
    if (db->list_entries == 0) {
      fprintf(stderr, "Warning: the database \"%s\" did not contain ANY "
          "loadable motifs!\n", db->name);
    }
    //destroy the set of ids for this database
    rbtree_destroy(dups);
    // destroy the motif file reader for this database
    mread_destroy(mread);
  }
  return motifs;
}



/*****************************************************************************
 * Tomtom program
 *****************************************************************************/
VERBOSE_T verbosity = NORMAL_VERBOSE;
#ifdef MAIN

/*****************************************************************************
 * Print a usage message with an optional reason for failure.
 *****************************************************************************/
static void usage(char *format, ...) {
  va_list argp;
  char *usage = 
    "\n"
    "Usage:\n"
    "  tomtom [options] <query file> <target file>+\n"
    "Options:\n"
    "  -o <output dir>  Name of directory for output files;\n"
    "                    will not replace existing directory\n"
    "  -oc <output dir> Name of directory for output files;\n"
    "                    will replace existing directory\n"
    "  -xalph           Convert the alphabet of the target motif databases\n"
    "                    to the alphabet of the query motif database\n"
    "                    assuming the core symbols of the target motif\n"
    "                    alphabet are a subset; default: reject differences\n"
    "  -bfile <background file>\n"
    "                   Name of background file;\n"
    "                    default: use the background from the query\n"
    "                    motif database\n"
    "  -motif-pseudo <pseudo count>\n"
    "                   Apply the pseudocount to the query and target motifs;\n"
    "                    default: apply a pseudocount of 0.1\n"
    "  -m <id>          Use only query motifs with a specified id;\n"
    "                    may be repeated\n"
    "  -mi <index>      Use only query motifs with a specifed index;\n"
    "                    may be repeated\n"
    "  -thresh <float>  Significance threshold; default: 0.5\n"
    "  -evalue          Use E-value threshold; default: q-value\n"
    "  -dist allr|ed|kullback|pearson|sandelin|blic1|blic5|llr1|llr5\n"
    "                   Distance metric for scoring alignments;\n"
    "                    default: pearson\n"
    "  -internal        Only allow internal alignments;\n"
    "                    default: allow overhangs\n"
    "  -min-overlap <int>\n"
    "                   Minimum overlap between query and target;\n"
    "                    default: 1\n"
    "  -incomplete-scores\n"
    "                   Ignore unaligned columns in computing scores\n"
    "                    default: use complete set of columns\n"
    "  -text            Output in text format (default is HTML)\n"
    "  -png             Create PNG logos; default: don't create PNG logos\n"
    "  -eps             Create EPS logos; default: don't create EPS logos\n"
    "  -no-ssc          Don't apply small-sample correction to logos;\n"
    "                    default: use small-sample correction\n"
    "  -verbosity [1|2|3|4]\n"
    "                   Set the verbosity of the program; default: 2 (normal)\n"
    "  -version         Print the version and exit\n"
    ;
  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fputs(usage, stderr);
    fflush(stderr);
  } else {
    fputs(usage, stderr);
  }
  exit(1);
}

// An easy way of assigning a number to each option.
// Be careful that there are not more than 62 options as otherwise the value
// of '?' will be used.
enum Opts {OPT_O, OPT_OC, OPT_BFILE, OPT_M, OPT_MI, OPT_TEXT, OPT_THRESH, 
  OPT_QTHRESH, OPT_EVALUE, OPT_DIST, OPT_INTERNAL, OPT_OVERLAP, OPT_PSEUDO,
  OPT_COLSCORES, OPT_PNG, OPT_EPS, OPT_NO_SSC,
  OPT_IC_SCORES, OPT_NO_RC, OPT_XALPH, OPT_VERBOSITY, OPT_VERSION};

/*****************************************************************************
 * process command line arguments and store the results in config.
 *****************************************************************************/
static void process_arguments(int argc, char** argv, TOMTOM_T *config) {
  char *endptr;
  long num;
  struct option tomtom_options[] = {
    {"o", required_argument, NULL, OPT_O},
    {"oc", required_argument, NULL, OPT_OC},
    {"bfile", required_argument, NULL, OPT_BFILE},
    {"m", required_argument, NULL, OPT_M},
    {"mi", required_argument, NULL, OPT_MI},
    {"text", no_argument, NULL, OPT_TEXT},
    {"thresh", required_argument, NULL, OPT_THRESH},
    {"q-thresh", required_argument, NULL, OPT_QTHRESH},
    {"evalue", no_argument, NULL, OPT_EVALUE},
    {"dist", required_argument, NULL, OPT_DIST},
    {"internal", no_argument, NULL, OPT_INTERNAL},
    {"min-overlap", required_argument, NULL, OPT_OVERLAP},
    {"motif-pseudo", required_argument, NULL, OPT_PSEUDO},
    {"column-scores", required_argument, NULL, OPT_COLSCORES},
    {"png", no_argument, NULL, OPT_PNG},
    {"eps", no_argument, NULL, OPT_EPS},
    {"no-ssc", no_argument, NULL, OPT_NO_SSC},
    {"incomplete-scores", no_argument, NULL, OPT_IC_SCORES},
    {"norc", no_argument, NULL, OPT_NO_RC},
    {"xalph", no_argument, NULL, OPT_XALPH},
    {"verbosity", required_argument, NULL, OPT_VERBOSITY}, //does nothing
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0} //boundary indicator
  };
  // set defaults
  config->now = time(NULL);
  config->q_files = arraylst_create(1);
  config->t_files = arraylst_create(1);
  config->clobber = true;
  config->outdir = "tomtom_out";
  config->bg_file = NULL;
  config->ids = rbtree_create(rbtree_strcmp, NULL, NULL, NULL, NULL);
  config->idxs = rbtree_create(rbtree_longcmp, rbtree_longcpy, free, NULL, NULL);
  config->text_only = false;
  config->html = true;
  config->sig_thresh = 0.5;
  config->sig_type_q = true;
  config->dist_name = "ed";
  config->dist_func = ed_scores;
  config->dist_dna_only = false;
  config->dist_allow_zeros = true;
  config->internal = false;
  config->min_overlap = 1;
  config->pseudo = 0.1;
  config->cs_file = NULL;
  config->rc = true;
  config->xalph = false;
  config->png = false;
  config->eps = false;
  config->ssc = true;
  config->complete_scores = true;

  // process arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", tomtom_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OPT_O:
      case OPT_OC:
        config->clobber = (opt == OPT_OC);
        config->outdir = optarg;
        break;
      case OPT_BFILE:
        config->bg_file = optarg;
        break;
      case OPT_M:
        rbtree_put(config->ids, optarg, NULL); 
        break;
      case OPT_MI:
        errno = 0;
        num = strtol(optarg, &endptr, 10);
        if ((errno != 0 && num == 0) || *endptr != '\0') {
          usage("Option -mi \"%s\" was not a number", optarg);
        } else if (num <= 0) {
          usage("Option -mi \"%s\" can not be zero or negative", optarg);
        }
        rbtree_put(config->idxs, &num, NULL); 
        break;
      case OPT_TEXT:
        config->text_only = true;
        break;
      case OPT_QTHRESH:
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(stderr, "Warning: -thresh should be used in preference to "
              "-q-thresh as -q-thresh is deprecated.\n");
        }
        // fall through
      case OPT_THRESH:
        config->sig_thresh = strtod(optarg, &endptr);
        if (*endptr != '\0') {
          usage("Option -thresh \"%s\" was not a number", optarg);
        } else if (config->sig_thresh < 0) {
          usage("Option -thresh \"%s\" can not be negative", optarg);
        }
        break;
      case OPT_EVALUE:
        config->sig_type_q = false;
        break;
      case OPT_DIST:
        config->dist_name = optarg;
        config->dist_dna_only = false;
        if (strcmp(optarg, "allr") == 0) {
          config->dist_func = allr_scores;
          config->dist_allow_zeros = false;
        } else if (strcmp(optarg, "ed") == 0) {
          config->dist_func = ed_scores;
        } else if (strcmp(optarg, "kullback") == 0) {
          config->dist_func = kullback_scores;
        } else if (strcmp(optarg, "pearson") == 0) {
          config->dist_func = pearson_scores;
        } else if (strcmp(optarg, "sandelin") == 0) {
          config->dist_func = sandelin_scores;
        } else if (strcmp(optarg, "blic1") == 0) {
          config->dist_func = blic1_scores;
          config->dist_dna_only = true;
        } else if (strcmp(optarg, "blic5") == 0) {
          config->dist_func = blic5_scores;
          config->dist_dna_only = true;
        } else if (strcmp(optarg, "llr1") == 0) {
          config->dist_func = llr1_scores;
          config->dist_dna_only = true;
        } else if (strcmp(optarg, "llr5") == 0) {
          config->dist_func = llr5_scores;
          config->dist_dna_only = true;
        } else {
          usage("Option -dist \"%s\" was not an allowed value", optarg);
        }
        break;
      case OPT_INTERNAL:
        config->internal = true;
        break;
      case OPT_OVERLAP:
        num = strtol(optarg, &endptr, 10);
        if ((errno != 0 && num == 0) || *endptr != '\0') {
          usage("Option -min-overlap \"%s\" was not a number", optarg);
        } else if (num < 0) {
          usage("Option -min-overlap \"%s\" can not be negative", optarg);
        }
        config->min_overlap = num;
        break;
      case OPT_PSEUDO:
        config->pseudo = strtod(optarg, &endptr);
        if (*endptr != '\0') {
          usage("Option -motif-pseudo \"%s\" was not a number", optarg);
        } else if (config->pseudo < 0) {
          usage("Option -motif-pseudo \"%s\" can not be negative", optarg);
        }
        break;
      case OPT_COLSCORES: // DEV only option
        config->cs_file = optarg;
        break;
      case OPT_PNG:
        config->png = true;
        break;
      case OPT_EPS:
        config->eps = true;
        break;
      case OPT_NO_SSC:
        config->ssc = false;
        break;
      case OPT_IC_SCORES:
        config->complete_scores = false;
        break;
      case OPT_NO_RC:
        config->rc = false;
        break;
      case OPT_XALPH:
        config->xalph = true;
        break;
      case OPT_VERBOSITY:
        num = strtol(optarg, &endptr, 10);
        if ((errno != 0 && num == 0) || *endptr != '\0') {
          usage("Option -verbosity \"%s\" was not a number", optarg);
        } else if (num < QUIET_VERBOSE || num > DUMP_VERBOSE) {
          usage("Option -verbosity \"%s\" was not an allowed value", optarg);
        }
        verbosity = num;
        break;
      case OPT_VERSION:
        fprintf(stdout, VERSION "\n");
        exit(EXIT_SUCCESS);
        break;
      case '?':
        usage(NULL);
    }
  }
  if (config->sig_type_q && config->sig_thresh > 1) {
    usage("Option -thresh (%g) must be less than 1 unless -evalue is specified",
        config->sig_thresh);
  }
  if (optind < argc) {
    arraylst_add(argv[optind], config->q_files);
    optind++;
  }
  while (optind < argc) {
    arraylst_add(argv[optind], config->t_files);
    optind++;
  }
  if (arraylst_size(config->q_files) == 0) {
    usage("No query motif database supplied");
  }
  if (arraylst_size(config->t_files) == 0) {
    usage("No target motif database supplied");
  }
}

static void cleanup_arguments(TOMTOM_T *config) {
  arraylst_destroy(NULL, config->q_files);
  arraylst_destroy(NULL, config->t_files);
  rbtree_destroy(config->ids);
  rbtree_destroy(config->idxs);
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
int main(int argc, char *argv[]) {
  TOMTOM_T opts;
  ARRAYLST_T *query_dbs, *target_dbs, *query_motifs, *target_motifs;
  MLEN_T *t_mlen;
  FILE *text_output; // where to write legacy output
  FILE *match_tmp = NULL;
  ALPH_T *alph;
  ARRAY_T *bg;
  bool stdin_used;
  int i, j, query_count, target_count, targets_len, qdb_i, qdb_remain;
  TOMTOM_MATCH_T *match, **match_list;
  RBTREE_T *seen;
  char *query_consensus, *target_consensus;
  bool output_error = false;
  // **********************************************
  // record cpu time
  // **********************************************
  (void) myclock();
  // **********************************************
  // process options
  // **********************************************
  process_arguments(argc, argv, &opts);
  // **********************************************
  // load background from a file
  // **********************************************
  bg = NULL;
  if (opts.bg_file != NULL && strcmp(opts.bg_file, "--motif--") != 0 && strcmp(opts.bg_file, "motif-file") != 0) {
    int order = 0;
    char *syms = NULL;
    bg = load_markov_model_without_alph(opts.bg_file, &order, &syms); // load 0-order bg
  }
  // **********************************************
  // check databases
  // **********************************************
  stdin_used = false;
  query_dbs = create_motif_dbs(opts.q_files, "query", opts.now, &stdin_used);
  target_dbs = create_motif_dbs(opts.t_files, "target", opts.now, &stdin_used);
  // **********************************************
  // load motifs
  // **********************************************
  alph = (opts.dist_dna_only ? alph_dna() : NULL);
  query_motifs = load_motifs(query_dbs, opts.ids, opts.idxs, opts.dist_allow_zeros, opts.pseudo, false, &alph, &bg);
  target_motifs = load_motifs(target_dbs, NULL, NULL, opts.dist_allow_zeros, opts.pseudo, opts.xalph, &alph, &bg);
  // check that the target motif count is large enough
  if (arraylst_size(target_motifs) < MIN_TARGET_DATABASE_SIZE) {
    fprintf(stderr,
        "Warning: Target database size too small (%d). "
        "Provide at least %d motifs for accurate p-value computation\n",
        arraylst_size(target_motifs), MIN_TARGET_DATABASE_SIZE);
  }
  // **********************************************
  // Reverse complement the target motifs
  // **********************************************
  if (!alph_has_complement(alph)) opts.rc = false; // check for alphabet support
  if (opts.rc) {
    add_reverse_complements(target_motifs); // add RC motifs at odd indexes
    for (i = 0; i < arraylst_size(target_dbs); ++i) { // double db sizes
      MOTIF_DB_T *db = (MOTIF_DB_T*)arraylst_get(i, target_dbs);
      db->list_index *= 2;
      db->list_entries *= 2;
    }
  }
  // *********************************************
  // Calculate motif length stats
  // *********************************************
  t_mlen = create_mlen(target_motifs);
  // **********************************************
  // create output files
  // **********************************************
  if (opts.text_only == true) {
    // Legacy: plain text output to standard out.
    text_output = stdout;
  } else {
    char *path;
    if (create_output_directory(opts.outdir, opts.clobber, (verbosity >= NORMAL_VERBOSE))) {
      // Failed to create output directory.
      die("Unable to create output directory %s.\n", opts.outdir);
    }
    // Create the name of the output files (text, XML, HTML) and open them
    path = make_path_to_file(opts.outdir, TXT_FILENAME);
    if ((text_output = fopen(path, "w")) == NULL) {
      die("Unable to open \"%s\" for writing.\n", path);
    }
    free(path);
    if ((match_tmp = tmpfile()) == NULL) {
      die("Unable to create temporary file for storing results.\n");
    }
  }
  // ********************************************
  // Start text output
  // ********************************************
  fprintf(text_output, "#Query ID\tTarget ID\tOptimal offset\tp-value\tE-value"
      "\tq-value\tOverlap\tQuery consensus\tTarget consensus\tOrientation\n");
  // ********************************************
  // Setup for the query loop
  // ********************************************
  query_count = arraylst_size(query_motifs);
  target_count = arraylst_size(target_motifs);
  targets_len = t_mlen->total;
  seen = rbtree_create(compare_tomtom_id, copy_tomtom_id, free, NULL, NULL);
  match_list = mm_malloc(target_count * sizeof(TOMTOM_MATCH_T*));
  qdb_i = -1; qdb_remain = 0;
  MOTIF_DB_T *query_db;
  query_db = NULL;
  // ********************************************
  // Loop over queries
  // ********************************************
  for (i = 0; i < query_count; i++, qdb_remain--) {
    double pssm_offset, pssm_scale;
    int query_len, target_len, col_offset, tdb_remain, tdb_i;
    int num_pv_lookup_array, num_array_lookups;
    MOTIF_T *query_motif, *target_motif, *original_motif;
    MATRIX_T* all_columnwise_scores, *reference_matrix, *pv_lookup_matrix, *pmf_matrix;
    // report progress
    if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Processing query %d out of %d \n", i + 1, query_count);
    // get the query database
    while (qdb_remain <= 0) {
      query_db = arraylst_get(++qdb_i, query_dbs);
      qdb_remain = query_db->list_entries;
    }
    // get the query motif
    query_motif = (MOTIF_T*)arraylst_get(i, query_motifs); 
    query_len = get_motif_length(query_motif);

    // Initialize an array to store the columnwise scores.
    all_columnwise_scores = allocate_matrix(query_len, targets_len);
    // ********************************************
    // Apply the distance function to each target
    // ********************************************
    col_offset = 0;
    for (j = 0; j < target_count; j++) {
      target_motif = (MOTIF_T*)arraylst_get(j, target_motifs);
      // Scores of query motif columns (the rows) vs all target columns.
      // Each row of the frequency array contains the scores of a single query
      // column vs. all possible target columns.
      opts.dist_func(query_motif, target_motif, bg, col_offset, &all_columnwise_scores);
      col_offset += get_motif_length(target_motif);
    }
    // If requested, store the columnwise scores in an external file.
    if (opts.cs_file != NULL) print_columnwise_scores(opts.cs_file, all_columnwise_scores);
    // UK: shift all column similarity scores
    if (opts.complete_scores) shift_all_pairwise_scores(all_columnwise_scores,
        query_len, targets_len, SHIFT_QUANTILE);
    // scale the columnwise scores
    scale_score_matrix(all_columnwise_scores, query_len, targets_len, NULL, 1.0,
        BINS, &pssm_offset, &pssm_scale);

    // ********************************************
    // Compute the pvalue distributions
    // ********************************************
    num_pv_lookup_array = (query_len * (query_len + 1))/2; // w * (w+1) / 2
    num_array_lookups = (BINS * query_len) + 1;
    // setup reference matrix
    reference_matrix = allocate_matrix(query_len, query_len);
    init_matrix(-1, reference_matrix);
    // setup pv lookup matrix
    pv_lookup_matrix = allocate_matrix(num_pv_lookup_array, num_array_lookups);
    init_matrix(0, pv_lookup_matrix);
    // setup pmf matrix
    pmf_matrix = allocate_matrix(num_pv_lookup_array, num_array_lookups);   // UK
    init_matrix(0, pmf_matrix);
    // calculate pvalues
    get_pv_lookup_new(all_columnwise_scores, query_len, targets_len, BINS,
        reference_matrix, pv_lookup_matrix, pmf_matrix);

    // ********************************************
    // Extract the scaled scores for each target
    // and compute the optimal offset for the
    // best p-value
    // ********************************************
    tdb_i = -1; tdb_remain = 0; col_offset = 0;
    MOTIF_DB_T *target_db;
    target_db = NULL;
    for (j = 0; j < target_count; ++j, --tdb_remain) {
      int optimal_offset, optimal_overlap, total_configurations;
      double optimal_pvalue, motif_pvalue, motif_evalue;
      // get the target database
      while (tdb_remain <= 0) {
        target_db = arraylst_get(++tdb_i, target_dbs);
        tdb_remain = target_db->list_entries;
      }
      // get the target motif
      target_motif = (MOTIF_T *) arraylst_get(j, target_motifs);
      // if we check reverse complements then get the normal motif
      if (!opts.rc || j % 2 == 0) {
        original_motif = target_motif;
      } else {
        original_motif = (MOTIF_T*)arraylst_get(j-1, target_motifs);
      }
      // get the target len
      target_len = get_motif_length(target_motif);

      // Compute the score and offset wrt the minimum pvalue CORRECTED.
      compare_motifs(all_columnwise_scores, query_len, target_len, col_offset,
          reference_matrix, pv_lookup_matrix,
          opts.internal, opts.min_overlap, opts.complete_scores, 
          pssm_scale, pssm_offset, &total_configurations,
          &optimal_offset, &optimal_pvalue, &optimal_overlap);

      if (opts.complete_scores) {
        motif_pvalue = optimal_pvalue; // UK: this is in fact the best shifted score
        motif_evalue = 0.0;
      } else {
        // Calculates p-values and e-values of unshfited scores.
        motif_pvalue = EV(optimal_pvalue, (opts.rc ? total_configurations * 2 : total_configurations));
        // we've already corrected for 2 strands in the p-value so don't correct twice!
        motif_evalue = motif_pvalue * (opts.rc ? target_count / 2 : target_count);
      }
      // Store this match.
      match_list[j] = new_tomtom_match(target_db, target_motif, original_motif,
          optimal_offset, optimal_overlap, motif_pvalue, motif_evalue);
      // update the column offset so we look at the right scores
      col_offset += get_motif_length(target_motif);
    }
    // ********************************************
    // Post process the matches to calculate the
    // p-values (in complete scoring mode) and q-values
    // ********************************************
    if (opts.complete_scores && target_count > 0) {
      calc_shifted_scores_pvalues_per_query(opts.rc, query_len, t_mlen, 
          target_count, reference_matrix, pmf_matrix, pssm_offset, 
          pssm_scale, match_list);
    }
    // Sort the match list
    qsort(match_list, target_count, sizeof(TOMTOM_MATCH_T*), (void*)tomtom_match_compare);
    // Compute q-values.
    convert_tomtom_p_to_q(target_count, match_list);
    // cleanup
    free_matrix(reference_matrix);
    free_matrix(pv_lookup_matrix);
    free_matrix(pmf_matrix);
    free_matrix(all_columnwise_scores);
    // ********************************************
    // Write out the results for a query motif
    // ********************************************
    query_consensus = get_cons(query_motif);
    rbtree_clear(seen); // clear seen motifs
    for (j = 0; j < target_count; j++) {
      TOMTOM_ID_T id;
      match = match_list[j];
      // copy the doubles into floats to try to replicate the old rounding behaviour
      float match_pvalue, match_evalue, match_qvalue;
      match_pvalue = match->pvalue;
      match_evalue = match->evalue;
      match_qvalue = match->qvalue;
      // only output significant matches
      if ((opts.sig_type_q ? match->qvalue : match->evalue) > opts.sig_thresh) break;
      // skip matches where we've already seen the opposite strand
      id.db_id = match->db->id;
      id.motif_id = get_motif_id(match->target);
      if (!rbtree_make(seen, &id, NULL)) continue;
      // record that the motif was used so we can output it later
      rbtree_put(match->db->matched_motifs, get_motif_id(match->original), match->original);
      // generate a consensus for the target motif
      target_consensus = get_cons(match->target);
      // Print text output. 
      // Note that this doesn't record the target database of the match so
      // it could be ambiguous.
      fprintf(text_output,
          "%s\t%s\t%d\t%g\t%g\t%g\t%d\t%s\t%s\t%c\n",
          get_motif_id(query_motif),
          get_motif_id(match->target),
          match->offset,
          match_pvalue,
          match_evalue,
          match_qvalue,
          match->overlap,
          query_consensus,
          target_consensus,
          get_motif_strand(match->target));
      if (opts.png || opts.eps) {
	      create_logo(
            match->db, 
            match->target, 
            query_db, 
            query_motif, 
            match->offset, 
            opts.png, 
            opts.eps, 
            opts.ssc, 
            opts.outdir
          );
      }
      // clean up
      free(target_consensus);
      // store match for XML output.
      if (!opts.text_only) {
        fprintf(match_tmp, "%d %s %d %s %c %d %g %g %g\n", 
            query_db->id, get_motif_id(query_motif),
            match->db->id, get_motif_id(match->target),
            get_motif_strand(match->target), match->offset,
            match->pvalue, match->evalue, match->qvalue);
      }
    }
    // clean up
    free(query_consensus);
    // record that the query motif was used
    rbtree_put(query_db->matched_motifs, get_motif_id(query_motif), query_motif);
    // clean up match list items but keep the list for reuse
    for (j = 0; j < target_count; j++) {
      free(match_list[j]);
      match_list[j] = NULL;
    }
  }
  // clean up
  rbtree_destroy(seen);
  free(match_list);

  // ********************************************
  // Write out the XML results
  // ********************************************
  if (!opts.text_only) {
    char *xml_path;
    FILE *xml_output;
    // print XML
    xml_path = make_path_to_file(opts.outdir, XML_FILENAME);
    if ((xml_output = fopen(xml_path, "w")) == NULL) {
      die("Unable to open \"%s\" for writing.\n", xml_path);
    }
    print_xml_results(
        xml_output,
        argc, argv,
        alph,
        opts.rc,
        opts.dist_name, opts.sig_type_q, opts.sig_thresh,
        opts.bg_file, bg,
        opts.now,
        query_dbs,
        target_dbs,
        match_tmp
        );
    // close matches temporary file
    fclose(match_tmp);
    // close xml output
    fclose(xml_output);
    //generate html from xml
    if (opts.html) {
      char *prog = get_meme_bin_file("tomtom_xml_to_html");
      if (prog != NULL) {
        STR_T *cmd;
        int ret;
        cmd = str_create(0);
        str_append2(cmd, prog);
        str_append(cmd, " ", 1);
        str_append2(cmd, xml_path);
        str_append(cmd, " ", 1);
        str_append_path(cmd, 2, opts.outdir, HTML_FILENAME);

        ret = system(str_internal(cmd));

        if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
          output_error = true;
          fprintf(stderr, "Warning: tomtom_xml_to_html exited abnormally and may "
              "have failed to create HTML output.\n");
        }

        str_destroy(cmd, false);
        free(prog);
      } else {
        output_error = true;
        fprintf(stderr, "Warning: could not find tomtom_xml_to_html. "
            "The HTML output could not be created.\n");
      }
    }
    myfree(xml_path);
  }
  // cleanup
  destroy_mlen(t_mlen);
  arraylst_destroy(destroy_motif_db, query_dbs);
  arraylst_destroy(destroy_motif_db, target_dbs);
  free_motifs(query_motifs);
  free_motifs(target_motifs);
  free_array(bg);
  alph_release(alph);
  cleanup_arguments(&opts);
  return output_error ? EXIT_FAILURE : EXIT_SUCCESS;
}/* Main tomtom*/
#endif

