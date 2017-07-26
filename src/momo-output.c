#include "motif-in.h"
#include "config.h"
#include "momo-output.h"
#include "momo-html-string.h"
#include "io.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "ceqlogo.h"
#include "momo-algorithm.h"
#include "momo-modl.h"

const int MAX_HTML_MATCHES = 1000;

/**********************************************************************
 * This function saves MOMO results as a tab-delimited text file
 *********************************************************************/
void print_momo_text_file(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  fprintf(momo_file, "MEME version %s\n\n", VERSION);
  fprintf(momo_file, "Alphabet= %.*s\n\n", 20, summary.alph->symbols + 1);
  fprintf(momo_file, "Background letter frequencies\n");
  
  const char* alph_letters = summary.alph_letters;

  ARRAY_T * bg_freqs = summary.bg_freqs;
  
  int i;
  int j;
  int k;
  int l;
  
  for (i = 0; i < strlen(alph_letters); ++i) {
    fprintf(momo_file, "%c %f ", alph_letters[i], get_array_item_defcheck(i, bg_freqs));
  }
  fprintf(momo_file, "\n\n");

  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  
  for (i = 0; i < arraylst_size(mod_table_keys); i++) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    unsigned long num_passing = (unsigned long) arraylst_size(mod_entry->seq_list);
    if (num_passing >= options.min_occurrences) {
      ARRAYLST_T* motifs = mod_entry->motifinfos;
      for (j = 0; j < arraylst_size(motifs); ++j) {
        MOTIF_INFO_T* motifinfo = arraylst_get(j, motifs);
        MOTIF_T* motif = motifinfo->motif;

        unsigned long motif_list_size = (unsigned long) arraylst_size(motifinfo->seqs);
        
        fprintf(momo_file, "MOTIF %s\n", get_motif_id(motif));
        fprintf(momo_file, "letter-probability matrix: alength= %d w= %d nsites= %lu E= 0\n", (int) (strlen(alph_letters)), options.width, motif_list_size);

        MATRIX_T* freqs = get_motif_freqs(motif);
        for (k = 0; k < options.width; k++) {
          for (l = 0; l < strlen(alph_letters); l++) {
            fprintf(momo_file, "%f\t", get_matrix_cell_defcheck(k, l, freqs));
          }
          fprintf(momo_file, "\n");
        }
        fprintf(momo_file, "\n");
      }
    }
  }
};

void momo_print_version(FILE *momo_file) {
  fprintf(momo_file, "MOMO version %s, (Release date: %s)", VERSION, ARCHIVE_DATE);
};

void momo_print_command_line(FILE *momo_file, MOMO_OPTIONS_T options) {
  fputs(options.command_line, momo_file);
};

void momo_print_parameters(FILE *momo_file, MOMO_OPTIONS_T options) {
  int i;
  
  fprintf(momo_file, "PARAMETERS:\n");
  
  // Algorithm
  ALGORITHM_T algorithm = options.algorithm;
  if (algorithm == simple) {
    fprintf(momo_file, "algorithm: simple\n");
  } else if (algorithm == motifx) {
    fprintf(momo_file, "algorithm: motifx\n");
  } else { // algorithm == modl
    fprintf(momo_file, "algorithm: modl\n");
  }
  
  // Arraylists
  ARRAYLST_T* phospho_filenames = options.phospho_filenames;
  fprintf(momo_file, "phospho filenames: \n");
  for (i = 0; i < arraylst_size(phospho_filenames); ++i) {
    char* phospho_filename = arraylst_get(i, phospho_filenames);
    fprintf(momo_file, "\tfile %d: %s\n", i+1, phospho_filename);
  }
  
  // Booleans
  fprintf(momo_file, "allow clobber: %s\n", (options.allow_clobber ? "true" : "false"));
  fprintf(momo_file, "eliminate repeats: %s\n", (options.eliminate_repeats ? "true" : "false"));
  fprintf(momo_file, "\teliminate repeat width: %d\n", options.eliminate_repeat_width);
  fprintf(momo_file, "filter: %s\n", (options.filter ? "true" : "false"));
  if (options.filter) {
    fprintf(momo_file, "\tfilter field: %s\n", options.filter_field);
    if (options.filter_type == le) {
      fprintf(momo_file, "\tfilter type: <=\n");
    } else if (options.filter_type == lt) {
      fprintf(momo_file, "\tfilter type: <\n");
    } else if (options.filter_type == eq) {
      fprintf(momo_file, "\tfilter type: =\n");
    } else if (options.filter_type == gt) {
      fprintf(momo_file, "\tfilter type: >\n");
    } else if (options.filter_type == ge) {
      fprintf(momo_file, "\tfilter type: >=\n");
    } else {
      fprintf(momo_file, "\tfilter type is unknown!\n");
    }
    fprintf(momo_file, "\tfilter threshold: %f\n", options.filter_threshold);
  }
  fprintf(momo_file, "hash fasta: %s\n", (options.hash_fasta ? "true" : "false"));
  fprintf(momo_file, "\thash fasta width: %d\n", options.hash_fasta_width);
  fprintf(momo_file, "remove unknowns: %s\n", (options.remove_unknowns ? "true" : "false"));
  fprintf(momo_file, "single motif per mass: %s\n", (options.single_motif_per_mass ? "true" : "false"));
  
  // char* values
  fprintf(momo_file, "command line: %s\n", options.command_line);
  fprintf(momo_file, "html path: %s\n", options.html_path);
  fprintf(momo_file, "output dirname: %s\n", options.output_dirname);
  fprintf(momo_file, "protein database filename: %s\n", options.protein_database_filename);
  fprintf(momo_file, "text path: %s\n", options.text_path);

  // const char* values
  fprintf(momo_file, "html filename: %s\n", options.HTML_FILENAME);
  fprintf(momo_file, "text filename: %s\n", options.TEXT_FILENAME);

  // doubles
  fprintf(momo_file, "score threshold: %f\n", options.score_threshold);
  
  // filetypes:
  FILETYPE_T fg_filetype = options.fg_filetype;
  FILETYPE_T bg_filetype = options.bg_filetype;
  if (fg_filetype == psm) {
    fprintf(momo_file, "fg filetype: psm\n");
  } else if (fg_filetype == prealigned) {
    fprintf(momo_file, "fg filetype: prealigned\n");
  } else { // fasta
    fprintf(momo_file, "fg filetype: fasta\n");
  }
  if (bg_filetype == psm) {
    fprintf(momo_file, "There was an error since bg filetype: psm\n");
  } else if (bg_filetype == prealigned) {
    fprintf(momo_file, "bg filetype: prealigned\n");
  } else { // fasta
    fprintf(momo_file, "bg filetype: fasta\n");
  }
  
  // ints
  fprintf(momo_file, "count threshold: %d\n", options.count_threshold);
  fprintf(momo_file, "min occurrences: %d\n", options.min_occurrences);
  fprintf(momo_file, "width: %d\n", options.width);

  fprintf(momo_file, "\n");

  
};

void momo_print_summary(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  fprintf(momo_file, "SUMMARY:\n");
  fprintf(momo_file, "<ul>\n");
  fprintf(momo_file, "<li>Number of Mods: %lu</li>\n", summary.num_mod);
  fprintf(momo_file, "<li>Number of Mod Types: %lu</li>\n", summary.num_modtype);
  fprintf(momo_file, "<li>Number of Mods Passing Filters: %lu</li>\n", summary.num_mod_passing);
  fprintf(momo_file, "<li>Number of Mod Types Passing Filters: %lu</li>\n", summary.num_modtype_passing);
  fprintf(momo_file, "</ul>\n");
  fprintf(momo_file, "<br>\n");
  
  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  
  int i;
  int j;
  int k;
  for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    ARRAYLST_T * motifs = mod_entry->motifinfos;
    for (j = 0; j < arraylst_size(motifs); ++j) {
      MOTIF_INFO_T* currmotifinfo = arraylst_get(j, motifs);
      MOTIF_T * currmotif = currmotifinfo->motif;
      char* motifid = currmotif->id + 1;
      char* motif_name = mm_malloc(strlen(motifid) + 5);
      motif_name[0] = '\0';
      strncat(motif_name, motifid, strlen(motifid));
      strncat(motif_name, ".png", 4);
      
      char* motif_file = mm_malloc(strlen(options.output_dirname) + strlen(motifid) + 2);
      motif_file[0] = '\0';
      strncat(motif_file, options.output_dirname, strlen(options.output_dirname));
      strncat(motif_file, "/", 1);
      strncat(motif_file, motifid, strlen(motifid));
      CL_create1(currmotif, FALSE, FALSE, "MOMO", motif_file, FALSE, TRUE);
      
      if (options.algorithm == motifx) {
        fprintf(momo_file, "final_pattern: %s score: %.12f foreground_matches: %d foreground_size: %d bg_matches: %d bg_size: %d<br>\n", motifid, currmotifinfo->score, currmotifinfo->fg_match, currmotifinfo->fg_size, currmotifinfo->bg_match, currmotifinfo->bg_size);
      } else if (options.algorithm == simple) {
        fprintf(momo_file, "final_pattern: %s foreground_size: %d<br>\n", motifid, currmotifinfo->fg_size);
      } else { // modl
        fprintf(momo_file, "final_pattern: %s score: %f foreground_size: %d<br>\n", motifid, currmotifinfo->score, currmotifinfo->fg_size);
      }
      
      fprintf(momo_file, "<img src=\"%s\">\n<br>\n", motif_name);
      fprintf(momo_file, "<ul>\n");
      for (k = 0; k < arraylst_size(currmotifinfo->seqs); ++k) {
        char* curr_motifinfo_seq = arraylst_get(k, currmotifinfo->seqs);
       fprintf(momo_file, "<li>%s</li>\n", curr_motifinfo_seq);
      }
      fprintf(momo_file, "</ul>\n");
      fprintf(momo_file, "<br>\n");
      
      // cleanup
      myfree(motif_name);
      myfree(motif_file);
    }
    // Print out MoDL log
    if (options.algorithm == modl) {
      fprintf(momo_file, "<b>MoDL Log</b><br>\n");
      fprintf(momo_file, "<ul>\n");
      MATRIX_T* bg_freqs = NULL;
      bg_freqs = get_count_matrix(bg_freqs, mod_entry->bg_seq_list, NULL, &options, &summary);
      for (j = 0; j < options.width; ++j) {
        for (k = 0; k < strlen(summary.alph_letters); ++k) {
          set_matrix_cell_defcheck(j,k,get_matrix_cell_defcheck(j,k,bg_freqs)/arraylst_size(mod_entry->bg_seq_list),bg_freqs);
        }
      }
      
      ARRAYLST_T* modl_ops = mod_entry->modl_ops;
      ARRAYLST_T* temp_list = arraylst_create();
      double minDL = INFINITY;
      for (j = 0; j < arraylst_size(modl_ops); ++j) {
        MODL_STEP_T* step = arraylst_get(j, modl_ops);
        fprintf(momo_file, "<li>STEP: %d, DL: %f<br>\n", j, step->score);
        do_step(step, temp_list, bg_freqs, options.max_motifs, &options, &summary, mod_entry);
        
        if (step->score < minDL) {
          minDL = step->score;
        }
        
        for (k = 0; k < arraylst_size(temp_list); ++k) {
          fprintf(momo_file, "%s<br></li>\n", regexmotif_to_string(arraylst_get(k, temp_list), mod_entry, &summary, &options));
        }
      }
      fprintf(momo_file, "</ul>\n");
      fprintf(momo_file, "Final DL: %f<br>\n", minDL);
      free_matrix(bg_freqs);
    }
  }
  
  fprintf(momo_file, "\n");
  
};


/**********************************************************************
 * This function saves MOMO results as an HTML file
 *********************************************************************/
void print_momo_html_file(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  const int MAX_TAG_SIZE = 1000;
  int html_string_size = strlen(momo_html_string);
  int i = 0;
  for (i = 0; i < html_string_size; ++i) {
    if (momo_html_string[i] != '@') {
      fputc(momo_html_string[i], momo_file);
    }
    else {
      char buffer[MAX_TAG_SIZE];
        ++i;
      int j = 0;
      while (momo_html_string[i] != '@' && j < (MAX_TAG_SIZE - 1)) {
        buffer[j] = momo_html_string[i];
        ++j;
        ++i;
      }
      if (momo_html_string[i] != '@') {
        die("MOMO tag buffer length exceeded\n");
      }
      buffer[j] = '\0';
      if (strcmp("version", buffer) == 0) {
        momo_print_version(momo_file);
      }
      else if (strcmp("command_line", buffer) == 0) {
        momo_print_command_line(momo_file, options);
      }
      else if (strcmp("summary", buffer) == 0) {
        momo_print_summary(momo_file, options, summary);
      }
      else if (strcmp("parameters", buffer) == 0) {
        momo_print_parameters(momo_file, options);
      }
    }
  }
  
}

void create_directory(MOMO_OPTIONS_T options) {
  
  const BOOLEAN_T PRINT_WARNINGS = FALSE;
  if (create_output_directory(
                              options.output_dirname,
                              options.allow_clobber,
                              PRINT_WARNINGS
                              )
      ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }
}

/**********************************************************************
 * This function saves the MOMO results as a set of files in a
 * directory:
 *
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *
 * allow_clobber will determine whether or not existing files will
 * be overwritten.
 *********************************************************************/
void print_momo_results(MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  // Create directory for motifs to have a location
  create_directory(options);
  
  // Print plain text.
  FILE *momo_file = fopen(options.text_path, "w");
  if (!momo_file) {
    die("Couldn't open file %s for output.\n", options.text_path);
  }
  print_momo_text_file(momo_file, options, summary);
  fclose(momo_file);
  
  // Print HTML
  momo_file = fopen(options.html_path, "w");
  if (!momo_file) {
    die("Couldn't open file %s for output.\n", options.html_path);
  }
  print_momo_html_file(momo_file, options, summary);
  fclose(momo_file);
}
