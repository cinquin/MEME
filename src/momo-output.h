#include "array-list.h"
#include "momo.h"
#include "alphabet.h"

/**********************************************************************
 * This function saves MOMO results as a tab-delimited text file
 *********************************************************************/
void print_momo_text_file(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary);

void momo_print_version(FILE *momo_file);

/**********************************************************************
 * This function saves MOMO results as an HTML file
 *********************************************************************/
void print_momo_html_file(
                          FILE *momo_file,
                          MOMO_OPTIONS_T options,
                          SUMMARY_T summary
                          );

void momo_print_command_line(FILE *momo_file, MOMO_OPTIONS_T options);
void momo_print_parameters(FILE *momo_file, MOMO_OPTIONS_T options);
void momo_print_summary(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary);

void create_directory(MOMO_OPTIONS_T options);

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
void print_momo_results(MOMO_OPTIONS_T options, SUMMARY_T);

