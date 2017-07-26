#include <stdio.h>
#include "string-list.h"
#include "utils.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

int main(int argc, const char *argv[])
{

  if (argc != 3) {
    char *usage = "Usage: import-html-template <variable-name> <file-name>\n";
    fputs(usage, stderr);
    return 0;
  }

  const char *var_name = argv[1];
  FILE *input = fopen(argv[2], "r");

  const int MAX_TAG_SIZE = 1000;
  char buffer[MAX_TAG_SIZE];
  int c;
  int prev_c = 0;

  STRING_LIST_T *template_tags = new_string_list();

  fprintf(stdout, "char *%s  = \n\"", var_name);
  while ((c = fgetc(input)) != EOF) {
    if (c == '"' || c == '\\') {
      fputc('\\', stdout);
    }
    if (c == '@') {
      int i = 0;
      while ((c = fgetc(input)) != '@') {
        if (i >= MAX_TAG_SIZE) {
          die("Template tag exceeded maximum allowed length.");
        }
        buffer[i] = c;
        ++i;
      }
      buffer[i] = 0;
      add_string(buffer, template_tags);
      fputc('@', stdout);
      fputs(buffer, stdout);
      fputc('@', stdout);
      continue;
    }
    if (c == '\n') {
      fputs("\\n\"\n\"", stdout);
    }
    else {
      fputc(c, stdout);
    }
    prev_c = c;
  }

  fputs("\";\n\n", stdout);

  int i = 0;
  int num_strings = get_num_strings(template_tags);
  fputs("/******************\n", stdout);
  for (i = 0; i < num_strings; ++i) {
    char *template_tag = get_nth_string(i, template_tags);
    fprintf(stdout, "void fimo_print_%s();\n", template_tag);
  }
  fputs("******************/", stdout);

  return 0;
}
