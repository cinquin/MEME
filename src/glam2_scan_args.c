#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>  /* non-ANSI */
#include <getopt.h>
#include "glam2_util.h"
#include "glam2_scan_args.h"
#include "config.h"

/* default values: */
static const char* default_output_dirname = "glam2scan_out";
static const int n_def = 25;
static const double D_def = 0.1;
static const double E_def = 2;
static const double I_def = 0.02;
static const double J_def = 1;

static void usage(void) {
  die("\
Usage: glam2scan [options] alphabet my_motif.glam2 my_seqs.fa\n\
Main alphabets: p = proteins, n = nucleotides\n\
Main options (default settings):\n\
-h: show all options and their default settings\n\
-o: output directory; will not clobber existing files\n\
-O: output directory (%s); allow clobbering\n\
-t: output text only to stdout\n\
-n: number of alignments to report (%d)\n\
-2: examine both strands - forward and reverse complement\n\
-v: print version and exit (also accepts --version)\n\
", default_output_dirname, n_def);
}

static void help(void) {
  fprintf(stderr, "\
Usage: glam2scan [options] alphabet my_motif.glam2 my_seqs.fa\n\
Alphabets: p = proteins, n = nucleotides, other = alphabet file\n\
Options (default settings):\n\
-h: show all options and their default settings\n\
-o: output directory; will not clobber existing files\n\
-O: output directory (%s); allow clobbering\n\
-t: output text only to stdout\n\
-n: number of alignments to report (%d)\n\
-2: examine both strands - forward and reverse complement\n\
-D: deletion pseudocount (%g)\n\
-E: no-deletion pseudocount (%.1f)\n\
-I: insertion pseudocount (%g)\n\
-J: no-insertion pseudocount (%.1f)\n\
-d: Dirichlet mixture file\n\
-v: print version and exit (also accepts --version)\n\
", default_output_dirname, n_def, D_def, E_def, I_def, J_def);
  exit(0);
}

void getargs(args *a, int argc, char **argv) {
  int c;
  bool o_specified = false;

  a->text_only = false;
  a->out_dir = (char *) default_output_dirname;
  a->hit_num = n_def;
  a->two_strands = 0;
  a->delete_pseudo = D_def;
  a->no_delete_pseudo = E_def;
  a->insert_pseudo = I_def;
  a->no_insert_pseudo = J_def;
  a->dirichlet_file = NULL;

  static struct option long_options[] = {
    {"version", no_argument,  0, 'v'},
    {0,         0,            0,  0 }
  };
  /* non-ANSI: */
  while ((c = getopt_long(argc, argv, "ho:O:tn:2D:E:I:J:d:v", long_options, NULL)) != -1) {
    switch (c) {
    case 'h':
      help();
    case 'o':
      a->out_dir = optarg;
      a->clobber = 0;
      o_specified = true;
      break;
    case 'O':
      a->out_dir = optarg;
      a->clobber = 1;
      o_specified = true;
      break;
    case 't':
      a->text_only = true;
      break;
    case 'n':
      a->hit_num = xatoi(optarg);
      if (a->hit_num < 0)  /* let's allow 0, even though it's stupid */
	die("%s: option -n should be at least 0\n", prog_name);
      break;
    case '2':
      a->two_strands = 1;
      break;
    case 'D':
      a->delete_pseudo = xatof(optarg);
      if (a->delete_pseudo <= 0)
	die("%s: option -D should be > 0\n", prog_name);
      break;
    case 'E':
      a->no_delete_pseudo = xatof(optarg);
      if (a->no_delete_pseudo <= 0)
	die("%s: option -E should be > 0\n", prog_name);
      break;
    case 'I':
      a->insert_pseudo = xatof(optarg);
      if (a->insert_pseudo <= 0)
	die("%s: option -I should be > 0\n", prog_name);
      break;
    case 'J':
      a->no_insert_pseudo = xatof(optarg);
      if (a->no_insert_pseudo <= 0)
	die("%s: option -J should be > 0\n", prog_name);
      break;
    case 'd':
      a->dirichlet_file = optarg;
      break;
    case 'v':
      fprintf(stderr, VERSION "\n");
      exit(EXIT_SUCCESS);
      break;
    case '?':
      usage();
    }
  }

  if (optind != argc-3)
    usage();

  if (o_specified && a->text_only)
    die("option -t cannot be used with options -o or -O\n", prog_name);

  a->alph_name = argv[optind++];
  a->motif_file = argv[optind++];
  a->seq_file = argv[optind++];
}

void printargs(FILE *fp, int argc, char **argv) {
  int i;
  for (i = 0; i < argc; ++i) {
    fputs(argv[i], fp);
    putc(i < argc-1 ? ' ' : '\n', fp);
  }
}
