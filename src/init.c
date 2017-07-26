/***********************************************************************
*                                                                      *
* MEME                                                                 *
* Copyright 1994-2015, The Regents of the University of California     *
* Author: Timothy L. Bailey                                            *
*                                                                      *
***********************************************************************/
// init.c 
/*
  Initialize meme.
*/

#define ABS_MIN_W 2

#include "meme.h"
#include "general.h"
#include "banner.h"
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "utils.h"
#include "io.h"
#include "psp.h"
#include "ushuffle.h"

#ifndef EXP
#define EXP 0
#else
#define EXP 1
#endif

// priors 
#define PROTEIN_PLIB "prior30.plib"

#define ROUNDERROR (1E-12)

// User input parameters 
static BOOLEAN check_syntax = FALSE; // exit after checking syntax if true 
static char *datafile = NULL; // positive examples 
static char *sf = NULL; // name to print for datafile 
static char *negfile = NULL; // negative examples 
static char *obj = "classic"; // objective function 
static BOOLEAN use_llr = FALSE; // use likelihood for starts in Classic mode
static char *stat = "mhg"; // statistical test
static char *bfile = NULL; // use default background Markov model file 
static char *pspfile = NULL; // use positional priors 
static BOOLEAN psp2 = FALSE; // true if 2-stranded positional priors 
static char *default_output_dirname = "meme_out";  /* default name of output
                                                   directory */
static BOOLEAN clobber = FALSE; // default is not to overwrite existing files 
static char *mod = "zoops"; // model type input string; default ZOOPS 
static char *alpha = "PROTEIN"; // default alphabet IUPAC protein 1-letter 
static char *alph_file = NULL; // default use built-in alphabet
static BOOLEAN revcomp = FALSE; // don't use reverse complement strand of DNA 
static int pal = 0;   /* = 0, no palindromes
           = 1, force DNA palindromes,
        */
static BOOLEAN ma_trim = TRUE; // trim width using multiple alignment method 
static double wg = 11; // default gap weight 
static double ws = 1; // default space weight 
static BOOLEAN endgaps = TRUE; // count end gaps in multiple alignment 
static double distance = 1e-5; // squared euclidean distance for convergence 
static char *prior = NULL; // prior type input string 
static double beta = -1; // scale factor for prior; defaults differ 
static int nmotifs = 1; // number of motifs to find 
static char *mfile = NULL; // name of known .motifs file
static int maxiter = 50; // max number iterations of EM on best start 
static double nsites = 0; // try one value of nsites0 only if > 0 
static int min_nsites = 0; // minimum nsites0 to try 
static int max_nsites = 0; // maximum nsites0 to try 
static double wnsites = 0.8; // weight on prior on nsites 
static int w = 0; // width of motifs 
static int min_w = MIN_W; // minimum W0 to try 
static BOOLEAN all_widths = FALSE; // all widths between min and max 
static BOOLEAN min_w_set; // if set on command line don't override 
static int max_w = MAX_W; // maximum W0 to try 
static MAP_TYPE map_type; // type of sequence to theta mapping 
static char *mapname = NULL; // map type input string 
static double map_scale=-1;   /* scale of sequence to theta mapping:
          Uni - size of add-n prior (n)
          Pam - PAM distance (120)
           Default set in init_em.
        */
static int n_spcons = 0; // number of specified start points 
static char *spcons[MAXG]; // starting point consensus strings 
static int main_hs = HSIZE; // size of heap at "main" w values 
static double hs_decrease = HS_DECREASE; // Rate of decrease for heap size 
static BOOLEAN x_branch = FALSE; // Use x_branch regardless of seq model 
static BOOLEAN no_x_branch = FALSE; /* Don't use x_branch, regardless of seq
                                   model */
static BOOLEAN w_branch = FALSE; // Controls whether width branching occurs 
static BOOLEAN print_heaps = FALSE; // Print heaps after branching rounds
static BOOLEAN print_pred = FALSE; /* Print out the predicted sites after each
                                   round of MEME search (eg subsequence, EM) */
static BOOLEAN print_pllr = FALSE;
 // print the LLR of the aligned planted sites
static int bfactor = BFACTOR; // branching factor for branching search 
static int maxsize= 100000; // dataset size limit 
static int seed = 0; // random number seed 
static double ctfrac = -1; // fraction of control sequences to use 
static int kmer = 2; // size of kmer to preserve frequency of when shuffling 
static double max_time = 0; // maximum allowed CPU time; ignore if 0 
//FIXME: what should the default be?
static int max_words = -1; // maximum number of words to test (no limit by default)

/***************************************************************************/
/*
        init_meme_background

        Read in the background Markov model or use the residue
        frequencies adjusted by add-one prior.

        Precalculate the log cumulative background probabilities.

        The log probability of any substring of any sequence will then
        be accessible via:
                Log_back(sample[i]->logcumback, j, w)
        where j is the start of the length-w substring of sequence i.
*/
/***************************************************************************/
static void init_meme_background (
  char *bfile, // background model file 
  BOOLEAN rc, // average reverse comps 
  DATASET *dataset // the dataset 
)
{
  int i, j;
  SAMPLE **samples = dataset->samples; // the sequences 
  int n_samples = dataset->n_samples; // # of sequences in dataset 
  ARRAY_T *back; // freq. of tuples 
  int order; // order of background model 

  // read in background Markov model or use dataset frequencies (0-order) 
  if (bfile) {
    // read in background model
    order = -1; // get the highest order available
    back = load_markov_model(dataset->alph, &order, bfile);
  } else {
    // use dataset frequencies + 1 as background model
    double *res_freq, wtr;
    int alen_core;
    res_freq = dataset->res_freq; // weighted residue freq 
    wtr = dataset->wgt_total_res; // weighted residue count 
    alen_core = alph_size_core(dataset->alph);
    order = 0;
    back = allocate_array(alen_core);
    for (i = 0; i < alen_core; i++) {
      // adjust residue frequencies with add-one prior to avoid 0s
      // and round to avoid system dependencies
      //set_array_item(i, (((res_freq[i] * wtr) + 1) / (wtr + alen_core)), back);
      double f; 
      RND(((res_freq[i] * wtr) + 1) / (wtr + alen_core), 8, f);
      set_array_item(i, f, back);
    }
  }
  if (rc) average_rc_markov_model(dataset->alph, order, back);
  // add x-tuples to model
  extend_markov_model(dataset->alph, true, SUM_FREQS, back);
  // normalize for each prefix if a bfile was used
  for (i = 0; i < get_array_length(back); i += alph_size_wild(dataset->alph)) {
    normalize_subarray(i, alph_size_core(dataset->alph), 1e-7, back); 
    set_array_item(i + alph_wild(dataset->alph), 1.0, back);
  }
  // store background model
  dataset->back = back;
  dataset->back_order = order;
 // precalculate the log cumulative background probabilities for each seq 
  dataset->log_total_prob = 0;
  for (i = 0; i < n_samples; i++) { // sequence 
    SAMPLE *s = samples[i]; // sequence 
    // compute probabilites 
    dataset->log_total_prob += calculate_log_cumulative_background(
        dataset->alph, true, dataset->back_order, dataset->back, s->seq,
        s->logcumback);
  } // sequence 
} // init_meme_background 

/**********************************************************************/
/*
  create_priors
*/
/**********************************************************************/
static PRIORS *create_priors(
  PTYPE ptype, // type of prior to use
  // beta for dirichlet priors;
  // < 0 only returns alphabet
  double beta,
  DATASET *dataset, // the dataset
  char *plib_name // name of prior library
)
{
  int i;
  PRIORS *priors;
  priors = (PRIORS *) mm_malloc(sizeof(PRIORS));
  memset(priors, 0, sizeof(PRIORS));
  priors->ptype = ptype;

  // set up the prior counts
  switch (ptype) {
    case Addone: // add one prior
      for (i = 0; i < alph_size_core(dataset->alph); i++) {
        priors->prior_count[i] = 1.0;
      }
      break;
    case Dirichlet: // simple dirichlet prior
      for (i = 0; i < alph_size_core(dataset->alph); i++) {
        priors->prior_count[i] = beta * get_array_item(i, dataset->back);
      }
      break;
    case Dmix: // mixture of dirichlet's
    case Mega: // megaprior heuristic
    case MegaP: // mod. megaprior heuristic
    {
      priors->plib = read_PriorLib(plib_name, beta, dataset->alph);
      // get b=0 prior for modified mega prior heuristic
      if (ptype == MegaP || ptype == Mega) { // used adj freq with Mega
        double b = 0;
        priors->plib0 = read_PriorLib(plib_name, b, dataset->alph);
      }
      break;
    }
  }
  return priors;
} // create_priors


/**********************************************************************/
/*
        init_meme

        Set up all the stuff for meme.
*/
/**********************************************************************/
void init_meme(
  int argc, // number of input arguments 
  char **argv, // input arguments 
  MODEL **model_p, // the model 
  MODEL **best_model_p, // the best model 
  MODEL **neg_model_p, // model of negative examples 
  DATASET **dataset_p, // the dataset 
  DATASET **neg_dataset_p, // dataset of negative examples 
  char *text_filename, // name of the text output file 
  char **output_dirname, // name of the output directory 
  FILE **text_output, // destination for text output
  ARRAYLST_T* seq_array,
  int width,
  BOOLEAN_T eliminate_repeats
)
{
  int i, j, len, pos;
  char cc;
  OBJTYPE objfun = Classic; // type of objective function 
  TESTTYPE test = mHG; // use Fisher exact test
  MOTYPE mtype; // type of model 
  PTYPE ptype; // type of prior 
  PRIORS *priors; // the prior probabilities model 
  P_POINT *p_point; // previously learned starting points 
  ALPH_T *alph; // alphabet for dataset 
  MODEL *model=NULL, *best_model=NULL, *neg_model=NULL;
  DATASET *dataset=NULL, *neg_dataset=NULL;
  double evt = BIG; // no E-value cutoff 
  BOOLEAN no_print = FALSE; // turn off printing if parallel and not main 
  char *plib_name = NULL; // use default library 
  BOOLEAN show_version = FALSE;
  char *np = NULL;

#ifdef PARALLEL
 // turn off printing if parallel and not the main processor 
  no_print = (mpMyID() != 0);
#endif

  *output_dirname = default_output_dirname;

 // get the command line arguments 
  i = 1;
#ifndef lint
 // print the command line 
  argv[0] = "";
  DO_STANDARD_COMMAND_LINE(1,
    USAGE(Usage:\tmeme\t<dataset> [optional arguments]\n);
    NON_SWITCH(1, <dataset> \t\tfile containing sequences in FASTA format\n,
      switch (i++) {
        case 1: datafile = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      });
     FLAG_OPTN(1, h, \t\t\tprint this message, USAGE_MESSAGE);
     DATA_OPTN(EXP, objfun, classic|nc|smhg|cv|nz|ll, \tobjective function (default: classic),
       obj = _OPTION_);
     FLAG_OPTN(EXP, use_llr, \t\tuse LLR (POP) in search for starts,
       use_llr = TRUE);
     DATA_OPTN(EXP, test, mhg|mrs, \t\tstatistical test type (default: mhg),
       stat = _OPTION_);
     DATA_OPTN(EXP, neg, <negdataset>,
       \tfile containing negative example sequences, negfile = _OPTION_);
     DATA_OPTN(1, o, <output dir>,
      \tname of directory for output files\n\t\t\t\twill not replace existing directory,
       *output_dirname = _OPTION_);
     DATA_OPTN(1, oc, <output dir>,
       \tname of directory for output files\n\t\t\t\twill replace existing directory,
       clobber=TRUE; *output_dirname = _OPTION_);
     FLAG_OPTN(1, text, \t\t\toutput in text format (default is HTML),
       TEXT_ONLY = TRUE);
     FLAG_OPTN(1, dna, \t\t\tsequences use DNA alphabet, alpha = "DNA");
     FLAG_OPTN(1, rna, \t\t\tsequences use RNA alphabet, alpha = "RNA");
     FLAG_OPTN(1, protein, \t\tsequences use protein alphabet,
       alpha = "PROTEIN");
     DATA_OPTN(1, alph, <alph file>, \tsequences use custom alphabet,
       alph_file = _OPTION_; alpha = "CUSTOM");
     DATA_OPTN(1, mod, oops|zoops|anr, \tdistribution of motifs, mod = _OPTION_);
     DATA_OPTN(1, nmotifs, <nmotifs>, \tmaximum number of motifs to find,
       nmotifs = atoi(_OPTION_));
     DATA_OPTN(1, evt, <ev>, \t\tstop if motif E-value greater than <evt>,
       evt = atof(_OPTION_));
     DATA_OPTN(1, nsites, <sites>, \tnumber of sites for each motif,
       nsites=atof(_OPTION_));
     DATA_OPTN(1, minsites, <minsites>,
       \tminimum number of sites for each motif, min_nsites=atoi(_OPTION_));
     DATA_OPTN(1,
       maxsites, <maxsites>, \tmaximum number of sites for each motif,
       max_nsites=atoi(_OPTION_));
     DATA_OPTN(1, wnsites, <wnsites>, \tweight on expected number of sites,
       wnsites=atof(_OPTION_));
     DATA_OPTN(1, w, <w>, \t\tmotif width, w = atoi(_OPTION_); min_w_set=TRUE);
     DATA_OPTN(1, minw, <minw>, \t\tminimum motif width,
               min_w = atoi(_OPTION_); min_w_set=TRUE);
     DATA_OPTN(1, maxw, <maxw>, \t\tmaximum motif width,
       max_w = atoi(_OPTION_));
     FLAG_OPTN(1, nomatrim,
       \t\tdo not adjust motif width using multiple\n\t\t\t\talignment,
       ma_trim = FALSE);
     DATA_OPTN(1, wg, <wg>, \t\tgap opening cost for multiple alignments,
       wg=atof(_OPTION_));
     DATA_OPTN(1, ws, <ws>, \t\tgap extension cost for multiple alignments,
       ws=atof(_OPTION_));
     FLAG_OPTN(1, noendgaps, \t\tdo not count end gaps in multiple alignments,
       endgaps = FALSE);
     DATA_OPTN(1, bfile, <bfile>, \tname of background Markov model file,
       bfile = _OPTION_);
     FLAG_OPTN(1, revcomp, \t\tallow sites on + or - DNA strands,
       revcomp = TRUE);
     FLAG_OPTN(1, pal, \t\t\tforce palindromes (requires -dna), pal = 1);
     DATA_OPTN(1, maxiter, <maxiter>, \tmaximum EM iterations to run,
       maxiter = atoi(_OPTION_));
     DATA_OPTN(1, distance, <distance>, \tEM convergence criterion,
       distance = atof(_OPTION_));
     DATA_OPTN(1, psp, <pspfile>, \tname of positional priors file,
       pspfile = _OPTION_);
     DATA_OPTN(EXP, psp2, <pspfile>, \tname of 2-stranded positional priors file\n\t\t\t\t(requires -revcomp),
       pspfile = _OPTION_;psp2 = TRUE);
     DATA_OPTN(1,
       prior, dirichlet|dmix|mega|megap|addone, \n\t\t\t\ttype of prior to use,
       prior = _OPTION_);
     DATA_OPTN(1, b, <b>, \t\tstrength of the prior, beta = atof(_OPTION_));
     DATA_OPTN(1, plib, <plib>, \t\tname of Dirichlet prior file,
       plib_name = strdup(_OPTION_));
     DATA_OPTN(1, maxwords, <maxwords>, \tmaximum number of words to test as EM starts, 
       max_words = atoi(_OPTION_));
     DATA_OPTN(1, spfuzz, <spfuzz>, \tfuzziness of sequence to theta mapping,
       map_scale = atof(_OPTION_));
     DATA_OPTN(1, spmap, uni|pam, \tstarting point seq to theta mapping type,
       mapname = _OPTION_);
     DATA_OPTN(1, cons, <cons>, \t\tconsensus sequence to start EM from,
       spcons[n_spcons++] = _OPTION_);
     DATA_OPTN(1, heapsize, <hs>,
               \tsize of heaps for widths where substring \n\t\t\t\tsearch occurs, 
               main_hs = atoi(_OPTION_));
     FLAG_OPTN(1, x_branch, \t\tperform x-branching, x_branch=TRUE);
     FLAG_OPTN(EXP, no_x_branch, \t\tdo not perform x-branching, 
               no_x_branch=TRUE);
     FLAG_OPTN(1, w_branch, \t\tperform width branching, w_branch=TRUE);
     FLAG_OPTN(1, allw, \t\t\tinclude all motif widths from min to max,
               all_widths=TRUE);
     DATA_OPTN(1, bfactor, <bf>,
       \t\tbranching factor for branching search, bfactor = atoi(_OPTION_));
     FLAG_OPTN(EXP, print_pred, \t\tprint out the sites predicted by meme,
       print_pred = TRUE);
     FLAG_OPTN(EXP, print_heaps, \t\tprint heaps after each branching round,
       print_heaps = TRUE);
     FLAG_OPTN(EXP, planted_LLR, \t\tprint the LLR of the aligned planted sites,
       print_pllr = TRUE);
     DATA_OPTN(1, maxsize, <maxsize>, \tmaximum dataset size in characters,
       maxsize = atoi(_OPTION_));
     FLAG_OPTN(1, nostatus, \t\tdo not print progress reports to terminal,
       NO_STATUS = TRUE);
     DATA_OPTN(1, p, <np>, \t\tuse parallel version with <np> processors, np = _OPTION_);
     DATA_OPTN(1, time, <t>, \t\tquit before <t> CPU seconds consumed,
       max_time = atof(_OPTION_));
     DATA_OPTN(1, sf, <sf>, \t\tprint <sf> as name of sequence file, sf = _OPTION_);
     FLAG_OPTN(2, check_syntax, \t\tcheck input syntax and exit,
       check_syntax = TRUE);
     FLAG_OPTN(1, V, \t\t\tverbose mode, VERBOSE = TRUE);
     FLAG_OPTN(1, version, \t\tdisplay the version number and exit, show_version = TRUE);
     DATA_OPTN(EXP, mfile, <mfile>, \tfile of known motifs, mfile = _OPTION_);
     DATA_OPTN(EXP, seed, <seed>, \t\tseed for random numbers for shuffling and \n\t\t\t\tsampling,
       seed = atoi(_OPTION_));
     DATA_OPTN(EXP, shuf, <kmer>, \t\tpreserve frequencies of k-mers of size <kmer> \n\t\t\t\twhen shuffling,
       kmer = atoi(_OPTION_));
     DATA_OPTN(EXP, ctfrac, <ctfrac>, \tfraction of control sequences to use,
       ctfrac= atof(_OPTION_));
     FLAG_OPTN(EXP, trace, \t\ttrace starting points, TRACE = TRUE);
     FLAG_OPTN(EXP, print_all, \t\tprint all debug information,
       PRINTALL = TRUE);
     FLAG_OPTN(EXP, print_w, \t\tprint erasure matrix, PRINT_W = TRUE);
     FLAG_OPTN(EXP, print_z, \t\tprint missing information matrix,
       PRINT_Z = TRUE);
     FLAG_OPTN(EXP, print_ll, \t\tprint log-likelihood during EM,
       PRINT_LL = TRUE);
     FLAG_OPTN(EXP, print_starts, \t\tprint starting points, 
       PRINT_STARTS = TRUE);
  )
#endif

#ifndef PARALLEL
 // check input arguments 
  if (np != NULL) {
    fprintf(stderr, "-p <np> given but Parallel MEME not configured. Refer to doc/install.html.\n");
    exit(1);
  }
#endif

 // exit if check_syntax is on 
  if (check_syntax) exit(0);

  if (show_version) {
    fprintf(stdout, VERSION "\n");
    exit(EXIT_SUCCESS);
  }

 // set random number generators 
  srand_mt(seed);
  set_randfunc((randfunc_t) random_mt); // for ushuffle

  if (TEXT_ONLY == TRUE) {
    // Legacy: plain text output to standard out.
    *text_output = stdout;
  }
  else {
    if (!no_print) {
      // allow clobbering of the default output directory
      if (*output_dirname == default_output_dirname) { 
        clobber = TRUE;
      } 
      if (create_output_directory(*output_dirname, clobber, !NO_STATUS)) {
        // Failed to create output directory.
        exit(1);
      }
      // Create the name of the text output file 
      // "<dir>/text_filename/" and open it for writing
      char *path = make_path_to_file(*output_dirname, text_filename);
      *text_output = fopen(path, "w"); //FIXME CEG check for errors
      myfree(path);
    }
  }
#ifdef PARALLEL
  // Send all text_output text from non-main node to bit bucket.
  if (mpMyID() != 0) *text_output = fopen("/dev/null", "w");
#endif

  // set all the print flags appropriately 
  if (PRINTALL) {
    PRINT_W = TRUE;
    PRINT_Z = TRUE;
    PRINT_LL = TRUE;
    PRINT_STARTS = TRUE;
  }

   // check that nmotifs >= 1 
  if (nmotifs < 1) {
    fprintf(stderr, "You must specify a minimum of 1 or more motifs if you use -nmotifs.\n");
    exit(1);
  }
  // check that nmotifs < MAXG
  if (nmotifs >= MAXG) {
    fprintf(stderr, "-nmotifs larger than MAXG-1.  Use smaller -nmotifs or recompile with larger MAXG.\n");
   exit(1);
  }

 // get the objective function type 
  if (!strcmp(obj, "classic")) {
    objfun = Classic;
  } else if (!strcmp(obj, "nc")) {
    objfun = NC;
  } else if (!strcmp(obj, "smhg")) {
    objfun = SmHG;
  } else if (!strcmp(obj, "nz")) {
    objfun = NZ;
  } else if (!strcmp(obj, "cv")) {
    objfun = CV;
//FIXME
  } else if (!strcmp(obj, "ll")) {
    objfun = LL;
  } else {
    fprintf(stderr, "Unknown objective function type %s. \n", obj);
    exit(1);
  }

 // Get the statistical test type.
  if (!strcmp(stat, "mhg")) {
    test = mHG;
  } else if (!strcmp(stat, "mrs")) {
    test = mRS;
  } else {
    fprintf(stderr, "Unknown statistical test type %s. \n", obj);
    exit(1);
  }

 // get the model type 
  if (!strcmp(mod, "anr") || !strcmp(mod, "tcm")) {
    mtype = Tcm;
    if (pspfile) { // PM FIXME 
      fprintf (stderr, "-mod anr (-mod tcm) not yet supported for with -psp.\n");
      exit(1);
    }
  } else if (!strcmp(mod, "oops")) {
    mtype = Oops;
  } else if (!strcmp(mod, "zoops")) {
    mtype = Zoops;
  } else {
    mtype = Zoops; // prevent warning 
    fprintf(stderr, "Unknown model type %s. \n", mod);
    exit(1);
  }

 // check ctfrac 
  if (ctfrac == 0 || ctfrac > 1) {
    fprintf(stderr, "ctfrac must be in (0, 1)\n");
    exit(1);
  }
 // ctfrac not given; set defaults
  if (ctfrac == -1) {
    if (objfun == CV || objfun == LL) {
      ctfrac = 0.5;
    } else if (objfun == NZ) {
      ctfrac = 1.0;
    }
  }

  // check the alphabet and set up default mappings and priors 
  if (strcmp(alpha, "DNA") == 0) {
    alph = alph_dna(); // builtin DNA alphabet
    if (!mapname) mapname = "uni"; // uniform prior mapping 
    if (!prior) prior = "dirichlet"; // simple dirichlet prior 
  } else if (strcmp(alpha, "RNA") == 0) {
    alph = alph_rna(); // builtin RNA alphabet
    if (!mapname) mapname = "uni"; // uniform prior mapping 
    if (!prior) prior = "dirichlet"; // simple dirichlet prior 
  } else if (strcmp(alpha, "PROTEIN") == 0) {
    alph = alph_protein(); // builtin Protein alphabet
    if (!mapname) mapname = "pam"; // PAM mapping 
    if (!prior) {
      switch (mtype) {
        case Oops: prior = "dmix"; break;
        case Zoops:
        case Tcm: prior = "megap"; break;
        default: prior = "dirichlet"; break;
      }
    }
  } else if (strcmp(alpha, "CUSTOM") == 0) {
    alph = alph_load(alph_file, !no_print); // load custom alphabet
    if (alph == NULL) exit(EXIT_FAILURE);
    if (!mapname) mapname = "uni";// uniform prior mapping
    if (!prior) prior = "dirichlet"; // simple dirichlet prior
  } else {
    fprintf(stderr, "Unknown alphabet type.\n");
    exit(1);
  }

  // check that the alphabet can perform requested features
  if (!alph_has_complement(alph)) {
    if (revcomp) {
      fprintf(stderr,
        "You must use a complementable alphabet if using -revcomp !\n");
      exit(1);
    }
    if (pal) {
      fprintf(stderr, "You must use a complementable alphabet if using -pal !\n");
      exit(1);
    }
  }

 // find out type of prior 
  if (!strcmp(prior, "dirichlet")) {
    ptype = Dirichlet;
    if (beta < 0) beta = 0.01; // default b = 0.01 
  } else if (!strcmp(prior, "dmix")) {
    ptype = Dmix;
    if (beta < 0) beta = 0; // default b = 0 for dmix 
  } else if (!strcmp(prior, "megadmix") || !strcmp(prior, "mega")) {
    ptype = Mega; // use mega prior heuristic 
  } else if (!strcmp(prior, "megap")) {
    ptype = MegaP; // discretization uses b=0 
  } else if (!strcmp(prior, "addone")) {
    ptype = Addone;
  } else {
    ptype = Dirichlet; // prevent warning 
    fprintf(stderr, "Unknown type of prior: %s!\n", prior);
    exit(1);
  }

  // Read the primary dataset and set up globals.
  dataset = (seq_array != NULL) ? create_meme_dataset_from_momo(seq_array, alph, width, eliminate_repeats) : read_seq_file(datafile, alph, revcomp, FALSE);
  if (!dataset) exit(1);

  // read in the negative dataset if a negative file given 
  if (negfile) {
    if (objfun!=SmHG && objfun!=CV && objfun!=NZ && objfun!=LL) {
      fprintf(stderr, "You must specify '-objfun' with 'nz', 'cv' or 'll' with '-neg'.\n");
      exit(1);
    }
    neg_dataset = read_seq_file(negfile, alph, revcomp, TRUE);
    if (!neg_dataset) exit(1);
  }

 // Create a control dataset if none provided by reading
 // in a new copy of the primary dataset and shuffling preserving 2-mers.
  if (! (objfun == Classic) && ! neg_dataset) {
      if (!NO_STATUS) fprintf(stderr, "Creating shuffled version of primary dataset as control...\n");
    neg_dataset = read_seq_file(datafile, alph, revcomp, TRUE);
 // Note: don't allow kmer<2 or runs of Ns will be broken up in repeatmasked
 // sequences, messing up the statistics.
    if (kmer < 2 || kmer > 6) {
      fprintf(stderr, "-shuf <kmer> : <kmer> must be in the range [2,..,6].\n");
      exit(1);
    }
    shuffle_dataset_letters(neg_dataset, kmer, revcomp);
  }

 // Sort the samples alphabetically by their integer-encoding,
 // accounting for the possibility of two strands, then randomize the order.
  if (neg_dataset) shuffle_dataset_order(neg_dataset);
//FIXME
//for (i=0;i<neg_dataset->n_samples;i++) printf("neg: %s %s\n", neg_dataset->samples[i]->sample_name, neg_dataset->samples[i]->seq);

 // initialize the background model 
  init_meme_background(bfile, revcomp, dataset);
  if (neg_dataset) init_meme_background(bfile, revcomp, neg_dataset);

 // Sort the samples alphabetically by their integer-encoding,
 // accounting for the possibility of two strands, then randomize the order.
  shuffle_dataset_order(dataset);
//FIXME
//for (i=0;i<neg_dataset->n_samples;i++) printf("pos: %s %s\n", dataset->samples[i]->sample_name, dataset->samples[i]->seq);

 // Put the last (after shuffling order) samples in the hold-out groups for CV and LL.
  if (objfun == CV || objfun == LL) {
    int n = dataset->n_samples;
    int h1 = (1-ctfrac) * n; // split 50:50
    int h2 = (1 - ctfrac/2) * n;
    int c1 = h2-h1;
    int c2 = n-h2;
 // Check that there are enough samples left for CV.
//FIXME: how large should the limit be?
#define CVMIN 5
      if (!NO_STATUS) fprintf(stderr, "CV: n %d size1 %d size2 %d\n", n, c1, c2);
    if (n < 3*CVMIN) {
      fprintf(stderr, 
        "You can only use '-cv' if the dataset has at least %d sequences.\n", 3*CVMIN);
      exit(1);
    }
    if (c1 < CVMIN || c2 < CVMIN) {
      double min_ctfrac = (2.0*CVMIN)/n;
      fprintf(stderr, 
        "You must specify '-ctfrac <ctfrac>' at least %g with this dataset.\n", min_ctfrac);
      exit(1);
    }
    dataset->n_group[0] -= c1+c2; // removed from Group 0
    dataset->n_group[1] = c1;
    dataset->n_group[2] = c2; 
    for (i=h1; i<h2; i++) dataset->samples[i]->group = 1;
    for (i=h2; i<n; i++) dataset->samples[i]->group = 2;
  } else if (objfun == NZ) {
 //
 // Add the first ctfrac of the control sequences to the 
 // primary dataset.
 // Note: this should be done after setting max_nsites if we
 // don't want to bias the search.  On the other hand? FIXME
 //
    int n = negfile ? neg_dataset->n_samples : dataset->n_samples;
    int nct = (ctfrac/2) * n; // split 50:50
 // Check that there are enough samples for NZ.
//FIXME: how large should the limit be?
#define NZMIN 2
    if (!NO_STATUS) fprintf(stderr, "NZ: n %d size %d\n", n, nct);
    if (n < 2*NZMIN) {
      if (negfile) {
	fprintf(stderr, 
	  "You can only use '-nz' if the control dataset has at least %d sequences.\n", 2*NZMIN);
      } else {
        fprintf(stderr, 
	  "You can only use '-nz' if the dataset has at least %d sequences.\n", 2*NZMIN);
      }
      exit(1);
    }
    if (nct < NZMIN) {
      double min_ctfrac = (2.0*NZMIN)/n;
      fprintf(stderr, 
        "You must specify '-ctfrac <ctfrac>' at least %g with this dataset.\n", min_ctfrac);
      exit(1);
    }
    add_control_samples(dataset, neg_dataset, nct, nct, revcomp);
    neg_dataset = NULL; // Prevent it from ever being freed!
  }

 // Set the number of primary samples.
  int n_primary_samples = (objfun==NZ) ? dataset->n_group[0] : dataset->n_samples;

 // Set state variables in dataset(s).
  if (objfun == CV || objfun == LL) {
    set_seq_groups_to_skip(dataset, FALSE, TRUE, TRUE);
  } else if (objfun == NZ) {
    set_seq_groups_to_skip(dataset, FALSE, FALSE, FALSE);
  } else {
    set_seq_groups_to_skip(dataset, FALSE, FALSE, FALSE);
  }
//FIXME: I think this is unneeded.
 //if (neg_dataset) {
 //  set_seq_groups_to_skip(neg_dataset, FALSE, FALSE, FALSE);
 //}

 // check we can actually use the dataset in this mode
  if (n_primary_samples == 1 && mtype != Tcm) {
    fprintf(stderr, "You must specify '-mod anr' to set the motif site model "
        "to 'Any Number of Repetitions [per sequence]' when only providing "
        "one sequence\n");
    exit(1);
  }

  // read in psp file if one given
  if (pspfile) {
    read_psp_file(pspfile, dataset, psp2, revcomp, mtype);

    // warn that we are using the W in the PSP file as minw 
    if (!min_w_set) {
      fprintf(stderr,"Setting minimum motif width to width of the prior in the PSP file: %d\n",
        dataset->psp_w);
      min_w = dataset->psp_w;
    }
  } else { // no PSP file 
    dataset->psp_w = min_w;
  }

  // Set the objective function.
  dataset->objfun = objfun;
  dataset->test = test;
  if (neg_dataset) neg_dataset->objfun = objfun;
  if (neg_dataset) neg_dataset->test = test;


  // prevent too long jobs 
  if (dataset->total_res > maxsize) {
    fprintf(stderr, "Dataset too large (> %d).  Rerun with larger -maxsize.\n",
      maxsize);
    exit(1);
  }

 // read in known motifs file 
  if (mfile) {
    FILE  *fdata = fopen(datafile, "r"); //FIXME CEG check for errors
    dataset->nkmotifs = read_motifs(alph, fdata, mfile, dataset->motifs, FALSE, NULL);
    nmotifs = dataset->nkmotifs;
    dataset->pal = 0; // no palindrome testing 
  } else {
    dataset->nkmotifs = 0;
  }

  // create the priors 
  if (ptype == Dmix || ptype == Mega || ptype == MegaP) {
    // make the name of the prior library 
    if (!plib_name) {
      if (alph_is_builtin_protein(alph)) { // default mixture prior for proteins
        plib_name = make_path_to_file(get_meme_etc_dir(), PROTEIN_PLIB);
      } else {
        fprintf(
          stderr, 
          "WARNING: When using DNA or a custom alphabet, "
          "and specifiying a prior type of\n"
          "'dmix', 'mega' or 'megap', a prior library must be provided.\n"
          "No prior library was provided, so a simple Dirichlet prior will be used.\n"
        );
        prior = "dirichlet";
        ptype = Dirichlet;
        if (beta <= 0) beta = 0.01; // default b = 0.01 for simple Dirichlet
      }
    }
  }
  if ((ptype == Mega || ptype == MegaP) && beta == -1) {
    // tlb 5-9-97; wgt_total_res 
    //beta = 10.0 * dataset->wgt_total_res; // size of mega prior 
    beta = 5.0 * dataset->wgt_total_res; // size of mega prior 
  }
  priors = create_priors(ptype, beta, dataset, plib_name);

  // set number of occurrences of sites 
  if (nsites != 0) {
    if (mtype == Oops) {
      fprintf(stderr, "You may not specify -sites with -mod oops\n");
      exit(1);
    }
    min_nsites = max_nsites = nsites;
  }

  // set search range for nsites 
  if (mtype == Oops) {
    min_nsites = max_nsites = n_primary_samples;
  } else if (mtype == Zoops) {
    if (min_nsites > n_primary_samples) {
      fprintf(stderr, "Minimum number of sites too large.  Setting to 2.\n");
      min_nsites = 2;
    }
    if (!min_nsites) min_nsites = 2; // default 
    if (max_nsites > n_primary_samples) {
      fprintf(stderr, "Maximum number of sites exceeded.  Setting to %d.\n",
        n_primary_samples);
      max_nsites = n_primary_samples;
    }
    if (!max_nsites) max_nsites = n_primary_samples; // default 
  } else { // TCM model 
    if (min_nsites<2) min_nsites = 2; // default 
    if (!max_nsites) max_nsites = MIN(5*n_primary_samples, 50);
  }

  // check that max number of sites >= min number of sites 
  if (min_nsites > max_nsites) {
    fprintf(
      stderr, 
      "The minimum number of sites is set to %d. "
      "It should be less than the max number of sites (%d).\n",
      min_nsites,
      max_nsites
    );
    exit(1);
  }
  // check that there are enough possible sites 
  if (min_nsites < 2) {
    fprintf(stderr, "You must specify a minimum of 2 or more sites.\n");
    exit(1);
  }
  if (max_nsites < 2) {
    fprintf(stderr, "It must be possible for at least 2 sites to fit.\n");
    exit(1);
  }

  // check weight on prior on nsites 
  if (wnsites >= 1 || wnsites < 0) {
    fprintf(stderr, "<wnsites> must be in range [0..1).\n"); exit(1);
  }

  // set up globals 
  if (w != 0) { // w specified; set min_w and max_w 
    max_w = min_w = w;
    fprintf(stderr,"w set, setting max and min widths to %d\n",w);
  }

  // check that no sequence too short 
  //if (dataset->min_slength < MIN_W) {
  if (dataset->min_slength < min_w) {
    fprintf(stderr,
     "All sequences must be at least %d characters long.  Set -w or -minw or remove ",
        //MIN_W);
        min_w);
    fprintf(stderr, "shorter sequences\nand rerun.\n");
    exit(1);
  }

  // oops model: limit max_w to shortest seq 
  if (mtype == Oops && max_w > dataset->min_slength) {
    max_w = dataset->min_slength;
    fprintf(stderr,
      "maxw > length of shortest sequence (%ld).", dataset->min_slength);
    fprintf(stderr, "  Setting maxw to %d.\n", max_w);
  }
  // all models: limit max_w to longest seq 
  if (max_w > dataset->max_slength) {
    max_w = dataset->max_slength;
    fprintf(stderr,
      "maxw > length of longest sequence (%ld).", dataset->max_slength);
    fprintf(stderr, "  Setting maxw to %d.\n", max_w);
  }
  if (max_w > MAXSITE) {
    fprintf(stderr,
      "maxw too large (> %d).  Recompile with larger MAXSITE.\n", MAXSITE);
    exit(1);
  }
  if (max_w < 0) { // use default 
    max_w = MIN(MAXSITE, dataset->min_slength); // maximum W0 
  }

  // check that min_w <= max_w 
  if (min_w > max_w) {
    if (pspfile) {
      fprintf(stderr, "PSP file w = %d > maxw = %d. Respecify larger -maxw.\n",
              dataset->psp_w, max_w);
      exit(1);
    }
     fprintf(stderr, "minw > maxw.  Setting minw to %d.\n", max_w);
     min_w = max_w;
  }

  // check that min_w and max_w are at least ABS_MIN_W 
  if (min_w < ABS_MIN_W) {
    fprintf(stderr,
      "Minimum width must be >= %d.  Respecify larger -w or -minw.\n",
      ABS_MIN_W);
    exit(1);
  } else if (max_w < ABS_MIN_W) {
    fprintf(stderr,
      "Maximum width must be >= %d.  Respecify larger -w or -maxw.\n",
      ABS_MIN_W);
    exit(1);
  }

  // must use TCM if only one sequence 
  if (mtype != Tcm && n_primary_samples==1) {
    fprintf(stderr,
      "You must specify '-mod anr' since your dataset contains only one sequence.\n"
    );
    fprintf(stderr,
      "Alternatively, you might wish to break your sequence into several sequences.\n"
    );
    exit(1);
  }

  // flag search for palindromes 
  dataset->pal = pal;
  if (neg_dataset) neg_dataset->pal = pal;

  // check that IUPAC alphabet if using PAM mapping 
  // mapname == "pam" && alphabet != PROTEIN0
  if (strcmp(mapname, "pam") == 0 && !(alph_is_builtin_protein(alph) || alph_is_builtin_dna(alph))) {
    fprintf(stderr,
     "Setting sequence to theta mapping type to `uni' since the alphabet is "
     "not built-in and 'pam' is only supported for the built-in alphabets.\n");
    mapname = "uni";
  }
  // get the type of mapping between sequences and thetas 
  if (strcmp(mapname, "uni") == 0) {
    map_type = Uni;
    if (map_scale == -1) map_scale = .5; // default add .5 
  } else if (strcmp(mapname, "pam") == 0) {
    map_type = Pam;
    if (map_scale == -1) map_scale = 120; // default PAM 120 
  } else {
    fprintf(stderr, "Unknown mapping type %s. \n", mapname);
    exit(1);
  }

  // set up the sequence to theta mapping matrix for starts 
  dataset->map = init_map(map_type, map_scale, alph, dataset->back, FALSE);
  dataset->lomap = init_map(map_type, map_scale, alph, dataset->back, TRUE);

  // set up p_point: previously learned components start points 
  p_point = (P_POINT *) mymalloc(sizeof(P_POINT));
  p_point->c = n_spcons; // number of specified starts 
  // the default starting points 
  for (i=0; i<MAXG; i++) {
    p_point->e_cons0[i] = NULL;
    p_point->w[i] = max_w;
    p_point->nsites[i] = nsites;
  }

  // setup user-specified start points 
  for (i = 0; i < n_spcons; i++) {
    uint8_t *e_cons = mm_malloc(sizeof(uint8_t) * MAXSITE); // encoded version 
    if (i >= MAXG) {
      fprintf(stderr, "Too many starting points.  Increase MAXG and recompile.");
      exit(1);
    }
    p_point->e_cons0[i] = e_cons;
    // encode as integer 
    for (j = 0; (cc = spcons[i][j]) != '\0'; j++) {
      if (alph_is_known(alph, cc)) {
        // encode to core + wildcard
        e_cons[j] = alph_encodec(alph, cc);
      } else {
        fprintf(stderr, "Illegal letter %c in consensus string!\n", cc);
        exit(1);
      }
    }
    // set width to length of consensus unless -w specified 
    if (w == 0) p_point->w[i] = j;
    // pad out the consensus sequence with A's 
    for ( ; j < MAXSITE; j++) e_cons[j] = 0;
  }

  // Setup heap size for storage of starting points:
  if (main_hs < 1) {
    fprintf(stderr, "Heap size must be >= 1.\n");
    exit(1);
  } else {
    dataset->main_hs = main_hs;
    dataset->hs_decrease = hs_decrease;
  }

  // Setup a struct recording the desired parameters for branching search:
  BRANCH_PARAMS *branch_params = NULL;
  branch_params = mymalloc(sizeof(BRANCH_PARAMS));

  // Setup branching factor for branching search:
  if (bfactor < 1) {
    fprintf(stderr, "bfactor must be >= 1.\n");
    exit(1);
  } else {
    branch_params->bfactor = bfactor;
  }

  // Record whether the user wants w-branching carried out:
  branch_params->w_branch = w_branch;

  // Record whether the user wants x-branching carried out:
  if (x_branch){
    branch_params->point_branch = X_ONLY;
  } else { 
    branch_params->point_branch = NO_POINT_B;
  }

  if (FALSE) {
  // FIXME: This code can be used to set branching as the default 
  // for different alphabets/sequence models using the -x_branch
  // and -no_x_branch switches. Currently the default is set to
  // no branching (in previous if statement).
  /* Record whether the user wants x_branching carried out. The result
     is recorded in the "point_branch" attribute of branch_params. Note
     that the user is only able to specify x_branching or no x_branching.
     This is because we found ACGT (ie regular) branching to be of little
     benefit. */
  if (x_branch) {
    // User has specified they want x-branching... 
    if (no_x_branch) {
      // Invalid to specify x_branch and no_branch simultaneously:
      fprintf(stderr, "x_branch and no_branch cannot be specified"\
                      "at the same time");
      exit(1);
    }
    else {
      branch_params->point_branch = X_ONLY;
    }
  } else if (no_x_branch) {
    branch_params->point_branch = NO_POINT_B;
  } else {
    // User did not specify x_branching => Decide based on sequence model:
    // NOTE: branching is the default for oops for DNA only
    if (mtype == Oops && !strcmp(alpha, "DNA")) {
      branch_params->point_branch = X_ONLY; // Only x_branch under oops by default.
    } else {
      branch_params->point_branch = NO_POINT_B; // Only x_branch under oops by default.
    }
  } // Deciding branch_params->x_branch
  } // end if 

  // Store the branching parameters for future reference:
  dataset->branch_params = branch_params;

  // Print sites predicted by MEME:
  dataset->print_pred = print_pred;

  // Record whether heaps are to be printed:
  dataset->print_heaps = print_heaps;

  // Record whether llr is to be printed for the alignment of planted sites:
  if (print_pllr) {
    if (w <= 0) {
      fprintf(stderr,
              "Not valid to request llr of planted sites unless length of sites"
              " has been specified.\n");
      exit(1);
    }
  }
  dataset->print_pllr = print_pllr;

  // make sure nmotifs is as large as the number of starting points 
  if (nmotifs < n_spcons) {
    nmotifs = n_spcons;
    fprintf(stderr, "Setting nmotifs to %d\n", nmotifs);
  }

  // create the model 
  model = create_model(mtype, revcomp, max_w, alph, objfun);
  best_model = create_model(mtype, revcomp, max_w, alph, objfun);

  // initialize log and exp lookup tables 
  init_log();
  init_exp();

  // Initialize the probability tables for the objective function.  
  // Get for 2 and up because weighting may cause < min_nsites sites 
  if (objfun == Classic) {
    // Initialize the probability tables for the objective function.
    init_llr_pv_tables(2, max_nsites, alph_size_core(alph), dataset->back, dataset->pal);
  }

  // Set the high-water mark for restricting the total number of seeds tested.
  // Note: Does not include holdout sets in CV mode.
  // Note: uses min_w as requested by user which may get reset below.
//FIXME:
  int n_seed_seqs = (objfun==CV || objfun==LL) ? dataset->n_group[0] : dataset->n_samples;
  dataset->last_seed_seq = n_seed_seqs - 1;
  dataset->last_seed_pos = dataset->samples[n_seed_seqs-1]->length;
  dataset->max_words = max_words;
  if (max_words > 0 && dataset->total_res > max_words) {
  // Find where the last seed word of the minimum width would start.
    int total_res = 0;
    for (i=0; i<n_seed_seqs; i++) {
      int slen = dataset->samples[i]->length;
      if (total_res + slen > max_words) {
        dataset->last_seed_seq = i;
        dataset->last_seed_pos = max_words - total_res - 1;
        break;
      }
      total_res += slen-min_w+1;
    }
    n_seed_seqs = dataset->last_seed_seq + 1;
  }
//FIXME
  if (!NO_STATUS) fprintf(stderr, "SEEDS: highwater mark: seq %d pos %d\n", 
    dataset->last_seed_seq, dataset->last_seed_pos);

  // Load balance the search for starting points in parallel mode.
  // Note: This must come after adding in control samples.
  balance_loop(dataset->samples, n_seed_seqs);

 // set up scratch models 
  model->pal = pal; //FIXME jj - this was uninitialised
  model->min_w = min_w;
  model->max_w = max_w;
  model->all_widths = all_widths;
  model->min_nsites = min_nsites;
  model->max_nsites = max_nsites;
  copy_model(model, best_model, alph);

 // put meme parameters in dataset 
  dataset->use_llr = use_llr;
  dataset->priors = priors;
  dataset->p_point = p_point;
  dataset->wg = wg;
  dataset->ws = ws;
  dataset->endgaps = endgaps;
  dataset->wnsites = wnsites;
  dataset->ma_adj = ma_trim;
  dataset->distance = distance;
  dataset->nmotifs = nmotifs;
  dataset->maxiter = maxiter;
  dataset->evt = evt;
  dataset->mod = mod;
  dataset->mapname = mapname;
  dataset->map_scale = map_scale;
  dataset->priorname = prior;
  dataset->beta = beta;
  dataset->seed = seed;
  dataset->ctfrac = ctfrac;
  dataset->shuffle = (! (objfun == Classic) && ! neg_dataset) ? kmer : -1;
 // save name of prior library 
  if (plib_name) {
    char *tmp, *base; // remove the directory, leaving just the file name
    for (tmp = base = plib_name; *tmp; tmp++) if (*tmp == '/') base = tmp + 1;
    dataset->plib_name = strdup(base); // copy the name
  } else {
    dataset->plib_name = NULL;
  }
  dataset->datafile = sf ? sf : datafile; // name to print 
  dataset->negfile = negfile;
  dataset->bfile = bfile;
  dataset->max_time = max_time;

  // save command line 
  dataset->command = NULL;
  argv[0] = "meme";
  /*for (i=pos=len=0; i<argc-2; i++) {*/ // don't save last arguments 
  // save all arguments; used to not save last 2; why??? 
  for (i=pos=len=0; i<argc; i++) {
    len += strlen(argv[i])+1; // +1 for space following 
    Resize(dataset->command, len+2, char); // +1 for null 
    strcpy((dataset->command)+pos, argv[i]);
    dataset->command[len-1] = ' ';
    dataset->command[len] = '\0';
    pos = len;
  }
  if (TEXT_ONLY) {
    dataset->output_directory = NULL;
  }
  else {
    dataset->output_directory = *output_dirname;
  }

  // set up return values 
  *model_p = model;
  *best_model_p = best_model;
  *neg_model_p = neg_model;
  *dataset_p = dataset;
  *neg_dataset_p = neg_dataset;

  // announce meme 
  banner("MEME", *text_output);
  fprintf(*text_output, "\n\n");

  // cleanup 
  free(plib_name);

} // init_meme 
