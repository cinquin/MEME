/*#define DEBUG*/
/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 2000-2015, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
#include "meme.h"
#include "dpalign.h"
#include "fisher_exact.h"
#include "ranksum_test.h"

/* rounding stuff */
#define RNDEPS 1e-12

//FIXME:
//BOOLEAN debug = TRUE;
BOOLEAN debug = FALSE;

static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */
  int w,                                        /* width of sites */
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off 					/* best offset */
);

// For sorting column indices by column RE.
typedef struct pair {int ind; double re;} PAIR;
static int pairCompare(PAIR *e1, PAIR *e2) {
  // compares the RE of two pairs (used for ascending qsort order)
  return (e1->re < e2->re) ? -1 : (e1->re > e2->re) ? 1 : e1->ind - e2->ind;
}

// Keep around between invocations.
static THETA saved_theta = NULL;

/***********************************************************************/
/*
	print_maxima
*/
/***********************************************************************/
static void print_maxima(
  int n_maxima,
  MODEL *model
)
{
#ifdef PARALLEL
#undef printf
#endif
  int i;
  printf("\n");
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    int y = model->maxima[i].y;
    BOOLEAN neg = model->maxima[i].negative;
    int rank = model->maxima[i].rank;
    double prob = model->maxima[i].prob;
    printf("max %d %15.12f x %d y %d set %s rank %d\n", 
      i+1, prob, x, y, neg ? "control" : "primary", rank);
  }
  fflush(stdout);
#ifdef PARALLEL
#define printf if (mpMyID() == 0) printf
#endif
} // print_maxima

/***********************************************************************/
/*
	get_best_nsites

	Find the value of nsites that minimizes the significance of the
	log likelihood ratio (classic) or Fisher exact test.
	If using Fisher exact test, break ties using largest nsites, then
        largest LLR.

	Returns nsites; sets the column scores, best log p-value,
        best log e-value and its LLR.
*/
/***********************************************************************/
static int get_best_nsites (
  MODEL *model,					/* the model */
  DATASET *dataset,				/* the dataset */
  int min_nsites,				/* minimum sites */
  int max_nsites,				/* minimum sites */
  int w, 	 				/* width of motif */
  int n_maxima,					/* number of maxima */
  int n_neg_maxima,				/* for use with ranks */
  P_PROB maxima,				/* the maxima positions */
  double *col_scores,				/* column scores */
  double *best_log_pv,				/* best p-value */
  double *best_log_ev,				/* best E-value */
  double *best_llr  				/* LLR of best p-value */
)
{
  int i;
  MOTYPE mtype = model->mtype;			/* type of model */
  int n_nsites; 				/* number of different nsites */
  int nsites;					/* number of sites */
  int best_nsites;				// best value of nsites
  S_POINT *s_points = NULL;			/* for use with align_top...*/
  BOOLEAN classic = (dataset->objfun==Classic);	// Classic meme
  BOOLEAN use_ranksum = (dataset->test==mRS);

  // Return p-value = 1 if min_nsites too large.
  if (min_nsites > n_maxima) {
    *best_log_pv = *best_log_ev = *best_llr = 0;
    return min_nsites;
  }

  /* limit maximum nsites */
  max_nsites = MIN(n_maxima, max_nsites);	
  n_nsites = max_nsites-min_nsites+1;	

  /* create array for use by align_top_subsequences */
  Resize(s_points, n_nsites, S_POINT);

  /* initialize the s_points for all widths in range [min_nsites,max_nsites] */
  for (i=0, nsites=min_nsites; i<n_nsites; i++, nsites++) {
    s_points[i].nsites0 = nsites;		/* number of sites */
    s_points[i].score = LITTLE;			/* no score yet */
    s_points[i].evaluate = TRUE;                /* Evaluate at every s_point */
  }

  /* get the probability that a site starting at position x_ij would
     NOT overlap a previously found motif; used in E_STEP.
  */
  get_not_o(dataset, model->w);

  // TLB: Works better without the following line on Yeast examples.
  //add_psp_to_log_not_o(dataset, model->w, model->invcomp, model->mtype);

  /* 
     align the top nsites sorted subsequences and compute the 
     log_pop or LLR function on each alignment with nsites in [min_nsites,max_nsites]
  */
  // NOTE: If not Classic, I think this call is only needed to get the LLR.
  BOOLEAN orig_use_llr = dataset->use_llr;	// What to compute on alignment?
  dataset->use_llr = FALSE;			// Causes log_pop to be used in Classic mode.
  (void) align_top_subsequences(mtype, w, dataset, 0, 0, 0, 0, n_nsites, 
    n_maxima, maxima, col_scores, s_points);
  dataset->use_llr = orig_use_llr;		// Restore setting.

  /* 
     determine the significance of the score for each number of 
     sites and chose the number of sites to minimize it 
  */
  *best_llr = -BIG;
  *best_log_ev = BIG;
  best_nsites = 0;
  for (i=0; i<n_nsites; i++) {			/* starting points */
    double score = s_points[i].score;		/* Classic: -log_pop; Otherwise: LLR */
    int N = s_points[i].nsites0;
    double wN = s_points[i].wgt_nsites; 
    double log_pv, log_ev;
    if (classic || dataset->objfun == NC) {
      log_ev = get_log_sig(score, model->mtype, w, wN, N, model->invcomp,
        model->pal, dataset);
      RND(log_ev, RNDDIG, log_ev);
      log_pv = log_ev;				// Not real p-value.
//FIXME
    } else if (dataset->objfun == LL && n_neg_maxima == 0) { 
      log_pv = log_ev = -s_points[i].score;	// -LLR
    } else {
      if (model->mtype == Oops || use_ranksum) {
        // OOPS model uses ranksum test
	int rank = model->maxima[N-1].rank;
        double ta_obs = model->maxima[N-1].ranksum;// sum of positive ranks
        int n = n_maxima + n_neg_maxima;
        int na = N;				// number of positives
        RSR_T *r = ranksum_from_stats(n, na, ta_obs);
        log_pv = log(RSR_get_p_left(r));
	RND(log_pv, RNDDIG, log_pv);
	log_ev = log_pv;			// Not real E-value.
//FIXME
//if (! NO_STATUS) printf("get_best_nsites: w %d n %d na %d ta_obs %g log_pv %g pv %g LLR %g log_ev %g\n", w, n, na, ta_obs, log_pv, exp(log_pv), score, log_ev);
      } else {
	// non-OOPS models use mHG unless ranksum requested
	int rank = model->maxima[N-1].rank;
	int pos = n_maxima;
	int neg = n_neg_maxima;
	int pos_succ = N;
	int neg_succ = rank - N;
	// Fisher Exact test.  Use slow version to get p-values close to 1.
	log_pv = getLogFETPvalue(pos_succ, pos, neg_succ, neg, FALSE);
	RND(log_pv, RNDDIG, log_pv);
	log_ev = log_pv;				// Not real E-value.
//FIXME
//if (! NO_STATUS) printf("get_best_nsites: w %d nsites %d rank %d pos_succ %d pos %d neg_succ %d neg %d log_pv %g pv %g LLR %g log_ev %g\n", w, N, rank, pos_succ, pos, neg_succ, neg, log_pv, exp(log_pv), score, log_ev);
      }
    }

    if (TRACE) printf("w %d N %d wN %f log_ev %f \n", w, N, wN, log_ev);

    /* save EV if best so far, breaking ties using width then LLR */
    if ( 
      (classic 
        && RNDEPS < *best_log_ev - log_ev)
      || (! classic
        && (
	  (log_ev < *best_log_ev)
	  || (log_ev == *best_log_ev && N > best_nsites)
	  || (log_ev == *best_log_ev && N == best_nsites && score > *best_llr)
        )
      )
    ) {
      best_nsites = N;				/* number of sites */
      *best_log_pv = log_pv;
      *best_log_ev = log_ev;
      *best_llr = score;
    }

  } /* nsites */

  myfree(s_points);

  return best_nsites;
} /* get_best_nsites */

/***********************************************************************/
/*
	set_z

	Set the z to 1/0 using list of sites (maxima).
*/
/***********************************************************************/
void set_z (
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
)
{
  int i, j;
  int nsites = model->nsites_dis;		/* new nsites */
  P_PROB maxima = model->maxima;		/* the maxima positions */
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */
  BOOLEAN invcomp = model->invcomp;

  /* set all z to 0 */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sample i */
    int lseq = s->length;			/* length of sequence */
    int min_j = invcomp ? -lseq : 0;		// minimum Z_i
    int max_j = lseq;				// maximum Z_i
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    for (j=min_j; j<=max_j; j++) {		// Z_i = j
      Zi(j) = 0;
    }
  }

  /* set z 1 for selected sites */
  for (i=0; i<nsites; i++) {
    SAMPLE *s = samples[maxima[i].x];		/* sample */
    if (maxima[i].negative){printf("NSITE %d is a negative !\n", i); exit(1);}
    int y = maxima[i].y;			/* position of site */
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    int j = maxima[i].ic ? -(y+1) : y+1;	// value of Z_i
    Zi(j) = 1.0;
  }

} /* set_z */

/***********************************************************************/
/*
	set_pY

	Initialize pY from z with given offset.
*/
/***********************************************************************/
static void set_pY(
  int w, 				/* motif width */
  BOOLEAN invcomp, 			/* use reverse complement strand, too */
  BOOLEAN pal,				/* force palindrome */
  DATASET *dataset,			/* the dataset */
  int off				/* offset to shift motif */
)
{
  int i, j;
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */

  /* 
    put integerized, weighted log z into pY array
  */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sequence i */
    int lseq = s->length;			/* sequence length */
    int last_j = lseq-w;			/* last start */
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    double sw = s->sw;				/* weight of sequence */
    int *pY = s->pY[0];				/* p(Y_j | theta_1) both */
    char *pYic = s->pYic;			/* site on - strand */

    if (lseq < w) continue;			/* sequence too short */

    /* initialize pY and pYic */
    for (j=0; j<=last_j; j++) {
      pY[j] = INT_LOG(0.0);			/* z == 0 */
      pYic[j] = '\0';				/* site on + strand */
    }
    for (j=0; j<lseq; j++) {			/* site start */
      int jj = j + off;				/* new site start */
      int k = jj+1;				// Z_i = k

      if (jj<0 || jj>last_j) continue;

      /* no z available? */
      if (j > last_j) {				/* no z available */
        pY[jj] = 0;				/* no site */
        pYic[jj] = '\0';			/* strand doesn't matter */
        continue;
      } 

      /* not using inverse strand, too? */
      if (!invcomp) {
        pY[jj] = INT_LOG(sw * Zi(k));
        pYic[jj] = '\0';			/* site on + strand */
        continue;
      } 

      /* using inverse complement strand, too */
      if (pal) {				// use sum of Zi(-k)+Zi(k)
        pY[jj] = INT_LOG(sw * MIN(1.0,(Zi(-k)+Zi(k))));	// FIXME??
      } else if (Zi(-k) > Zi(k)) {		// - strand
        pY[jj] = INT_LOG(sw * Zi(-k));
      } else {					// + strand
        pY[jj] = INT_LOG(sw * Zi(k));
      }
      pYic[jj] = (Zi(-k) > Zi(k)) ? '\1' : '\0';	/* choose strand */

    } /* site start */

  } /* sequence */

} /* set_pY */

/***********************************************************************/
/*
	maxima_compare

	Compare maxima based on "prob".
        Return >0 if second has the higher probability, then
	>0 if second is a control and first is a positive, then
	>0 if second has larger X, then
	>0 if second has larger Y.
*/
/***********************************************************************/
static int maxima_compare(
  const void *v1,
  const void *v2
)
{
  const struct p_prob *s1 = (const struct p_prob *) v1;
  const struct p_prob *s2 = (const struct p_prob *) v2;

  if (s1->prob != s2->prob) {
    return(s2->prob > s1->prob ? +1 : -1);
  } else if (s1->negative != s2->negative) {
    return(s2->negative ? +1 : -1);
  } else if (s1->x != s2->x) {
    return(s2->x - s1->x);
  } else {
    return(s2->y - s1->y);
  }
} // maxima_compare

/***********************************************************************/
/*
	sort_maxima

	In order to sort the maxima, we assign their "prob"
	field either Z (classic) or the LLR of the site (otherwise).
  	Break ties to sort negatives first.  That way the last
  	positive tied with them will yield the best p-value.
  	This should make p-values (slightly) conservative if nsites is limited.

        Sets the ranks of the positive maxima and removes the control 
	maxima if there are any.
*/
/***********************************************************************/
static void sort_maxima(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the primary dataset */
  DATASET *control,			/* the control dataset */
  int n_pos_maxima,			/* number of possible sites */
  int n_neg_maxima			/* number of possible control sites */
)
{
  int i, j;
  BOOLEAN invcomp = model->invcomp;     	/* reverse complement strand */
  int n_maxima = n_pos_maxima + n_neg_maxima;	// total possible real/control sites

  //
  // Set the "prob" field in the maxima.  Use:
  //   objfun == Classic: z_i of site
  //   otherwise: erased LLR of the site
  //
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    BOOLEAN neg = model->maxima[i].negative;
    SAMPLE *s = neg ? control->samples[x] : dataset->samples[x];
    int y = model->maxima[i].y;
    double score;
//FIXME: Classic might be improved if it used erased LLR, too.
// TLB: I left this alone to keep Classic identical with v4.10.2
// if sequence order shuffling is not done.
    if (dataset->objfun == Classic) {
      double *zi = s->z;			// Zi[j], j in [-lseq...+lseq]
      int j = y + 1;				// Z_i = j
      score = invcomp ? MIN(1.0,Zi(-j)+Zi(j)) : Zi(j);
    } else {
      // Use the (erased) LLR of sites for sorting maxima.
      int w = model->w;
      double *lcb = s->logcumback;		/* log cumulative bkg. probability */
      THETA logtheta1 = model->logtheta;	/* motif log(theta) */
      double *not_o = s->not_o;			/* Pr(V_ij = 1) */
      double init = -Log_back(lcb, y, w) + LOG(not_o[y]);
      double llr = init;
      for (j=0; j<w; j++) llr += logtheta1(j, (int) s->res[y+j]);
      // negative strand
      if (invcomp) {
        double llr_rc = init;
        THETA logtheta1_rc = model->logtheta_rc;/* motif log(theta) rev. comp. */
        for (j=0; j<w; j++) llr_rc += logtheta1_rc(j, (int) s->res[y+j]);
        llr = MAX(llr, llr_rc);
      }
      score = llr;
    }
    RND(score, 11, model->maxima[i].prob);
  }

  //
  // Sort the maxima by "prob" which contains Z or LLR of maxima.
  // Break ties to sort negatives first.  That way the last
  // positive tied with them will have the best p-value.
  // This should make p-values (slightly) conservative if nsites is limited.
  //
  qsort((char *) model->maxima, n_maxima, sizeof(p_prob), maxima_compare);

  //
  // Set the ranks of the positive maxima and remove the control maxima if any
  //
  if (n_neg_maxima > 0) {
    for (i=j=0; i<n_maxima; i++) {
//FIXME
//printf("maxima: i %d %s %s %d score %g\n", i, model->maxima[i].negative?"neg":"pos", dataset->samples[model->maxima[i].x]->sample_name, model->maxima[i].y, model->maxima[i].prob);
      if (! model->maxima[i].negative) {
        int rank = i+1;
	model->maxima[i].rank = rank;
	model->maxima[i].ranksum = (j==0) ? rank : model->maxima[j-1].ranksum + rank;
	model->maxima[j++] = model->maxima[i];
      }
    }
  }

  if (debug) print_maxima(n_pos_maxima, model);

} // sort_maxima


/***********************************************************************/
/*
	get_column_re_list

  Return a list of (col_num, col_re) ordered by col_re,
  where col_re is the relative entropy of the motif column.
*/
/***********************************************************************/
  PAIR *get_column_re_list(
    MODEL *model, 
    DATASET *dataset
)
{
  int i;

  // Get the relative entropy (RE) of each motif column.
  calc_entropy(model, dataset);
  PAIR *pairs = NULL;
  Resize(pairs, model->w, PAIR);

  // Set up the index/RE pairs for sorting.
  for (i=0; i< model->w; i++) {
    pairs[i].ind = i;
    pairs[i].re = model->rentropy[i];
  }

  // Sort column indices by increasing column RE.
  qsort((void *)pairs, model->w, sizeof(PAIR), (void *)pairCompare);

  return(pairs);

} // get_column_re_list

/***********************************************************************/
/*
	diff_set_pY

	Determine the set of columns to mask, masking
	some number of lowest relative entropy columns
	to the uniform distribution.
	Shift the motif (and shorten it) so that 
	there are no masked columns on the flanks.
	Do an E_STEP on both the primary and control datasets.
	Initialize pY from Z for both datasets.

	Sets model->theta to the selected range of columns
	with masked internal columns containing their original values.
	Sets model->w to the new width.

	Returns the best offset of the original motif.
*/
/***********************************************************************/
static int diff_set_pY(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, j;
  int min_w = model->min_w;			/* minimum width */
  BOOLEAN invcomp = model->invcomp;     	/* reverse complement strand */
  BOOLEAN pal = model->pal;     		/* force DNA palindromes */
  ALPH_T *alph = dataset->alph;			// alphabet

  // Save a copy of the motif.
  THETA theta = model->theta;			//theta of motif
  if (saved_theta == NULL) {
    create_2array(saved_theta, double, model->max_w+1, alph_size_wild(alph));
  }
  copy_theta(theta, saved_theta, model->w, alph_size_wild(alph));

  //
  // Get a list of (col_num, col_re) ordered by col_re,
  // where col_re is the relative entropy of the motif column.
  //
  PAIR *pairs = get_column_re_list(model, dataset);

  //
  // Decide how many hi-RE columns to keep and hence the motif width.
  //
  int n_keep;
  // Tanaka 2014 code
  n_keep = 6 + floor((model->w-6)/3);
  //
  //n_keep = ceil(2 * sqrt(model->w));
  //n_keep = min_w + ceil(sqrt(model->w));
  //n_keep = ceil(sqrt(model->w));
  //n_keep = min_w + 2*ceil(sqrt(model->w));
  //n_keep = min_w + ceil(sqrt(model->w - min_w));
  //
  // Progression code:
  // Use sqrt(2) progression just like start widths in starts.c.
  // n_keep = min_w + floor((model->w - min_w)/3);
  // int prev_w = (int) (model->w/sqrt(2) + 0.5);	// prev w in progression
  // Make sure all smaller widths are reachable.
  // if (n_keep > prev_w + 1) n_keep = prev_w + 1;
  //
  // Make sure not keeping too many.
  if (n_keep > model->w) n_keep = model->w;
  // Make sure not keeping too few.
  if (n_keep < model->min_w) n_keep = model->min_w;	
  // How many columns to mask:
  int n_mask = model->w - n_keep;
  //FIXME
  //printf("n_keep = %d n_mask = %d model->w %d\n", n_keep, n_mask, model->w);

  // Mask non-selected, low-RE columns.
  double f = 1.0/alph_size_core(alph);		// uniform frequency
  for (i=0; i<n_mask; i++) {
    int index = pairs[i].ind;
    for (j=0; j<alph_size_core(alph); j++) theta(index, j) = f;
    theta(index, j) = 0;			// not sure if necessary
  }

  // Get the start and end indices of the non-masked columns.
  int small_index = model->w - 1;
  int big_index = 0;
  for (i=model->w-1; i>=n_mask; i--) {
    int index = pairs[i].ind;
    if (index < small_index) small_index = index;
    if (index > big_index) big_index = index;
  }
  myfree(pairs);

  // Shift model and set new width (contains only non-masked columns).
  if (n_mask > 0) {
    int new_w = big_index - small_index + 1;
    // Increase width of motif to at least min_w, centering.
    if (new_w < min_w) {
      int left_margin = small_index;
      int right_margin = model->w - big_index - 1;
      while (new_w < min_w) {
	// Increase and center.
	if (left_margin >= right_margin) {
	  left_margin = --small_index;
	} else {
	  right_margin = model->w - (++big_index) - 1;
	}
	new_w++;
      }
    }
    if (small_index > 0) copy_theta(model->theta+small_index, model->theta, new_w, alph_size_wild(alph));
    model->w = new_w;
  }

  // E_STEP on the primary dataset.
  get_not_o(dataset, model->w);			// Already set during EM.
  model->ll = E_STEP(model, dataset);		// set z using masked model
  set_pY(model->w, invcomp, pal, dataset, 0);

  // E_STEP on the control dataset.
  get_not_o(neg_dataset, model->w);		// Not set yet so initialize.
  E_STEP(model, neg_dataset);
  set_pY(model->w, invcomp, pal, neg_dataset, 0);

  // Restore the unmasked span of the original model's theta matrix.
  copy_theta(saved_theta+small_index, model->theta, model->w, alph_size_wild(alph));

  return(small_index);
} // diff_set_pY

/***********************************************************************/
/*
	meme_mHG

	Do an E_STEP on both a holdout set.
	Do an E_STEP on the control set.
	Initialize pY from Z for both datasets.
	Find the best value of sites.

	Returns the best number of sites and the log(p-value)
	and LLR of that p-value.
*/
/***********************************************************************/
static int meme_mHG(
  MODEL *model,				/* the model */
  DATASET *primary,			/* the primary dataset */
  DATASET *control,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *),	// E_STEP function
  BOOLEAN primary0,			// include group 0 in primary
  BOOLEAN primary1,			// include group 1 in primary
  BOOLEAN primary2,			// include group 2 in primary
  BOOLEAN control0,			// include group 0 in control
  BOOLEAN control1,			// include group 1 in control
  BOOLEAN control2,			// include group 2 in control
  int min_nsites,			/* minimum nsites */
  int max_nsites,		 	/* maximum nsites */
  double score_thresh,			// restrict to maxima
					// with scores >= score_thresh
					// set to -BIG to ignore
  double *log_pv,			/* log p-value of score */
  double *llr  				// LLR of best p-value
)
{
  int i, j;
  int nsites;
  BOOLEAN invcomp = model->invcomp;     	/* reverse complement strand */
  BOOLEAN pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  double col_scores[MAXSITE];			/* column scores */

  // E_STEP on primary dataset.
  get_not_o(primary, model->w);
  model->ll = E_STEP(model, primary);		// set z using masked model
  set_pY(model->w, invcomp, pal, primary, 0);

  // E_STEP on the control dataset.
  if (control && control != primary) {
    get_not_o(control, model->w);		// Not set yet so initialize.
    E_STEP(model, control);
    set_pY(model->w, invcomp, pal, control, 0);
//FIXME
//for (i=0;i<control->n_samples;i++) {printf("%s ", control->samples[i]->sample_name); for(j=0;j<control->samples[i]->length-model->w+1;j++) printf("%d ", control->samples[i]->pY[0][j]); printf("\n");}
  }

  //
  // get the possible sites (maxima) in the primary and control datasets
  //
  // Primary:
  set_seq_groups_to_include(primary, primary0, primary1, primary2);
  int n_pos_maxima = get_max(mtype, primary, FALSE, model->w, model->maxima, invcomp, FALSE);
  restore_seq_groups_to_include(primary);
  // Control:
  int n_neg_maxima = 0;
  if (control) {
    set_seq_groups_to_include(control, control0, control1, control2);
    n_neg_maxima = get_max(mtype, control, TRUE, model->w, (model->maxima)+n_pos_maxima, invcomp, FALSE);
    restore_seq_groups_to_include(control);
//FIXME
//printf("n_neg_maxima: %d\n", n_neg_maxima);
  }

  //
  // * Sort all maxima (in both datasets together) by site LLR (in Z).
  // 	Break ties to sort negatives first. That way the last
  // 	positive tied with them will yield the best p-value.
  // 	This should make p-values (slightly) conservative if nsites is limited.
  // * Rank the positive maxima.
  // * Remove the control maxima.
  //
  sort_maxima(model, primary, control, n_pos_maxima, n_neg_maxima);

  //
  // Find the value of nsites corresponding to the score threshold.
  // (Don't do this for OOPS model).
  //
  if (mtype != Oops && score_thresh != -BIG) {
    // Binary search of maxima for maxima[i-1].prob < score_threshold <= maxima[i].prob
    int lo = 0;
    int hi = n_pos_maxima-1;
    while (hi-lo > 1) {			// binary search
      int mid = (lo+hi)/2;		// midpoint
      if (model->maxima[mid].prob >= score_thresh) { lo = mid; } else { hi = mid; }
    }
    nsites = model->maxima[hi].prob >= score_thresh ? hi+1 : lo+1;
    max_nsites = MIN(max_nsites, nsites);
    max_nsites = MAX(min_nsites, max_nsites);
    min_nsites = max_nsites;
  }

  //
  // Get the best number of sites (if score threshold not given) and get the p-value.
  //
  nsites = get_best_nsites(model, primary, min_nsites, max_nsites,
      model->w, n_pos_maxima, n_neg_maxima, model->maxima, col_scores, log_pv, log_pv, llr);
  return(nsites);
} // meme_mHG

/***********************************************************************/
/*
	classic_get_width_and_nsites

	Adjust the width of the motif by
		1) using the multiple alignment trim procedure	
		2) optimizing E-value over all subsets of columns

	Updates the best motif information into the model.

	Returns the best starting offset.
*/
/***********************************************************************/
static int classic_get_width_and_nsites(
  MODEL *model,				/* the model */
  DATASET *dataset, 			/* the dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, w, ini_w, best_w, min_w, max_w, ma_w, ma_off;
  int n_maxima; 				/* number of possible sites */
  int min_nsites = model->min_nsites;		/* minimum nsites */
  int max_nsites = model->max_nsites;		/* maximum nsites */
  BOOLEAN invcomp = model->invcomp;     	/* reverse complement strand */
  BOOLEAN pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  double col_scores[MAXSITE];			/* column scores */
  double score;					// score of motif (-log_pop or LLR)
  double log_pv, best_log_pv;			/* log p-value of motif */
  double log_ev;				/* log E-value of score */
  double llr;					// LLR of best E-value

  //
  // initialize pY from Z
  //
  set_pY(model->w, invcomp, pal, dataset, 0);

  // Don't include the PSP probabilities in log_not_o (used in get_max())
  // because z already takes them into account.
  get_not_o(dataset, model->w);

  //
  // Get the maxima.
  //
  n_maxima = get_max(mtype, dataset, FALSE, model->w, model->maxima, invcomp, FALSE);

  if (n_maxima < min_nsites) return 0; 		// Exit if not enough sites.

  //
  // Sort the maxima.
  //
  sort_maxima(model, dataset, NULL, n_maxima, 0);

  //
  // Set the minimum and maximum widths.
  //
  ini_w = best_w = model->w;			/* initial motif width */
  min_w = model->min_w;				/* minimum width */

  //
  // Get the best nsites using the full motif width.
  // 
  int best_nsites = (min_nsites==max_nsites) ? min_nsites :
    get_best_nsites(model, dataset, min_nsites, max_nsites,
      model->w, n_maxima, 0, model->maxima, col_scores, &log_ev, &log_ev, &llr);

  //
  // Trim the alignment to include the minimum possible number of gaps
  //
  max_w = ma_w = ini_w;
  ma_off = 0;
  w = 0;
  if (dataset->ma_adj) {
    int flank = w/2;				/* amt. of flank for mult. a. */
    ma_w = ma_adjust(model->maxima, best_nsites, max_w, flank, min_w, dataset, &ma_off);
    max_w = ma_w;				/* new maximum width */
    /*
      update the maxima positions by shifting them
    */
    for (i=0; i<best_nsites; i++) {
      BOOLEAN ic = model->maxima[i].ic;			/* on - strand */
      model->maxima[i].y += (ic ? ini_w-ma_w-ma_off : ma_off);
    }
  } /* ma_adj */

  //
  // Shorten the motif based on p-value of alignments.
  //
  /* 
    get the p-values of the columns given the best number of sites
    and the gap-trimmed width
  */
  (void) get_best_nsites(model, dataset, best_nsites, best_nsites,
    max_w, n_maxima, 0, model->maxima, col_scores, &log_ev, &log_ev, &llr);

  /* 
    find subsequence of columns with best p-value
  */
  best_log_pv = log_ev = BIG;
  int best_off = 0;
  for (w=min_w; w<=max_w; w++) {		/* width */
    int l_off;					/* left offset */

    for (l_off=0; l_off<=max_w-w; l_off++) {	/* left edge */

      /* get the product of column p-values */
      for (i=score=0; i<w; i++) score += col_scores[i+l_off];

      /* get the p-value of the pop */
      if (dataset->objfun == Classic) {		// score is log(product of p-values)
	log_pv = get_log_sig(-score, mtype, w, best_nsites, 0, invcomp, pal, dataset);
      } else if (dataset->objfun == NC) {	// score is LLR
	log_pv = get_log_sig(score, mtype, w, best_nsites, 0, invcomp, pal, dataset);
      } else {
	fprintf(stderr, "\nUnknown type of objective function.\n");
	exit(1);
      }

      if (TRACE) 
	printf(
	"ini_w %d ma_w %d w %d ma_off %d off %d log_pv %f init cons %*.*s\n", 
	ini_w, ma_w, w, ma_off, l_off, log_pv, w, w, model->cons0+l_off+ma_off);

      /* save if best so far: better log_pv */
      if (RNDEPS < best_log_pv - log_pv) {
	if (TRACE) printf("better: w %d best log_pv %g\n", w, log_pv);
	best_w = w;
	best_off = l_off;
	best_log_pv = log_pv;
      }
    } /* l_off */
  } /* w  */

  /*
    update the maxima positions by shifting them
  */
  for (i=0; i<best_nsites; i++) {
    BOOLEAN ic = model->maxima[i].ic;			/* on - strand */
    model->maxima[i].y += (ic ? ma_w-best_w-best_off : best_off);
  }
  
  /* 
    get the best number of sites for the shortened motif and the final E-value 
  */
  best_nsites =
    get_best_nsites(model, dataset, min_nsites, best_nsites,
      best_w, n_maxima, 0, model->maxima, col_scores, &log_ev, &log_ev, &llr);

  /* 
    set the best motif info in the model
  */
  model->w = best_w;				/* best width */
  model->nsites_dis = best_nsites;		/* after discretization */
  model->logpv = best_log_pv;			/* p-value */
  model->logev = log_ev;			/* E-value */
  model->llr = llr;				// LLR of best E-value

  return(best_off);
} // classic_get_width_and_nsites

/***********************************************************************/
/*
	smhg_get_width_and_nsites

	Adjust the width of the motif by using a number
	of high relative entropy columns determined based
	on motif width.
	Determine the optimal nsites using the selective mHG test,
	which masks non-selected columns within the window
	spanned by the selected high-RE columns.

	Updates the best motif information into the model.

	Returns the best starting offset.
*/
/***********************************************************************/
static int smhg_get_width_and_nsites(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int min_nsites = model->min_nsites;		/* minimum nsites */
  int max_nsites = model->max_nsites;		/* maximum nsites */
  double log_pv;				/* log p-value of score */
  double llr;					// LLR of best E-value

  //
  // Initialize pY (and set best motif width in model).
  //
  int offset = diff_set_pY(model, dataset, neg_dataset, E_STEP);

  //
  // Use the mHG test to determine the best number of sites.
  //
  int nsites = meme_mHG(model, dataset, neg_dataset, E_STEP, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    min_nsites, max_nsites, -BIG, &log_pv, &llr);

  //
  // Update the model.
  //
  model->nsites_dis = (min_nsites==max_nsites) ? min_nsites : nsites;
  model->logpv = log_pv;			/* p-value */
  model->logev = log_pv;			/* E-value */
  model->llr = llr;				// LLR of best E-value

  return(offset);
} // smhg_get_width_and_nsites

/***********************************************************************/
/*
	nz_cv_ll_get_width_and_nsites

	Find the optimal number of sites and widths for the motif by
	optimizing:
		NZ: mHG on groups 0,2 vs. 1
		CV: mHG on group 1 vs. control
		LL: LLR on primary sequences (groups 1,2,3)
	Note: In OOPS mode, the Wilcoxon rank-sum test replaces
	mHG.

	Widths are searched by using the motif columns
	spanned by the first 2, 3, ..., w highest-RE columns.
	The optimal score threshold is the one that gives the
	best value of the objective function (mHG, or LLR).
	Sets the best score (mHG or LLR) in model->logev, 
        sets the best LLR in model->llr, and the best nsites
	in model->nsites_dis. These are to be used for optimization 
	over starting points.

	Find the final p-value for the motif using the optimal
	score threshold determined above and:
		NZ: Fisher exact test on group 0 vs. 2 
		CV: Fisher exact test on group 2 vs. control
		LL: Fisher exact test on primary vs. control
	Sets that p-value in model->logpv.

	Returns the best starting offset.
*/
/***********************************************************************/
static int nz_cv_ll_get_width_and_nsites(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset (including noise) */
  DATASET *neg_dataset,			// the control dataset
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, w;
  int ini_w = model->w;				// initial motif width
  DATASET *control = NULL;
  OBJTYPE objfun = dataset->objfun;		// objective function
  BOOLEAN pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  int min_nsites = model->min_nsites;		/* minimum nsites */
  int max_nsites = model->max_nsites;		/* maximum nsites */
  ALPH_T *alph = dataset->alph;			// alphabet
  double log_pv;				/* log p-value of score */
  double llr;					// LLR of best p-value
  int nsites;
  BOOLEAN p0=0, p1=0, p2=0, c0=0, c1=0, c2=0;	// skip flags for primary/control
  double thresh;				// score threshold
  double scale=0;				// for adjusing min and max nsites

  // Save a copy of the motif.
  THETA theta = model->theta;			// theta of motif
  if (saved_theta == NULL) {
    create_2array(saved_theta, double, model->max_w+1, alph_size_wild(alph));
  }
  copy_theta(theta, saved_theta, model->w, alph_size_wild(alph));

  //
  // Get a list of (col_num, col_re) ordered by increasing col_re,
  // where col_re is the relative entropy of the motif column.
  //
  PAIR *pairs = get_column_re_list(model, dataset);

  //
  // Find the best motif using mHG applied to each span of
  // high-RE columns and the holdout1 + control sets.
  //
  int small_index = model->w - 1;
  int big_index = 0;
  int best_start = 0;
  int best_end = 0;
  int best_nsites = 0;
  double best_log_pv = BIG;
  double best_llr = -BIG;
  double best_score = -BIG;
  int best_w = model->w;
  int maxsites=0, minsites=0;

  //
  // Add motif columns in order of decreasing relative entropy.
  // Stop when they span the maximum width.
  //
  for (i=model->w-1; i>=0; i--) {
    int index = pairs[i].ind;
    if (index < small_index) small_index = index;
    if (index > big_index) big_index = index;
    int new_w = big_index - small_index + 1;

    // Motif width must be >= mininimum width
    if (new_w < model->min_w) continue;

    // Replace the motif with the new span of high-RE columns.
    if (small_index > 0) copy_theta(saved_theta+small_index, model->theta, new_w, alph_size_wild(alph));
    model->w = new_w;

    //
    // Set the parameters for scoring this motif.
    //
    if (objfun == NZ) {
      // Set parameters to score motif using mHG on primary group 0,2 vs primary group 1.
      control = dataset;
      p0=TRUE; p1=FALSE; p2=TRUE;
      c0=FALSE; c1=TRUE; c2=FALSE;
      thresh = -BIG;
      scale = (double) (dataset->n_group[0] + dataset->n_group[2])/dataset->n_group[0];
    } else if (objfun == CV) {
      // Set parameters to score motif using mHG on group 1 (holdout 1) vs. control.
      control = neg_dataset;
      p0=FALSE; p1=TRUE; p2=FALSE;
      c0=TRUE; c1=TRUE; c2=TRUE; 
      thresh = -BIG;
      scale = (double) dataset->n_group[1]/dataset->n_samples;
    } else if (objfun == LL) {
      // Set parameters to score motif by computing its optimum LLR on group 0.
      control = NULL;
      p0=TRUE; p1=FALSE; p2=FALSE;
      c0=FALSE; c1=FALSE; c2=FALSE;
      thresh = -BIG;
      scale = (double) dataset->n_group[0]/dataset->n_samples;
    }

    //
    // Score motif using the parameters set above.
    //
    minsites = MAX(2,scale*min_nsites);
    maxsites = MAX(2,scale*max_nsites);
    nsites = meme_mHG(model, dataset, control, E_STEP, p0, p1, p2, c0, c1, c2, minsites, maxsites, thresh, &log_pv, &llr);

    //
    // Save motif offset, nsites and p-value information if best found so far.
    // Break ties using 1) width 2) LLR
    //
    if (
      (log_pv < best_log_pv)
      || (log_pv == best_log_pv && new_w < best_w) 
      || (log_pv == best_log_pv && new_w == best_w && llr > best_llr)
    ) {
      best_start = small_index;
      best_end = big_index;
      best_w = best_end - best_start + 1; 
      best_nsites = nsites;
      best_log_pv = log_pv;
      best_llr = llr;
      best_score = model->maxima[nsites-1].prob;
    }

    // Done if span is same as initial motif width;
    if (new_w == ini_w) break;
  }
  myfree(pairs);

  //
  // Extract the best columns from the motif and re-set its width.
  //
  model->w = best_w;
  if (best_start > 0) copy_theta(saved_theta+best_start, model->theta, best_w, alph_size_wild(alph));

  //
  // Store the SCORE and LLR of the best motif found above.
  //
  model->logev = best_log_pv;			// SCORE for comparing models
  model->llr = best_llr;			// LLR of best SCORE

  //
  // Set up parameters for estimating the final p-value of the motif.
  //
  if (objfun == NZ) {
    // Set up parameters to use the Fisher exact test with best_score threshold on group 0 vs. 2.
    control = dataset;
    p0=TRUE; p1=FALSE; p2=FALSE;
    c0=FALSE;c1=FALSE; c2=TRUE;
    //FIXME: test1 on crp0.s with 20 runs with seeds [0,..,19] showed using best_score rather than best_nsites has 15 wins, 3 ties, 2 losses.
    thresh = best_score;
    scale = 1;
//FIXME:
if (debug) printf("Calling mHG with min_nsites %d max_nsites %d thresh %g p0 %d p1 %d p2 %d c0 %d c1 %d c2 %d\n", minsites, maxsites, thresh, p0, p1, p2, c0, c1, c2);
  } else if (objfun == CV) {
    // Set up parameters to use the Fisher exact test on group 2 (holdout2) vs. control.
    control = neg_dataset;
    p0=FALSE; p1=FALSE; p2=TRUE;
    c0=TRUE; c1=TRUE; c2=TRUE;
    scale = (double) dataset->n_group[2]/dataset->n_samples;
    thresh = best_score;
  } else if (objfun == LL) {
    // Use mHG on group 1 vs. control to get best_score.
    // Then set up parameters to use mHG on group 2 vs. control.
    // 1) Get best_score threshold using mHG (group 1 vs. control):
    control = neg_dataset;
    p0=FALSE; p1=TRUE; p2=FALSE;
    c0=TRUE; c1=TRUE; c2=TRUE;
    thresh = -BIG;
    scale = (double) dataset->n_group[1]/dataset->n_samples;
    minsites = MAX(2, scale*min_nsites);
    maxsites = MAX(2, scale*max_nsites);
    nsites = meme_mHG(model, dataset, control, E_STEP, p0, p1, p2, c0, c1, c2, minsites, maxsites, thresh, &log_pv, &llr);
    best_score = model->maxima[nsites-1].prob;
    // 2) Set up to get final p-value using Fisher exact test (group 2 vs. control):
    control = neg_dataset;
    p0=FALSE; p1=TRUE; p2=FALSE;
    p0=FALSE; p1=FALSE; p2=TRUE;
    c0=TRUE; c1=TRUE; c2=TRUE;
    thresh = best_score;
    scale = (double) dataset->n_group[2]/dataset->n_samples;
  }

  //
  // Get the final p-value using type of test determined
  // by the values set above.
  //
  minsites = MAX(2, scale*min_nsites);
  maxsites = MAX(2, scale*max_nsites);
  nsites = meme_mHG(model, dataset, control, E_STEP, p0, p1, p2, c0, c1, c2, minsites, maxsites, thresh, &log_pv, &llr);
  model->logpv = log_pv;			/* FINAL p-value */

//FIXME
if (debug) printf("FINAL p-value: ini_w %d w %d offset %d best_nsites %d new_nsites %d final logpv %g final pv %g final llr %g selection pv %g\n", ini_w, model->w, best_start, best_nsites, nsites, log_pv, exp(log_pv), llr, exp(model->logev));

  //
  // Set up parameters to get the optimum value of nsites
  // and to set the maxima for use by em().
  //
  if (objfun == NZ) {
    // Set up to use mHG on group 0 vs. groups 1,2. 
    control = dataset;
    p0=TRUE; p1=FALSE; p2=FALSE;
    c0=FALSE; c1=TRUE; c2=TRUE;
    thresh = -BIG;
  } else if (objfun == CV || objfun == LL) {
    // Set up to use mHG on group 0,1,2 vs. control.
    control = neg_dataset;
    p0=TRUE; p1=TRUE ; p2=TRUE;
    c0=TRUE; c1=TRUE; c2=TRUE;
    thresh = -BIG;
  }

  //
  // Get the final nsites and set the maxima using the parameters set above.
  //
  minsites = min_nsites;
  maxsites = max_nsites;
  model->nsites_dis = meme_mHG(model, dataset, control, E_STEP, p0, p1, p2, c0, c1, c2, minsites, maxsites, thresh, &log_pv, &llr);
//FIXME
if (debug) printf("TUNED NSITES using ALL: ini_w %3d w %3d offset %4d selection score %g final pv %1g full score %g new llr %g selection nsites %3d full nsites %3d\n", ini_w, model->w, best_start, model->logev, exp(model->logpv), exp(log_pv), llr, best_nsites, model->nsites_dis);

  // Return offset of final motif from start of original motif.
  return(best_start);
} // nz_cv_ll_get_width_and_nsites

/***********************************************************************/
/*
	discretize	

	Search over width and offset of motif to minimize E-value of
	log likelihood ratio.  

		1) get best nsites using E-value
		2) calculate p-value of each column of motif
		3) shorten using p-value
		4) get best nsites using E-value

	Sets z to 1.0 for sites, 0 for non-sites in dataset.
	Sets w, nsites_dis, maxima and sig in model.

	Returns the optimum number of sites.
*/
/***********************************************************************/
void discretize(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  DATASET *neg_dataset,			/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *)	// E_STEP function
)
{
  int i, j;
  MOTYPE mtype = model->mtype;			/* type of model */
  OBJTYPE objfun = dataset->objfun;		// objective function
  BOOLEAN classic = (objfun == Classic);	// Classic meme

  /* 
    set initial and minimum allowed widths
  */
  int ini_w = model->w;				/* initial motif width */
  int min_w = model->min_w;			/* minimum width */

  /* 
    create space for the maxima from both datasets
  */
  int n_maxima = (mtype==Tcm) ? ps(dataset, min_w) : dataset->n_samples;
  if (objfun != Classic && objfun != NZ) {
    n_maxima += (mtype==Tcm) ? ps(neg_dataset, min_w) : neg_dataset->n_samples;
  }
  Resize(model->maxima, n_maxima, p_prob);

  /* 
    Get best width and nsites for motif.
  */
  int best_off = 0;
  if (classic || objfun==NC) {
    best_off = classic_get_width_and_nsites(model, dataset, E_STEP);
  } else if (objfun==SmHG) {
    best_off = smhg_get_width_and_nsites(model, dataset, neg_dataset, E_STEP);
  } else if (objfun==CV) {
    best_off = nz_cv_ll_get_width_and_nsites(model, dataset, neg_dataset, E_STEP);
  } else if (objfun==NZ) {
    best_off = nz_cv_ll_get_width_and_nsites(model, dataset, NULL, E_STEP);
  } else if (objfun==LL) {
    best_off = nz_cv_ll_get_width_and_nsites(model, dataset, neg_dataset, E_STEP);
  } else {
    fprintf(stderr, "\nUnknown type of objective function.\n");
    exit(1);
  }
  
//FIXME
    if (!NO_STATUS) {
#if defined(PARALLEL) && defined(DEBUG_PARALLEL)
      fprintf(stdout, "\nnode %d FINAL: ini_w %d w %d nsites %d selection %.17g pv %g llr %g\n", mpMyID(), ini_w, model->w, model->nsites_dis, exp(model->logev), exp(model->logpv), model->llr);
#else
      //printf("\nFINAL: ini_w %d w %d nsites %d selection %g pv %g llr %g\n", ini_w, model->w, model->nsites_dis, exp(model->logev), exp(model->logpv), model->llr);
#endif
    }

  //
  // Discretize the sites in Z from the maxima stored in the model.
  //
  set_z(model, dataset);

  if (TRACE) {
    int w = model->w;
    printf( 
      "ini_w %d w %d off %d nsites %d pv %9.2e EV %9.2e cons %*.*s %s\n",
      ini_w, w, best_off, model->nsites_dis, exp(model->logpv), 
      exp(model->logev), w, w, (model->cons0)+best_off, model->cons0);
  }

} /* discretize */

/**********************************************************************/
/*
	ma_adj

	Shorten a motif to the longest g-alignment of width at least
	min_w.  A g-alignment is an alignment with no more than g
	gapped sequences per column.  Values of g in [0..] are tried
	until an aligment of width min_w or greater is found.

	Returns best width and offset.
*/
/**********************************************************************/
static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */ 
  int w,                                        /* width of sites */ 
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off					/* best offset */
)
{
  char **ma;					/* multiple aligment */
  int left = MIN(flank, sites[0].y);		/* left edge after algnmnt. */
  int right = left+w-1;				/* right edge after algnmnt. */

  /* get the multiple alignment */
  ma = dp_multi_align(sites, nsites, w, flank, dataset);

  /* get longest g-alignment of width at least min_w */
  (void) g_align(ma, nsites, strlen(ma[0]), left, right, min_w, off, &w);

  /* free space */
  free_2array(ma, nsites);

  return w;					/* return the width */
} /* ma_adjust */

