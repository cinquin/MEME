/*
 * $Id$
 * 
 * $Log$
 * Revision 1.1.1.1.2.1  2006/01/31 08:34:53  tbailey
 * Put FIXMEs where we branching search is messing up parallel MEME.
 *
 * Revision 1.1.1.1  2005/07/29 17:23:03  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1996, The Regents of the University of California    *
*	Author: Tim Bailey & Bill Grundy			        *
*								       *
***********************************************************************/

#ifdef PARALLEL
 
#include "macros.h"
#include <mpi.h>
#include "mp.h"
#include "message.h"

/**********************************************************************
 *
 * void balance_loop1
 *
 * This function does load-balancing on the search for starting
 * points. It divides the dataset by the number of processors at the
 * character level. The work assigned to each processor is stored in a
 * global struct.
 *
 * IN: SAMPLE **samples             the dataset
 *     int n_samples                size of the dataset
 * OUT: SEQ_PARAMS *start_n_end     (global) start and end points for this node
 *
 * This function is called in 'init.c'.
 **********************************************************************/
void balance_loop1(SAMPLE **samples, int n_samples)
{
  int n_residues;
  float residues_per_node, target;
  int iseq, inode, seq_start, seq_end;

  /* Allocate storage for the sequence parameters. */
  start_n_end = (SEQ_PARAMS *)mymalloc(sizeof(SEQ_PARAMS));

  /* First count how many residues we have. */
  for (iseq = 0, n_residues = 0; iseq < n_samples; iseq++)
    n_residues += samples[iseq]->length;
  
  /* Calculate the number of residues per node. */
  residues_per_node = MAX(1, n_residues / mpNodes());
//FIXME
  if (!NO_STATUS) fprintf(stderr, "BALANCE: samples %d residues %d nodes %d residues/node %.0f\n", n_samples, n_residues, mpNodes(), residues_per_node);
  
  /*
  fprintf(stderr, "n_residues=%d residues_per_node=%g\n", n_residues,
    residues_per_node);
  */
  
  /* Set the starting points for node 0. */
  if (mpMyID() == 0) {
    start_n_end->start_seq = 0;
    start_n_end->start_off = 0;
  }
  
  /* Calculate per-node starting and ending points. */
  for (iseq = 0, seq_start = 0, seq_end = samples[0]->length - 1,
       inode = 0, target = residues_per_node;
       (target < n_residues) || (iseq < n_samples);
       /* Incrementing is done in the loop */) {
    
    /* See whether we've moved past the target. */      
    if (seq_end >= target) {
      if (mpMyID() == inode) {
	start_n_end->end_seq = iseq;
	start_n_end->end_off = (int)target - seq_start - 1;
      }
      if (mpMyID() == inode + 1) {
	start_n_end->start_seq = iseq;
	start_n_end->start_off = (int)target - seq_start;
      }

      /*
      fprintf(stderr, "end %d %d\nnode %d start %d %d ",
	iseq, (int)target-seq_start-1, inode+1, iseq, (int)target-seq_start);
      */

      /* Move to a new node. */
      inode++;
      /* Calculate the target for this node. */
      target = residues_per_node * (inode+1);

    } else {
      
      /* Move to a new sequence. */
      if (++iseq < n_samples) {
	/* Calculate the starting and ending points of this sequence. */
	seq_start = seq_end;
	seq_end += samples[iseq]->length; /* 1 more than the index. */
      }
    }
  }
  
  /* Set the ending points for the last node. */
  if (mpMyID() == mpNodes() - 1) {
    start_n_end->end_seq = n_samples - 1;
    start_n_end->end_off = samples[n_samples - 1]->length - 1;
  }
  
  /* Search for (and correct) fencepost problems. */
  if (start_n_end->start_off >= samples[start_n_end->start_seq]->length) {
    start_n_end->start_seq++;
    start_n_end->start_off = 0;
  } else if (start_n_end->start_off <= -1) {
    start_n_end->start_seq--;
    start_n_end->start_off = samples[start_n_end->start_seq]->length - 1;
  }
  
  if (start_n_end->end_off >= samples[start_n_end->end_seq]->length) {
    start_n_end->end_seq++;
    start_n_end->end_off = 0;
  } else if (start_n_end->end_off <= -1) {
    start_n_end->end_seq--;
    start_n_end->end_off = samples[start_n_end->end_seq]->length - 1;
  }
  
#ifdef DEBUG
  fprintf(stderr, "%d: (%d,%d) -> (%d,%d)\n", mpMyID(), 
	  start_n_end->start_seq, start_n_end->start_off,
	  start_n_end->end_seq, start_n_end->end_off);
#endif /* DEBUG */
  
  /* Do some error checking. */
  /*
     if (inode != mpNodes()) {
       fprintf(stderr, "Warning! inode(%d) != mpNodes(%d)\n",
	     inode, mpNodes());
     }
  */
  if (iseq != n_samples)
    fprintf(stderr, "Warning! iseq(%d) != n_samples(%d)\n",
	    iseq, n_samples);
}

/**********************************************************************
 *
 * void store_consensus
 *
 * Because the reduction across models does not include the
 * candidates, we must transfer the human-readable consensus string
 * into the model so that it can be printed later.
 *
 * This function is called by 'meme.c'.
 **********************************************************************/
void store_consensus(MODEL *model, CANDIDATE *candidates)
{
    int w = model->w;				/* width of motif */
    S_POINT *s_point = candidates[w].s_point; 	/* starting point for model */

    /* Make sure we have at least one candidate. */
    if (s_point == NULL) {
      /* Guarantee that this model won't get chosen. */
      model->logev = log(BIG);
    } else {
      /*fprintf(stderr, "%d: BEFORE consensus=%s\n", mpMyID(),
	s_point->cons0);*/
      strcpy(model->cons0, s_point->cons0);
    }
} /* store_consensus */

/**********************************************************************
 *
 * void save_theta_ptrs
 *
 * When we broadcast the winning model, the array of pointers in that
 * model gets sent to everyone. Since those pointers are nonsense to
 * other nodes, those nodes have to store their own pointers and then
 * restore them after the broadcast.
 *
 * The maxima array is also saved and restored here.
 *
 **********************************************************************/

void save_theta_ptrs(MODEL *model, int save_or_restore)
{
  static THETA saved_theta;
  static THETA saved_logtheta;
  static THETA saved_logtheta_rc;
  static THETA saved_obs;
  static P_PROB saved_maxima;

  if (save_or_restore == 1) { 		/* Save the theta matrix pointer */
    saved_theta = model->theta;
    saved_logtheta = model->logtheta;
    saved_logtheta_rc = model->logtheta_rc;
    saved_obs = model->obs;
    saved_maxima = model->maxima;
  } else {
    model->theta = saved_theta;
    model->logtheta = saved_logtheta;
    model->logtheta_rc = saved_logtheta_rc;
    model->obs = saved_obs;
    model->maxima = saved_maxima;
    Resize(model->maxima, model->nsites_dis, p_prob);
  }
} /* save_theta_ptrs */

/**********************************************************************
 *
 * void theta_packer
 *
 * This function packs (or unpacks) the 2-D theta matrix into a 1-D
 * array so that it can be sent in a single broadcast.
 * 
 **********************************************************************/
void theta_packer(
   MODEL *model,
   double *theta_matrix,
   double *obs_matrix,
   int pack_or_unpack,
   int alength
)
{
  THETA theta, obs;

  int width, i, j;

  /* An index for stepping through the linear theta matrix. */
  int i_theta = 0;

  /* Get the theta matrix for this component. */
  theta = model->theta;
  obs = model->obs;

  /* Find the width of motif. */
  width = model->w;
  
  /*
  fprintf(stderr, "%d: alength=%d width=%d\n", mpMyID(), alength, width);
  */

  /* Iterate across positions. */
  for (i = 0; i < width; i++) {
    /* Iterate through the letters in the dataset. */
    for (j = 0; j < alength; j++) {

      /* Copy from theta into the theta matrix or vice versa. */
      if (pack_or_unpack == 1) {
	theta_matrix[i_theta] = theta(i, j);
	obs_matrix[i_theta++] = obs[i][j];
      } else {
	theta(i, j) = theta_matrix[i_theta];
	obs[i][j] = obs_matrix[i_theta++];
      }
    }
  }

} /* void theta_packer */

/**********************************************************************
 *
 * void max_packets
 *
 * Given two packets, select the one with the best sig.
 * Ties are resolved in favor of:
 *   1) lowest start width
 *   2) lowest start nsites
 * in Classic mode.  Otherwise, we first break ties using
 *   1) model width
 *   2) model LLR
 * and then the previous two rules are applied.
 * This function is used in the reduction across models.
 *
 **********************************************************************/
void max_packets(void *f_data, void *f_result, int *f_length,
		  MPI_Datatype *datatype)
{
  double r_sig = ((REDUCE_PACKET *)f_result)->sig;
  double d_sig = ((REDUCE_PACKET *)f_data)->sig;
  double r_sw = ((REDUCE_PACKET *)f_result)->s_width;
  double d_sw = ((REDUCE_PACKET *)f_data)->s_width;
  double r_sn = ((REDUCE_PACKET *)f_result)->s_nsites;
  double d_sn = ((REDUCE_PACKET *)f_data)->s_nsites;
  double r_w = ((REDUCE_PACKET *)f_result)->width;
  double d_w = ((REDUCE_PACKET *)f_data)->width;
  double r_llr = ((REDUCE_PACKET *)f_result)->llr;
  double d_llr = ((REDUCE_PACKET *)f_data)->llr;
  BOOLEAN classic = ((REDUCE_PACKET *)f_result)->classic == 1;

  /* Compare the data in the two packets. */
  // These rules must agree with discretize() and
  // the logic of save_candidate() or serial and
  // parallel results can differ.
  BOOLEAN data_wins;
  if (
    (d_sig < r_sig)
    || (!classic 			// not Classic
      && (				// and data better
	(d_sig == r_sig && d_w < r_w) 
	|| (d_sig == r_sig && d_w == r_w && d_llr > r_llr)
	|| (d_sig == r_sig && d_w == r_w && d_llr == r_llr && d_sw < r_sw)
	|| (d_sig == r_sig && d_w == r_w && d_llr == r_llr && d_sw == r_sw && d_sn < r_sn)
      )
    )
    || (classic				// Classic
      && (				// and data better
	(d_sig == r_sig && d_sw < r_sw) 
	|| (d_sig == r_sig && d_sw == r_sw && d_sn < r_sn )
      )
    )
  ) {
    data_wins = TRUE;
  } else {
    data_wins = FALSE;
  }

  // Now do the reduction if data wins.
  if (data_wins) {
    /* Put the data into the result position. */
    ((REDUCE_PACKET *)f_result)->sig = d_sig;
    ((REDUCE_PACKET *)f_result)->s_width = d_sw;
    ((REDUCE_PACKET *)f_result)->s_nsites = d_sn;
    ((REDUCE_PACKET *)f_result)->width = d_w;
    ((REDUCE_PACKET *)f_result)->llr = d_llr;
    ((REDUCE_PACKET *)f_result)->classic = ((REDUCE_PACKET *)f_data)->classic;
    ((REDUCE_PACKET *)f_result)->ID = ((REDUCE_PACKET *)f_data)->ID;
  }
}

/**********************************************************************
 *
 * void reduce_across_models
 *
 * Do a reduction across an entire model. The model with the lowest
 * score in the 'logev' field gets propagated to every node.
 * Ties are resolved in favor of:
 *   1) lowest start width
 *   2) lowest start nsites
 *
 * This function is called by 'meme.c'.
 **********************************************************************/
void reduce_across_models(
   MODEL *model,
   int alength
)
{
  static int init;
  static MPI_Datatype reduction_packet_type;
  static MPI_Op max_packets_op;
  REDUCE_PACKET a_packet, best_packet;
  double theta_matrix[MAXALPH * (MAXSITE + 1)];
  double obs_matrix[MAXALPH * (MAXSITE + 1)];

  /* Add up the number of EM iterations from each node. */
  mpReduceAdd(&(model->iter));

  /* Package the best model's sig together with the processor ID. */
  // With the Classic objective function, ties are broken in
  // by selecting the model that would be evaluated *first* in serial mode.
  // For other objective functions, model width and LLR are used
  // before that to break ties.
  a_packet.sig = model->logev;
  a_packet.s_width = model->pw;
  a_packet.s_nsites = model->psites;
  a_packet.width = (double)model->w;
  a_packet.llr = model->llr;
  a_packet.classic = model->objfun==Classic ? 1.0 : 0.0;
  a_packet.ID = (double)mpMyID();
    
  /*fprintf(stdout, "%d: Prior score=%g\n", mpMyID(), a_packet.sig);*/

  /* Reduce across all the sig-ID pairs. */
  /* Make sure we have a handle for the reduction function. */
  if (init == 0) {
    init = 1;
    MPI_Type_contiguous(8, MPI_DOUBLE, &reduction_packet_type);
    MPI_Type_commit(&reduction_packet_type);
    MPI_Op_create(max_packets, TRUE, &max_packets_op);
  }

  /* Do the reduction. */
  MPI_Allreduce((void *)&a_packet, (void *)&best_packet, 1,
		reduction_packet_type, max_packets_op, MPI_COMM_WORLD);

  /*fprintf(stdout, "%d: Winner=%d score=%g\n", mpMyID(), (int)best_packet.ID,
	  best_packet.sig);*/

  /* The losers store their theta matrix and maxima pointers so they don't get
     overwritten in the upcoming broadcast. */
  if (mpMyID() != (int)best_packet.ID) save_theta_ptrs(model, 1);

  /* Broadcast the best model from the winning ID. */
  mpBroadcast((void *)model, sizeof(MODEL), (int)best_packet.ID);

  /*fprintf(stderr, "%d: Past model broadcast.\n", mpMyID());
  fflush(stderr);*/

  /* After the broadcast, the losers restore their original pointers. */
  if (mpMyID() != (int)best_packet.ID) save_theta_ptrs(model, 0);

  /* The winner packages the theta matrix up into a linear array. */
  if (mpMyID() == (int)best_packet.ID)
    theta_packer(model, theta_matrix, obs_matrix, 1, alength);

  /* Broadcast the theta and obs matrix arrays and discrete maxima. */
  mpBroadcast((void *)theta_matrix, sizeof(double) * alength * model->w,
	      (int)best_packet.ID);
  mpBroadcast((void *)obs_matrix, sizeof(double) * alength * model->w,
	      (int)best_packet.ID);
  mpBroadcast((void *)model->maxima, sizeof(p_prob) * model->nsites_dis, 
	      (int)best_packet.ID);

  /*fprintf(stderr, "%d: Past theta broadcast.\n", mpMyID());
  fflush(stderr);*/

  /* Everyone else unpacks the theta matrix. */
  if (mpMyID() != (int)best_packet.ID)
    theta_packer(model, theta_matrix, obs_matrix, 0, alength);
}

/**********************************************************************
 *
 * void max_s_packets
 *
 * Find the maximum from two s_point packets.
 * This function is used in the reduction across starting points.
 *
 * Ties are broken by taking first starting point in dataset.
 *
 **********************************************************************/
void max_s_packets(void *f_data, void *f_result, int *f_length,
		   MPI_Datatype *datatype)
{
  int i_nsites0;
  double rs, ds, ri, di, rj, dj;

  /* Compare each s_point in the array. */
  for (i_nsites0 = 0; i_nsites0 < *f_length; i_nsites0++) {

    /* Get the two scores. */
    rs = ((S_POINT_PACKET *)f_result + i_nsites0)->score;
    ds = ((S_POINT_PACKET *)f_data + i_nsites0)->score;
    ri = ((S_POINT_PACKET *)f_result + i_nsites0)->iseq;
    di = ((S_POINT_PACKET *)f_data + i_nsites0)->iseq;
    rj = ((S_POINT_PACKET *)f_result + i_nsites0)->ioff;
    dj = ((S_POINT_PACKET *)f_data + i_nsites0)->ioff;

    /* Compare the two scores breaking ties so that starting points
       tested first in serial mode win. */
    int data_wins = FALSE;
    if (ds > rs || (ds == rs && di < ri) || (ds == rs && di == ri && dj < rj)) {
      /* Copy the new data into the result location. */
      ((S_POINT_PACKET *)f_result + i_nsites0)->score = ds;
      ((S_POINT_PACKET *)f_result + i_nsites0)->iseq = di;
      ((S_POINT_PACKET *)f_result + i_nsites0)->ioff = dj;
      data_wins = TRUE;
    }

  }
} // max_s_packets

/**********************************************************************
 *
 * void reduce_across_s_points
 *
 * Do a reduction across an array of starting points. For each
 * position in the array, the starting point with the largest value in
 * the 'score' field gets propagated to every node.
 *
 * This function is called by subseq7.c.
 **********************************************************************/
void reduce_across_s_points(
  S_POINT *s_points, 
  SAMPLE **samples, 
  int n_nsites0,
  int n_starts
)
{
  static int init;
  static MPI_Datatype s_point_packet_type;
  static MPI_Op max_s_packets_op;
  int i_packet;
  S_POINT_PACKET packets[1000], best_packets[1000];
  
  /* 
     NOTE: sizeof(S_POINT) = 344 bytes
     sizeof(S_POINT_PACKET) = 16 bytes
  */
  
  /* Package the scores into an array of packets. */
  /* Don't need w0 & nsites0 because they're invariant across processors.*/
  /* Don't want e_cons0 & cons0 because they are pointers
     (We recalculate e_cons0 and contents of cons0 later). */
  for (i_packet = 0; i_packet < n_nsites0; i_packet++) {
    packets[i_packet].score = s_points[i_packet].score;
    packets[i_packet].iseq = (double)s_points[i_packet].iseq;
    packets[i_packet].ioff = (double)s_points[i_packet].ioff;
  }

  /*
  printf("BEFORE\n");
  for (i_packet = 0; i_packet < n_nsites0; i_packet++)
    fprintf(stdout, "node %d, packet %d: score=%g iseq=%d ioff=%d\n",
	    mpMyID(), i_packet, packets[i_packet].score,
	    (int)packets[i_packet].iseq, (int)packets[i_packet].ioff);
  fflush(stdout);
  */

  /* MPI */
  /* Reduce across the list of packets. */
  if (init == 0) {
    init = 1;
    MPI_Type_contiguous(3, MPI_DOUBLE, &s_point_packet_type);
    MPI_Type_commit(&s_point_packet_type);
    MPI_Op_create(max_s_packets, TRUE, &max_s_packets_op);
  }

  /* Do the reduction. */
  MPI_Allreduce((void *)&packets, (void *)&best_packets, n_nsites0,
		s_point_packet_type, max_s_packets_op, MPI_COMM_WORLD);

  /*
  printf("AFTER\n");
  for (i_packet = 0; i_packet < n_nsites0; i_packet++)
    fprintf(stdout, "node %d, packet %d: score=%g iseq=%d ioff=%d\n",
	    mpMyID(), i_packet, best_packets[i_packet].score,
	    (int)best_packets[i_packet].iseq,
	    (int)best_packets[i_packet].ioff);
  fflush(stdout);
  */
  
  /* Unpack the data. */
  for (i_packet = 0; i_packet < n_nsites0; i_packet++) {
    s_points[i_packet].score = best_packets[i_packet].score;
    s_points[i_packet].iseq = (int)best_packets[i_packet].iseq;
    s_points[i_packet].ioff = (int)best_packets[i_packet].ioff;

    /* Set the e_cons0 fields to the proper memory locations. */
    s_points[i_packet].e_cons0 =
      samples[s_points[i_packet].iseq]->res + s_points[i_packet].ioff;
  }
  
  /* Add up all of the values of n_starts. */
  mpReduceAdd(&n_starts);
}

/**********************************************************************/
/*
	get_start_n_end

	Get the precalculated astarting and ending sequences and offsets.
*/
/**********************************************************************/
void get_start_n_end (
  int *start_seq,			/* starting sequence number */
  int *start_off,			/* offset in starting sequence */
  int *end_seq,				/* ending sequence number */
  int *end_off 				/* offset in ending sequence */
)
{
  *start_seq = start_n_end->start_seq;
  *start_off = start_n_end->start_off;
  *end_seq = start_n_end->end_seq;
  *end_off = start_n_end->end_off;
} /* get_start_n_end */

/**********************************************************************
 * void print_model
 *
 * Print some of the contents of a model.
 * (Debugging only.)
 **********************************************************************/
void print_model(
   char *label,
   MODEL *model
)
{		       
  fprintf(stderr, "%d: %s -- sig=%1.2e", mpMyID(), label, model->logev);
  fflush(stdout);
}

#endif /* PARALLEL */
