/*
 * $Id: dpalign.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:37:27  nadya
 * Initial revision
 *
 */

#ifndef DPALIGN_H
#  define DPALIGN_H

char *dp_align(
  ALPH_T *alph, // alphabet
  int n, // length of logodds matrix 
  double **logodds, // logodds matrix 
  char *seq1, // consensus of logodds 
  int m, // length of sequence 
  uint8_t *eseq, // integer encoded sequence 
  double wg, // gap cost (initialization)  
  double ws, // space cost (extension)  
  BOOLEAN endgaps // penalize end gaps if TRUE 
);

char **dp_multi_align(
  P_PROB sites, // the sites 
  int nsites, // the number of sites 
  int w, // width of sites 
  int flank, // add flank cols on lft+rgt 
  DATASET *dataset // the dataset 
);

int g_align(
  char **ma, // the multiple alignment 
  int nseq, // number of sequences in alignment 
  int ncol, // number of columns in alignment 
  int left, // leftmost allowed end 
  int right, // rightmost allowed end 
  int minw, // minimum width 
  int *off, // offset of g-alignment rel. to left 
  int *w // width of g-alignment 
);

#endif
