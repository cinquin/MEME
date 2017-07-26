/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995-2014, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
 
#ifndef REGRESS_H 
#define REGRESS_H

extern double regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

double w_regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *w,                    /* weights */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

#endif
 
