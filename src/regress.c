#include "macros.h"

/*
	regress

	Least squares regression on points (x,y) to give
		y = mx + b

	Returns the root mean squared error of the fit.
*/
double regress(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *m,			/* slope */
  double *b 			/* y intercept */
)
{
  int i;
  double sx=0, sy=0, sxx=0, sxy=0;
  double mse=0;

  for (i=0; i<n; i++) {
    sx += x[i];
    sy += y[i];
    sxx += x[i]*x[i];
    sxy += x[i]*y[i];
  }

  *m = (n*sxy - sy*sx) / (n*sxx - sx*sx);
  *b = (sy - *m*sx)/n;

  for (i=0; i<n; i++) {
    double err = y[i] - (*m*x[i] + *b);
    mse += err * err;
  }
  mse = sqrt(mse);
  mse /= n;

  return mse;
}
/*
	w_regress

	Weighted least squares regression on points (x,y,w) to give
		y = mx + b

	Returns the weighted root mean squared error of the fit.
*/
double w_regress(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y,			/* y values */
  double *w,			/* weights */
  double *m,			/* slope */
  double *b 			/* y intercept */
)
{
  int i;
  double s=0, sx=0, sy=0, sxx=0, sxy=0;
  double mse=0;

  for (i=0; i<n; i++) {
    double ww = w[i];
    s += ww;
    sx += x[i]*ww;
    sy += y[i]*ww;
    sxx += x[i]*x[i]*ww;
    sxy += x[i]*y[i]*ww;
  }

  *m = (s*sxy - sx*sy) / (s*sxx - sx*sx);
  *b = (sxx*sy - sx*sxy)/(s*sxx - sx*sx);

  for (i=0; i<n; i++) {
    double err = (y[i] - (*m*x[i] + *b)) * w[i];
    mse += err * err;
  }
  mse = sqrt(mse);
  mse /= n;

  return mse;
}
