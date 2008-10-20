/*****************************************************************
 **
 ** update ystar
 **
 ** simon jackman, dept of political science, stanford university
 ** mar 2001
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Rmath.h>
#include "util.h"
#include "ideal.h"

/* putting stuff in y star */

double updatey(double **ystar, double **y, double **x, double **beta,
	       double **xHat, double **bHat, double **z,
	       double sd,
	       int n, int m, int d, int iter)
{
  int i,j,k;
  double *xrow, *brow,
    *xHatRow, *bHatRow,
    mu, s, e, muhat;
  //  float z;

  s = 0.0;
  e = 0.0;
  for(i=0;i<n;i++){               /* loop over legislators */
    xrow = x[i];
    xHatRow = xHat[i];

    for(j=0;j<m;j++){             /* loop over proposals */
      brow = beta[j];
      bHatRow = bHat[j];
      
      mu = -1.0*sd*brow[d];          /* intercept */ /*d+1*/
      muhat = -1.0*bHatRow[d];
      for(k=0;k<d;k++){
	mu += brow[k]*sd*xrow[k];
	muhat += bHatRow[k]*xHatRow[k];
      }
      if (y[i][j]==9.0){         /* sample untruncated, missing responses */
	ystar[i][j] = rnorm(mu,1.0);
      }
      else{                      /* sample from truncated normals */
	ystar[i][j] = dtnorm(&mu,&sd,&y[i][j]);  /* try two */
      }
      z[i][j] = ystar[i][j]/sd;
      e = z[i][j] - muhat;
      s += (e*e);
    }
  }
  return(s);
}
