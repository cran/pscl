#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double r_sd(double s, double df)
{
  double root, r, g, sd;
  
  g = rchisq(df);
  r = s/g;
  root = sqrt(r);

  // Rprintf("r_sd: s,g, s/g = %14.4lf %14.4lf %14.4lf\n",s,g,r);

  //Rprintf("r_sd: how do I take the square root of something? %14.10lf\n",
  //    sqrt(2.0));

  // Rprintf("r_sd: how do I take the square root of something? %14.10lf\n",
  //	  sqrt(r));
  
  // root = sqrt(exp(rnorm(0.0,4.0)));

  // root = runif(0.5,5.0);

  return(root);
}
