#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "ideal.h"

double**dvecTOdmat(double *vtr, double **dmtrx, int rows, int columns);
double *dmatTOdvec(double *vtr, double **dmtrx, int rows, int columns);

double *xxprod, **xxchol;
double *xz;
double *bxprod, **bchol;
double *bz;
double *bbp, **bba;
double *xxp, **xxa;
double **bpb, *xprior, **xpriormat, *xbar, **xvpost, *bpw, **w;
double **xpx, **bvpost, **bpriormat, *bprior, *bbar, *xpy;

void IDEAL(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	   int *impute1, int *meanzero1, double *xpriormeans1, 
	   double  *xpriorprec1, double *bpriormeans1, double *bpriorprec1, 
	   double *xstart1, double *bstart1, double *xoutput, double *boutput,
	   int *burnin1, int *usefile, int *bsave, char **filename1)
{
  int e, xocursor, bocursor, xlength, blength, q, nm, iter;
  int inloop, **ok, burnin,prn,n,m,d,maxiter,thin,impute,meanzero;
  double **ystar, **x, **xreg, **y, **beta, **bp, **bpv;
  double **xp, **xpv, *xtemp, *btemp;
  FILE *ofp;
   // extern double **bpb, *xprior, **xpriormat, *xbar, **xvpost, *bpw, **w;
   // extern double **xpx, **bvpost, **bpriormat, *bprior, *bbar, *xpy;


  n=*n1;
  m=*m1;
  d=*d1;
  maxiter=*maxiter1;
  thin=*thin1;
  impute=*impute1;
  meanzero=*meanzero1;
  burnin=*burnin1;

  prn=0;

  /*Creating the matrices we'll need*/
  iter = 0;                  /* initialize iter count */ 
  nm = n * m;                
  q = d + 1;                 /* item parameters, per item */ 
  y = dmatrix(n,m);      /* roll call data */
  ystar = dmatrix(n,m);  /* latent utility differential */ 
  x = dmatrix(n,d);      /* latent indicators */ 
  xreg = dmatrix(n,q);   /* regressors for updates of beta */
  beta = dmatrix(m,q);   /* item parameters */ 
  bp  = dmatrix(m,q);    /* initialize prior means, item parameters */
  bpv = dmatrix(m,q);    /* initialize prior variances, item parameters */
  xp  = dmatrix(n,d);    /* initialize prior means, latent traits */
  xpv = dmatrix(n,d);    /* initialize prior variances, latent traits */
  ok = imatrix(n,m);     /* initialize ok indicator matrix */  

  if (*usefile == 1) {
    ofp = fopen(R_ExpandFileName(*filename1), "a");
    
    if (ofp == NULL) {
      calcerror("Can't open outfile file!\n");
    }
  }

  GetRNGstate();



  /*for error checking: the parameters*/
  /*printf("Checking parameters\n");
  printf("n: %d, m: %d, d: %d, maxiter: %d, thin: %d, impute: %d, meanzero: %d\n", 
  		  n, m, d, maxiter, thin, impute, meanzero);
		  
  printf("\ny vector\n");
  for (a=0; a < nm; a++) {
  		printf("y1[%d] %g\n", a, y1[a]);
  		} 
  printf("\nbpriormeans1 vector\n");
  for(a=0; a <m*q; a++) {
  		printf("bpriormean1[%d] %g\n",a, bpriormeans1[a]);
  		}
  printf("\nbpriorprec1 vector\n");
  for(a=0; a <m*q; a++) {
  		printf("bpriorprec1[%d] %g\n",a, bpriorprec1[a]);
  		}
  printf("\nxpriormeans1 vector\n");
  for(a=0; a <n*d; a++) {
  		printf("xpriormean1[%d] %g\n",a, xpriormeans1[a]);
  		}	
  printf("\nxpriorprec1 vector\n");
  for(a=0; a <n*d; a++) {
  		printf("xpriorprec1[%d] %g\n",a, xpriorprec1[a]);
  		}
  */
  
  /*populate dmatrices with passed in parameters*/
  dvecTOdmat(y1, y, n, m);
  dvecTOdmat(bpriormeans1, bp, m, q);
  dvecTOdmat(bpriorprec1, bpv, m, q);
  dvecTOdmat(xpriormeans1, xp, n, d);
  dvecTOdmat(xpriorprec1, xpv,n ,d);
  dvecTOdmat(xstart1, x, n, d);
  dvecTOdmat(bstart1, beta, m, q);

  //for error checking: printing all the matrices
  /*
  Rprintf("\nChecking matrices\n");
  Rprintf("y matrix\n");
  printmat(y, n, m);
  Rprintf("\nbp matrix\n");
  printmat(bp, m, q);
  Rprintf("\nbpv\n");
  printmat(bpv, m, q);
  
  Rprintf("\nxp\n"); 
  printmat(xp, n, d); 
  Rprintf("\nxpv\n"); 
  printmat(xpv, n ,d); 
   */
   
  /* W R I T E   I N I T I A L   V A L U E S   T O   O U T P U T   V E C T O R S  */
  xtemp = dvector(n*d);
  xlength = n * d;
  if (burnin == 0){
    if (*usefile == 1) {
    }
    else {
      xocursor = n * d -1; /*gives location in xoutput vector*/
      dmatTOdvec(xoutput, x, n, d);
    }
  }
  else {
    xocursor = -1;
  }
  
  btemp = dvector(m*q);
  blength = m * q;
  if (burnin == 0) {
    if (*bsave == 1) {
      if (*usefile == 1) {
      }
      else {
	bocursor = m * q -1;  /*gives location in boutput vector*/
	dmatTOdvec(boutput, beta, m, q);
      }
    }
  }
  else {
    bocursor = -1;
  }
  
  check(y,ok,n,m);


  /*******************************
   * INITIALIZE REUSED VARIABLES *
   *******************************/

  bpb = dmatrix(d,d);
  bpw  = dvector(d);
  xbar = dvector(d);
  xvpost = dmatrix(d,d);
  xprior = dvector(d);
  xpriormat = dmatrix(d,d);
  w = dmatrix(n,m);

  xpy = dvector(q);
  xpx = dmatrix(q,q);
  bbar = dvector(q);
  bprior = dvector(q);
  bvpost = dmatrix(q,q);
  bpriormat = dmatrix(q,q);

  xz = dvector(d);
  xxprod = dvector(d);
  xxchol = dmatrix(d,d);

  bz = dvector(q);
  bxprod = dvector(q);
  bchol = dmatrix(q,q);

  xxp = dvector(d);
  xxa = dmatrix(d,d);

  bbp = dvector(q);
  bba = dmatrix(q,q);



  
  /**********************************************************************/
  /* C O M M E N C E   I T E R A T I O N S                              */
  /**********************************************************************/

                           
  while(iter<maxiter){                       /* Gibbs sampler loop */
    for(inloop=0;inloop<thin;inloop++){     /* thinning loop */
      iter++;                     /* increment iter count */
      if(iter>maxiter)                    /* are we done? */
	break;
      //Rprintf("\niter: %d\n",iter);
      updatey(ystar,y,x,beta,n,m,d,iter);   
      //Rprintf("past update y\n");
      
      makexreg(xreg,x,n,d,q);
      //Rprintf("past makexreg\n");
      
      updateb(ystar,ok,beta,xreg,bp,bpv,n,m,d,impute);
      //Rprintf("past updateb\n");
      
      updatex(ystar,ok,beta,x,xp,xpv,n,m,d,impute,meanzero);
      //Rprintf("past updatex\n"); 

      R_CheckUserInterrupt();               /* check for user interrupt */
      
    }
    
  
    /**********************************************************************/
    /*I N P U T I N G   N E W   V A L U E S  I N T O  E X P O R T  V E C S*/
    /**********************************************************************/ 
    if (++prn==3) {
      Rprintf("\nCurrent Iteration: %d",iter);
      prn=0;
    } 
    if (iter>=burnin) {
      // the x matrix into a vector for this iteration
      if (*usefile == 1) {
	dmatTOdvec(xtemp, x, n, d); //replace with function call
	fprintf(ofp, "%d", iter);
	for (e = 0; e < xlength; e++) {
	  fprintf(ofp, ",%f", xtemp[e]);
	} 
	if (*bsave != 1) {
	  fprintf(ofp,"\n");
	}
      }
      else {
	dmatTOdvec(xtemp, x, n, d); 
	for (e = 0; e < xlength; e++) {
	  xocursor++;
	  xoutput[xocursor] = xtemp[e];	 
	}
      }
      //the b matrix into vector form for this iteration	
      if (*bsave == 1) {
	if (*usefile == 1) { // this is not the most efficient way to do this
	  dmatTOdvec(btemp, beta, m, q); //replace with function call
	  for (e = 0; e < blength; e++) {
	    fprintf(ofp, ",%f", btemp[e]);
	  } 
	  fprintf(ofp,"\n");
	}
	else {
	  dmatTOdvec(btemp, beta, m, q);
	  for (e = 0; e < blength; e++) {
	    bocursor++;
	    boutput[bocursor] = btemp[e];	 
	  } 
	}
      }
    }
  }
  
  PutRNGstate();

  if (*usefile == 1) {
    fclose(ofp);
  }

  return;
  
}
