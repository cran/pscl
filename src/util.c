#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Print.h>
#include <R_ext/Error.h>

void memallocerror()
{
  error("Memory allocation error.\n");
}

void calcerror(char error_text[])
{
  error(error_text);
}

int *ivector(long n)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc(n*sizeof(int));
	if (!v) memallocerror();
	return v;
}
double *dvector(long n)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc(n*sizeof(double));
	if (!v) memallocerror();
	return v;
}

double **dmatrix(long nr, long nc)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i;
	double **m;

	m=(double **) malloc(nr*sizeof(double*));
	if (!m) memallocerror();
	for(i=0; i<nr; i++) {
	  m[i]=(double *) malloc(nc*sizeof(double));
	}
	return m;
}

int **imatrix(long nr, long nc)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i;
	int **m;

	m=(int **) malloc(nr*sizeof(int*));
	if (!m) memallocerror();
	for(i=0; i<nr; i++) {
	  m[i]=(int *) malloc(nc*sizeof(int));
	}
	return m;
}

void free_dmatrix(double **m, long nr)
/* free a double matrix allocated by dmatrix() */
{
  long i;

  for(i=0; i<nr; i++) {
    free(m[i]);
  }
  free(m);
}

void free_imatrix(int **m, long nr)
/* free an int matrix allocated by imatrix() */
{
  long i;

  for(i=0; i<nr; i++) {
    free(m[i]);
  }
  free(m);
}
