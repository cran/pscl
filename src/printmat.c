#include <stdio.h>
#include "util.h"

void printmat(double **mat, int nr, int nc)
{
  int i,j;

  for(i=0;i<nr;i++){
    for(j=0;j<nc;j++){
      fprintf(stdout,
	      "mat[%d][%d]=%2.3lf ",
	      i,j,mat[i][j]);
    }
    fprintf(stdout,"\n");
  }
  return;
}
