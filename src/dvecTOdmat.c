double**dvecTOdmat(double *vtr, double **dmtrx, int rows, int columns)
{

    int i, j, counter;
	
	
	counter = 0;
	for (i=0; i < rows; i++) {
		for(j=0; j < columns; j++) {
				dmtrx[i][j]= vtr[counter];
				counter++;
		}
	}
	return dmtrx;


}
