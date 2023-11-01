#include <stdio.h>
#include <stdlib.h>

void smacofReadInputFile(char *fname, int *irow, int *icol, double *delta, double *weights) {
	FILE *stream = fopen(fname, "r");
	if (stream == NULL) {
		printf("Error: cannot open %s\n", fname);
		exit(EXIT_FAILURE);
	}
	int k = 0;
	fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k], &weights[k]);
	while(!feof(stream)) {
		k++;
	    fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k], &weights[k]);
	}
	fclose(stream);
	return;
}

int main() {
	double delta[3] = {0}; 
	double weights[3] = {0}; 
	int irow[3] = {0}; 
	int icol[3] = {0};
	char fname[] = "ekman.txt";
	(void)smacofReadInputFile(fname, irow, icol, delta, weights);
	for (int i = 0; i < 3; i++) {
		printf("%3d %3d %4.0f %4.0f\n", irow[i], icol[i], delta[i], weights[i]);
	}
	return(EXIT_SUCCESS);
}
