#include "smacofWeighted.h"

void smacofWeightedReadInputFile(char *fname, unsigned *irow, unsigned *icol,
                                 double *delta, double *weights) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    unsigned k = 0;
    fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k], &weights[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k],
               &weights[k]);
    }
    fclose(stream);
    return;
}

