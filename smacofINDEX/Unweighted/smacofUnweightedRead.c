#include "smacofUnweighted.h"

// read dissimilarities in standard order from a file

void smacofUnweightedReadInputFile(char *fname, double *delta) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    unsigned k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    fclose(stream);
    return;
}
