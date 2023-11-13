#include "../Include/smacof.h"

void smacofPrintAnyMatrix(const size_t n, const size_t p, const int width,
                          const int precision, double (*x)[n]) {
    for (size_t i = 0; i < n; i++) {
        for (size_t s = 0; s < p; s++) {
            printf("%+*.*f", width, precision, x[i][s]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

// read data from a file into a one-dimensional array

void smacofUnweightedReadInputFile(char *fname, double *delta) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    int k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    fclose(stream);
    return;
}
