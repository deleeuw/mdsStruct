#include "smacofCommon.h"

void smacofPrintAnyMatrix(const unsigned n, const unsigned p,
                          const unsigned width, const unsigned precision,
                          const double (**x) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned s = 0; s < p; s++) {
            printf("%+*.*f", width, precision, x[i][s]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

// read data from a file into a one-dimensional array

void smacofUnweightedReadInputFile(const char *fname, double *delta) {
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

void smacofReadParameterFile(char *name, unsigned *n, unsigned *p, unsigned *itmax, unsigned *init, unsigned *feps,
                          unsigned *ceps, unsigned *width, unsigned *precision, unsigned *verbose,
                          unsigned *relax) {
    FILE *fp;
    char *line = (char *)malloc(80 * sizeof(char));
    assert(line != NULL);
    char aux1[20], aux2[10];
    unsigned par = 0;
    strcat(name, ".par");
    size_t len = 0;
    ssize_t read;
    fp = fopen(name, "r");
    assert(fp != NULL);
    while ((read = getline(&line, &len, fp)) != -1) {
        sscanf(line, "%s %s %d", aux1, aux2, &par);
        if (strcmp("n", aux1) == 0) {
            *n = par;
        }
        if (strcmp("p", aux1) == 0) {
            *p = par;
        }
        if (strcmp("itmax", aux1) == 0) {
            *itmax = par;
        }
        if (strcmp("init", aux1) == 0) {
            *init = par;
        }
        if (strcmp("feps", aux1) == 0) {
            *feps = par;
        }
        if (strcmp("ceps", aux1) == 0) {
            *ceps = par;
        }
        if (strcmp("width", aux1) == 0) {
            *width = par;
        }
        if (strcmp("precision", aux1) == 0) {
            *precision = par;
        }
        if (strcmp("verbose", aux1) == 0) {
            *verbose = par;
        }
        if (strcmp("relax", aux1) == 0) {
            *relax = par;
        }
    }
    fclose(fp);
    free(line);
    return;
}

