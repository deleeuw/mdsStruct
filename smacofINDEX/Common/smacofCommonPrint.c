#include "smacofCommon.h"

void smacofPrintAnyMatrix(const unsigned n, const unsigned p, const double *x,
                          const unsigned width, const unsigned precision) {
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned s = 1; s <= p; s++) {
            printf(" %+*.*f", width, precision, x[MINDEX(i, s, n)]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSymmetricMatrix(const unsigned n, const double *x,
                                const unsigned width,
                                const unsigned precision) {
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned j = 1; j <= n; j++) {
            unsigned ij = UINDEX(i, j, n);
            printf(" %+*.*f", width, precision, x[ij]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSHMatrix(const unsigned n, const double *x,
                         const unsigned width, const unsigned precision) {
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned j = 1; j <= n; j++) {
            if (i == j) {
                printf(" %+*.*f", width, precision, 0.0);
            } else {
                printf(" %+*.*f", width, precision, x[PINDEX(i, j, n)]);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

/*
void smacofPrintSHMatrixIJ(const double *d, const unsigned *pn, unsigned *pm,
                           const unsigned *irow, const unsigned *icol, const
unsigned *pw, const unsigned *pr) { unsigned n = *pn, m = *pm, nn = n * (n - 1)
/ 2; double *dnan = (double *)calloc((size_t)nn, (size_t)sizeof(double)); for
(unsigned k = 1; k <= nn; k++) { dnan[VINDEX(k)] = NaN;
    }
    for (unsigned k = 1; k <= m; k++) {
        unsigned i = irow[VINDEX(k)], j = icol[VINDEX(k)];
        dnan[SINDEX(i, j, n)] = d[VINDEX(k)];
    }
    (void)smacofPrintSHMatrix(dnan, pn, pw, pr);
    free(dnan);
    return;
}
*/

void smacofPrintLTMatrix(const unsigned n, const double *x,
                         const unsigned width, const unsigned precision) {
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned j = 1; j <= n; j++) {
            if (i >= j) {
                printf(" %*.*f", width, precision, x[TINDEX(i, j, n)]);
            } else {
                printf(" %*.*f", width, precision, 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSLTMatrix(const unsigned n, const double *x,
                          const unsigned width, const unsigned precision) {
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned j = 1; j <= n; j++) {
            if (i > j) {
                printf(" %*.*f", width, precision, x[SINDEX(i, j, n)]);
            } else {
                printf(" %*.*f", width, precision, 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

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
