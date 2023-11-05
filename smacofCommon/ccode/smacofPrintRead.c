#include "smacof.h"

void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp,
                          const int *pw, const int *pr) {
    int n = *pn, p = *pp, width = *pw, precision = *pr;
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            printf(" %+*.*f", width, precision, x[MINDEX(i, s, n)]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSymmetricMatrix(const double *x, const int *pn, const int *pw,
                                const int *pr) {
    int n = *pn, width = *pw, precision = *pr;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            int ij = UINDEX(i, j, n);
            printf(" %+*.*f", width, precision, x[ij]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSHMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr) {
    int n = *pn, width = *pw, precision = *pr;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                printf(" %+*.*f", width, precision, 0.0);
            } else {
                printf(" %+*.*f", width, precision, d[PINDEX(i, j, n)]);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSHMatrixIJ(const double *d, const int *pn, int *pm,
                           const int *irow, const int *icol, const int *pw,
                           const int *pr) {
    int n = *pn, m = *pm, nn = n * (n - 1) / 2;
    double *dnan = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    for (int k = 1; k <= nn; k++) {
        dnan[VINDEX(k)] = NaN;
    }
    for (int k = 1; k <= m; k++) {
        int i = irow[VINDEX(k)], j = icol[VINDEX(k)];
        dnan[SINDEX(i, j, n)] = d[VINDEX(k)];
    }
    (void)smacofPrintSHMatrix(dnan, pn, pw, pr);
    free(dnan);
    return;
}

void smacofPrintLTMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr) {
    int n = *pn, width = *pw, precision = *pr;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i >= j) {
                printf(" %*.*f", width, precision, d[TINDEX(i, j, n)]);
            } else {
                printf(" %*.*f", width, precision, 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr) {
    int n = *pn, width = *pw, precision = *pr;
    ;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                printf(" %*.*f", width, precision, d[SINDEX(i, j, n)]);
            } else {
                printf(" %*.*f", width, precision, 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofWeightedReadInputFile(char *fname, int *irow, int *icol, double *delta,
                         double *weights) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    int k = 0;
    fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k], &weights[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%d %d %lf %lf", &irow[k], &icol[k], &delta[k],
               &weights[k]);
    }
    fclose(stream);
    return;
}

void smacofUnweightedReadInputFile(char *fname, int *irow, int *icol, double *delta) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    int k = 0;
    fscanf(stream, "%d %d %lf", &irow[k], &icol[k], &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%d %d %lf", &irow[k], &icol[k], &delta[k]);
    }
    fclose(stream);
    return;
}