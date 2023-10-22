#include "smacof.h"

void smacofPrintSDCMatrix(const double *v, const int *pn, const int *pw,
                          const int *pr) {
    int n = *pn;
    double *d = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        d[VINDEX(i)] = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            d[VINDEX(i)] += -v[PINDEX(i, j, n)];
        }
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                printf(" %+6.4f", d[VINDEX(j)]);
            } else {
                printf(" %+6.4f", v[PINDEX(i, j, n)]);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    free(d);
    return;
}

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

void smacofPrintSymmetricHollowMatrix(const double *d, const int *pn,
                                      const int *pw, const int *pr) {
    int n = *pn;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                printf(" %+6.4f", 0.0);
            } else {
                printf(" %+6.4f", d[PINDEX(i, j, n)]);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}