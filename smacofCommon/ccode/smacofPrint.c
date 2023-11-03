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

void smacofPrintLTMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr) {
    int n = *pn;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i >= j) {
                printf(" %+6.4f", d[TINDEX(i, j, n)]);
            } else {
                printf(" %+6.4f", 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofPrintSLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr) {
    int n = *pn;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                printf(" %+6.4f", d[SINDEX(i, j, n)]);
            } else {
                printf(" %+6.4f", 0.0);
            }
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}
