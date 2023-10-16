#include "smacof.h"

// Computes the Moore-Penrose inverse of an SDC
// matrix of order n and rank n-1 from only the elements
// below the diagonal

void smacofMPInverseSDCMatrix(const double *w, double *vinv, const int *pn) {
    int n = *pn, m = n * (n - 1) / 2, ik = 0, jk = 0, ij = 0;
    double add = 1.0 / ((double)n), piv = 0.0;
    double *d = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= m; i++) {
        vinv[VINDEX(i)] = -w[VINDEX(i)];
    }
    for (int i = 1; i <= n; i++) {
        d[VINDEX(i)] = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            d[VINDEX(i)] += vinv[PINDEX(i, j, n)];
        }
        d[VINDEX(i)] = -d[VINDEX(i)] + add;
    }
    for (int i = 1; i <= m; i++) {
        vinv[VINDEX(i)] += add;
    }
    for (int k = 1; k <= n; k++) {
        piv = d[VINDEX(k)];
        for (int j = 1; j <= (n - 1); j++) {
            if (j == k) {
                continue;
            }
            jk = PINDEX(j, k, n);
            for (int i = (j + 1); i <= n; i++) {
                if (i == k) {
                    continue;
                }
                ik = PINDEX(i, k, n);
                ij = SINDEX(i, j, n);
                vinv[ij] = vinv[ij] - vinv[ik] * vinv[jk] / piv;
            }
        }
        for (int i = 1; i <= n; i++) {
            if (i == k) {
                continue;
            }
            ik = PINDEX(i, k, n);
            d[VINDEX(i)] = d[VINDEX(i)] - vinv[ik] * vinv[ik] / piv;
        }
        for (int i = 1; i <= n; i++) {
            if (i == k) {
                continue;
            }
            ik = PINDEX(i, k, n);
            vinv[ik] = vinv[ik] / piv;
        }
        d[VINDEX(k)] = -1 / piv;
    }
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            ij = SINDEX(i, j, n);
            vinv[ij] = vinv[ij] + add;
        }
    }
    for (int i = 1; i <= n; i++) {
        d[VINDEX(i)] = d[VINDEX(i)] + add;
    }
    free(d);
    return;
}

void smacofPrintSDCMatrix(const double *v, const int *pn, const int *pw,
                          const int *pr) {
    int n = *pn, m = n * (n - 1) / 2;
    double *d = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        d[VINDEX(i)] = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            d[VINDEX(i)] += v[PINDEX(i, j, n)];
        }
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                printf(" %+6.4f", d[VINDEX(j)]);
            } else {
                printf(" %+6.4f", -v[PINDEX(i, j, n)]);
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
    int n = *pn, m = n * (n - 1) / 2;
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

void smacofMultiplySDCMatrix(const double *a, const double *x, double *y,
                             const int *pn, const int *pp) {
    int n = *pn, p = *pp;
    for (int s = 1; s <= p; s++) {
        for (int i = 1; i <= n; i++) {
            double sum = 0.0;
            for (int j = 1; j <= n; j++) {
                if (j == i) {
                    continue;
                }
                int ij = PINDEX(i, j, n);
                sum += a[ij] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            y[MINDEX(i, s, n)] = sum;
        }
    }
    return;
}

void smacofRootMeanSquare(const double *x, const double *y, const double *w,
                          const int *pn, const int *pp, double *change) {
    int n = *pn, p = *pp;
    double sum = 0.0, sumsum = 0.0, dis = 0.0, djs = 0.0;
    for (int s = 1; s <= p; s++) {
        for (int i = 1; i <= n; i++) {
            sum = 0.0;
            dis = x[MINDEX(i, s, n)] - y[MINDEX(i, s, n)];
            for (int j = 1; j <= n; j++) {
                if (i == j) {
                    continue;
                }
                djs = x[MINDEX(j, s, n)] - y[MINDEX(j, s, n)];
                sum += w[PINDEX(i, j, n)] * (dis - djs);
            }
            sumsum += sum * x[MINDEX(i, s, n)];
        }
    }
    *change = sqrt(sumsum / ((double)n));
    return;
}

void smacofMaxDifference(const double *x, const double *y, double *maxdiff,
                         const int *pn, const int *pp) {
    *maxdiff = 0.0;
    int n = *pn, p = *pp;
    for (int s = 1; s <= p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (int i = 1; i <= n; i++) {
            posdiff =
                MAX(posdiff, fabs(x[MINDEX(i, s, n)] - y[MINDEX(i, s, n)]));
        }
        for (int i = 1; i <= n; i++) {
            negdiff =
                MAX(negdiff, fabs(x[MINDEX(i, s, n)] + y[MINDEX(i, s, n)]));
        }
        *maxdiff = MAX(*maxdiff, MIN(posdiff, negdiff));
    }
    return;
}

/*

int main() {
    double xold[8] = {1.0, 2.0, 3.0, 4.0, 1.0, 4.0, 9.0, 16.0};
    double xnew[8] = {1.0, 2.0, 3.0, 4.0, 1.0, 4.0, 9.0, 16.0};
    double r[2] = {0.0, 0.0}, eps = 1e-6;
    double cross[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int n = 4, p = 2, itmax = 100, width = 10, precision = 6;
    (void)smacofPrintSDCMatrix(cross, &n, &width, &precision);
    (void)smacofGramSchmidt(xold, r, &n, &p);
    (void)smacofPrintAnyMatrix(xold, &n, &p, &width, &precision);
    (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
    (void)smacofGramSchmidt(xnew, r, &n, &p);
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    (void)memcpy(xold, xnew, (size_t) (n * p * sizeof(double)));
    printf("got here\n\n");
    (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
    (void)smacofGramSchmidt(xnew, r, &n, &p);
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    (void)memcpy(xold, xnew, (size_t) (n * p * sizeof(double)));
    printf("and here\n\n");
    (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
    (void)smacofGramSchmidt(xnew, r, &n, &p);
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    (void)memcpy(xold, xnew, (size_t) (n * p * sizeof(double)));
    printf("and here\n\n");
    (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
    (void)smacofGramSchmidt(xnew, r, &n, &p);
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    (void)memcpy(xold, xnew, (size_t) (n * p * sizeof(double)));
    printf("and here\n\n");
    (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
    (void)smacofGramSchmidt(xnew, r, &n, &p);
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    (void)memcpy(xold, xnew, (size_t) (n * p * sizeof(double)));
    (void)smacofSimultaneousIteration(cross, r, &n, &p, &itmax, &eps);
    return (EXIT_SUCCESS);
}

 */