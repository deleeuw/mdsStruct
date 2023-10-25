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
                sum += -a[ij] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            y[MINDEX(i, s, n)] = sum;
        }
    }
    return;
}

void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      double *maxdiff, const int *pn,
                                      const int *pp) {
    *maxdiff = 0.0;
    int n = *pn, p = *pp;
    for (int s = 1; s <= p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (int i = 1; i <= n; i++) {
            posdiff =
                MAX(posdiff, fabs(x[MINDEX(i, s, n)] - y[MINDEX(i, s, n)]));
            negdiff =
                MAX(negdiff, fabs(x[MINDEX(i, s, n)] + y[MINDEX(i, s, n)]));
        }
        *maxdiff = MAX(*maxdiff, MIN(posdiff, negdiff));
    }
    return;
}

void smacofMaxDistanceDifference(const double *dold, const double *dnew,
                                 double *pdchange, const int *pm) {
    double dchange = 0.0;
    int m = *pm;
    for (int k = 1; k <= m; k++) {
        dchange = MAX(dchange, fabs(dold[VINDEX(k)] - dnew[VINDEX(k)]));
    }
    *pdchange = dchange;
    return;
}

void smacofAddSDCDiagonal(const double *a, double *b, const int *pn) {
    int n = *pn;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            sum += a[PINDEX(i, j, n)];
        }
        rsum[VINDEX(i)] = sum;
    }
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            int ij = TINDEX(i, j, n);
            if (i == j) {
                b[ij] = -rsum[VINDEX(i)];
            } else {
                b[ij] = a[SINDEX(i, j, n)];
            }
        }
    }
    free(rsum);
    return;
}

void smacofMakeIIandJJ(const int *pn, int *ii, int *jj) {
    int n = *pn;
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            int ij = SINDEX(i, j, n);
            ii[ij] = i;
            jj[ij] = j;
        }
    }
    return;
}

/*
 int main() {
    int n = 4, p = 2, width = 10, precision = 6;
    double v[6] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0};
    double x[8] = {1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0};
    double y[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    (void)smacofPrintSDCMatrix(v, &n, &width, &precision);
    (void)smacofPrintAnyMatrix(x, &n, &p, &width, &precision);
    (void)smacofMultiplySDCMatrix(v, x, y, &n, &p);
    (void)smacofPrintAnyMatrix(y, &n, &p, &width, &precision);
    double sum = 0.0;
    for (int i = 1; i <= 4; i++) {
      for (int s = 1; s <= 2; s++) {
        sum += x[MINDEX(i, s, 4)] * y[MINDEX(i, s, 4)];
      }
    }
    printf("%10.6f\n", sum);
}
*/