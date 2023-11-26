#include "smacofAlt.h"

void smacofJacobi(double *a, double *evec, double *eval, const int *pn,
                  const int *pm, const int *pitmax, const int *peps,
                  const bool *pverbose) {
    int n = *pn, m = *pm, itel = 1, itmax = *pitmax;
    bool verbose = *pverbose;
    double d = 0.0, s = 0.0, t = 0.0, u = 0.0, v = 0.0, p = 0.0, q = 0.0,
           r = 0.0;
    double fold = 0.0, fnew = 0.0, eps = pow(10.0, -(double)*peps);
    double *oldi = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double *oldj = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            evec[MINDEX(i, j, n)] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 1; i <= m; i++) {
        fold += SQUARE(a[TINDEX(i, i, n)]);
    }
    while (true) {
        for (int j = 1; j <= m; j++) {
            for (int i = j + 1; i <= n; i++) {
                p = a[TINDEX(i, j, n)];
                q = a[TINDEX(i, i, n)];
                r = a[TINDEX(j, j, n)];
                if (fabs(p) < 1e-10) continue;
                d = (q - r) / 2.0;
                s = (p < 0) ? -1.0 : 1.0;
                t = -d / sqrt(SQUARE(d) + SQUARE(p));
                u = sqrt((1 + t) / 2);
                v = s * sqrt((1 - t) / 2);
                for (int k = 1; k <= n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    oldi[VINDEX(k)] = a[TINDEX(ki, ik, n)];
                    oldj[VINDEX(k)] = a[TINDEX(kj, jk, n)];
                }
                for (int k = 1; k <= n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    a[TINDEX(ki, ik, n)] =
                        u * oldi[VINDEX(k)] - v * oldj[VINDEX(k)];
                    a[TINDEX(kj, jk, n)] =
                        v * oldi[VINDEX(k)] + u * oldj[VINDEX(k)];
                }
                for (int k = 1; k <= n; k++) {
                    oldi[VINDEX(k)] = evec[MINDEX(k, i, n)];
                    oldj[VINDEX(k)] = evec[MINDEX(k, j, n)];
                    evec[MINDEX(k, i, n)] =
                        u * oldi[VINDEX(k)] - v * oldj[VINDEX(k)];
                    evec[MINDEX(k, j, n)] =
                        v * oldi[VINDEX(k)] + u * oldj[VINDEX(k)];
                }
                a[TINDEX(i, i, n)] =
                    SQUARE(u) * q + SQUARE(v) * r - 2 * u * v * p;
                a[TINDEX(j, j, n)] =
                    SQUARE(v) * q + SQUARE(u) * r + 2 * u * v * p;
                a[TINDEX(i, j, n)] =
                    u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
            }
        }
        fnew = 0.0;
        for (int i = 1; i <= m; i++) {
            fnew += SQUARE(a[TINDEX(i, i, n)]);
        }
        if (verbose) {
            printf("itel %3d fold %15.10f fnew %15.10f\n", itel, fold, fnew);
        }
        if (((fnew - fold) < eps) || (itel == itmax)) {
            break;
        }
        fold = fnew;
        itel++;
    }
    for (int i = 1; i <= n; i++) {
        eval[VINDEX(i)] = a[TINDEX(i, i, n)];
    }
    free(oldi);
    free(oldj);
    return;
}

void smacofGramSchmidt(double *x, double *r, int *pn, int *pp) {
    int n = *pn, p = *pp, s = 0;
    while (s < p) {
        for (int t = 0; t < s; t++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += x[i][t] * x[i][t];
            }
            for (int i = 0; i < n; i++) {
                x[i][t] -= sum * x[i][t];
            }
        }
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += SQUARE(x[i][s]); 
        }
        sum = sqrt(sum);
        r[s] = sum;
        for (int i = 1; i <= n; i++) {
            x[i][s] /= sum;
        }
        s++;
    }
    return;
}

void smacofCenter(double *x, const int *np, const int *pp) {
    int n = *np, p = *pp;
    for (int s = 0; s < p; s++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i][s];
        }
        sum /= (double)n;
        for (int i = 0; i < n; i++) {
            x[i][s] -= sum;
        }
    }
    return;
}

void smacofMultiplySymmetricMatrixAlt(const double *a, const double *x, double *y,
                                   const int *pn, const int *pp) {
    int n = *pn, p = *pp;
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += a[i][j] * x[j][s];
            }
            y[i][s] = sum;
        }
    }
    return;
}

void smacofDistanceAlt(const double *x, double *d, const int *pn, const int *pp) {
    int n = *pn, p = *pp;
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double sum = 0.0;
            for (int s = 0; s < p; s++) {
                sum += SQUARE(x[i][s] - x[j][s]);
            }
            double dij = sqrt(fabs(sum))
            d[i][j] = dij;
            d[j][i] = dij;
        }
    }
    return;
}
