#include "smacofCommon.h"

void smacofJacobi(const unsigned n, const unsigned m, double *a, double *evec,
                  double *eval, const unsigned itmax, const unsigned peps,
                  const bool verbose) {
    unsigned itel = 1;
    double d = 0.0, s = 0.0, t = 0.0, u = 0.0, v = 0.0, p = 0.0, q = 0.0,
           r = 0.0;
    double fold = 0.0, fnew = 0.0, eps = pow(10.0, -(double)peps);
    double *oldi = (double *)calloc((size_t)n, (size_t)sizeof(double));
    assert(oldi != NULL);
    double *oldj = (double *)calloc((size_t)n, (size_t)sizeof(double));
    assert(oldj != NULL);
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned j = 1; j <= n; j++) {
            evec[MINDEX(i, j, n)] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (unsigned i = 1; i <= m; i++) {
        fold += SQUARE(a[TINDEX(i, i, n)]);
    }
    while (true) {
        for (unsigned j = 1; j <= m; j++) {
            for (unsigned i = j + 1; i <= n; i++) {
                p = a[TINDEX(i, j, n)];
                q = a[TINDEX(i, i, n)];
                r = a[TINDEX(j, j, n)];
                if (fabs(p) < 1e-10) continue;
                d = (q - r) / 2.0;
                s = (p < 0) ? -1.0 : 1.0;
                t = -d / sqrt(SQUARE(d) + SQUARE(p));
                u = sqrt((1 + t) / 2);
                v = s * sqrt((1 - t) / 2);
                for (unsigned k = 1; k <= n; k++) {
                    unsigned ik = MIN(i, k);
                    unsigned ki = MAX(i, k);
                    unsigned jk = MIN(j, k);
                    unsigned kj = MAX(j, k);
                    oldi[VINDEX(k)] = a[TINDEX(ki, ik, n)];
                    oldj[VINDEX(k)] = a[TINDEX(kj, jk, n)];
                }
                for (unsigned k = 1; k <= n; k++) {
                    unsigned ik = MIN(i, k);
                    unsigned ki = MAX(i, k);
                    unsigned jk = MIN(j, k);
                    unsigned kj = MAX(j, k);
                    a[TINDEX(ki, ik, n)] =
                        u * oldi[VINDEX(k)] - v * oldj[VINDEX(k)];
                    a[TINDEX(kj, jk, n)] =
                        v * oldi[VINDEX(k)] + u * oldj[VINDEX(k)];
                }
                for (unsigned k = 1; k <= n; k++) {
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
        for (unsigned i = 1; i <= m; i++) {
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
    for (unsigned i = 1; i <= n; i++) {
        eval[VINDEX(i)] = a[TINDEX(i, i, n)];
    }
    free(oldi);
    free(oldj);
    return;
}

void smacofGramSchmidt(const unsigned n, const unsigned p, double *x,
                       double *q) {
    for (unsigned s = 1; s <= p; s++) {
        if (s > 0) {
            for (unsigned t = 1; t <= s; t++) {
                double sum = 0.0;
                for (unsigned i = 1; i <= n; i++) {
                    sum += x[MINDEX(i, s, n)] * x[MINDEX(i, t, n)];
                }
                q[MINDEX(t, s, p)] = sum;
                for (unsigned i = 1; i <= n; i++) {
                    x[MINDEX(i, s, n)] -= sum * x[MINDEX(i, t, n)];
                }
            }
        }
        double sum = 0.0;
        for (unsigned i = 1; i <= n; i++) {
            sum += SQUARE(x[MINDEX(i, s, n)]);
        }
        sum = sqrt(sum);
        q[MINDEX(s, s, p)] = sum;
        for (unsigned i = 1; i <= n; i++) {
            x[MINDEX(i, s, n)] /= sum;
        }
    }
}

void smacofCenter(const unsigned n, const unsigned p, double *x) {
    for (unsigned s = 1; s <= p; s++) {
        double sum = 0.0;
        for (unsigned i = 1; i <= n; i++) {
            sum += x[MINDEX(i, s, n)];
        }
        sum /= (double)n;
        for (unsigned i = 1; i <= n; i++) {
            x[MINDEX(i, s, n)] -= sum;
        }
    }
    return;
}

void smacofMultiplySymmetricMatrix(const unsigned n, const unsigned p,
                                   const double *a, const double *x,
                                   double *y) {
    for (unsigned s = 1; s <= p; s++) {
        for (unsigned i = 1; i <= n; i++) {
            double sum = 0.0;
            for (unsigned j = 1; j <= n; j++) {
                sum += a[UINDEX(i, j, n)] * x[MINDEX(j, s, n)];
            }
            y[MINDEX(i, s, n)] = sum;
        }
    }
    return;
}

void smacofDistance(const unsigned n, const unsigned p, const double *x,
                    double *d) {
    for (unsigned j = 1; j < (n - 1); j++) {
        for (unsigned i = (j + 1); i <= n; i++) {
            double sum = 0.0;
            for (unsigned s = 1; s <= p; s++) {
                sum += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            double dij = sqrt(fabs(sum));
            d[SINDEX(i, j, n)] = dij;
        }
    }
    return;
}
