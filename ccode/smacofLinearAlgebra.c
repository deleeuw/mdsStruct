#include "smacof.h"

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
                    int ik = IMIN(i, k);
                    int ki = IMAX(i, k);
                    int jk = IMIN(j, k);
                    int kj = IMAX(j, k);
                    oldi[VINDEX(k)] = a[TINDEX(ki, ik, n)];
                    oldj[VINDEX(k)] = a[TINDEX(kj, jk, n)];
                }
                for (int k = 1; k <= n; k++) {
                    int ik = IMIN(i, k);
                    int ki = IMAX(i, k);
                    int jk = IMIN(j, k);
                    int kj = IMAX(j, k);
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
        if (((fnew - fold) < eps) || (itel == itmax)) break;
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


/*
    int main() {
    double a[10] = {758.0, 504.0, 261.0, 510.0, 374.0,
                    193.0, 382.0, 115.0, 212.0, 427.0};
    double evec[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double eval[4] = {0.0, 0.0, 0.0, 0.0};
    int n = 4, m = 2, itmax = 100, width = 15, precision = 10, nr = 1;
    int eps = 10;
    bool verbose = true;
    (void)smacofJacobi(a, evec, eval, &n, &m, &itmax, &eps, &verbose);
    (void)smacofPrintLTMatrix(a, &n, &width, &precision);
    (void)smacofPrintAnyMatrix(evec, &n, &n, &width, &precision);
    (void)smacofPrintAnyMatrix(eval, &nr, &n, &width, &precision);
}
*/