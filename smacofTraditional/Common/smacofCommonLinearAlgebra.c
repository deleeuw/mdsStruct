#include "smacofCommon.h"

void smacofJacobi(const int n, const int ndim, double **a, double **evec,
                  double *eval, const int itmax, const int ieps,
                  const bool verbose) {
    int itel = 1;
    double d = 0.0, s = 0.0, t = 0.0, u = 0.0, v = 0.0, p = 0.0, q = 0.0,
           r = 0.0;
    double fold = 0.0, fnew = 0.0, eps = pow(10.0, -(double)ieps);
    double *oldi = smacofMakeAnyVector(n);
    double *oldj = smacofMakeAnyVector(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            evec[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < ndim; i++) {
        fold += SQUARE(a[i][i]);
    }
    while (true) {
        for (int j = 0; j < ndim; j++) {
            for (int i = j + 1; i < n; i++) {
                p = a[i][j];
                q = a[i][i];
                r = a[j][j];
                if (fabs(p) < 1e-10) {
                    continue;
                }
                d = (q - r) / 2.0;
                s = (p < 0) ? -1.0 : 1.0;
                t = -d / sqrt(SQUARE(d) + SQUARE(p));
                u = sqrt((1 + t) / 2);
                v = s * sqrt((1 - t) / 2);
                for (int k = 0; k < n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    oldi[k] = a[ki][ik];
                    oldj[k] = a[kj][jk];
                }
                for (int k = 0; k < n; k++) {
                    int ik = MIN(i, k);
                    int ki = MAX(i, k);
                    int jk = MIN(j, k);
                    int kj = MAX(j, k);
                    a[ki][ik] = u * oldi[k] - v * oldj[k];
                    a[kj][jk] = v * oldi[k] + u * oldj[k];
                }
                for (int k = 0; k < n; k++) {
                    oldi[k] = evec[k][i];
                    oldj[k] = evec[k][j];
                    evec[k][i] = u * oldi[k] - v * oldj[k];
                    evec[k][j] = v * oldi[k] + u * oldj[k];
                }
                a[i][i] = SQUARE(u) * q + SQUARE(v) * r - 2 * u * v * p;
                a[j][j] = SQUARE(v) * q + SQUARE(u) * r + 2 * u * v * p;
                a[i][j] = u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
                a[j][i] = u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
            }
        }
        fnew = 0.0;
        for (int i = 0; i < ndim; i++) {
            fnew += SQUARE(a[i][i]);
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
    for (int i = 0; i < n; i++) {
        eval[i] = a[i][i];
    }
    (void)smacofFreeAnyVector(oldi);
    (void)smacofFreeAnyVector(oldj);
    return;
}

// for now without pivoting

void smacofGramSchmidt(const int n, const int p, double **x, double **q) {
    for (int s = 0; s < p; s++) {
        if (s > 0) {
            for (int t = 0; t < s; t++) {
                double sum = 0.0;
                for (int i = 0; i < n; i++) {
                    sum += x[i][s] * x[i][t];
                }
                q[t][s] = sum;
                for (int i = 0; i < n; i++) {
                    x[i][s] -= sum * x[i][t];
                }
            }
        }
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += SQUARE(x[i][s]);
        }
        sum = sqrt(sum);
        q[s][s] = sum;
        for (int i = 0; i < n; i++) {
            x[i][s] /= sum;
        }
    }
}

void smacofCenter(const int n, const int p, double **x) {
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

/*
void smacofInvertPDMatrix(const double *x, double *xinv, const int *pn) {
    int n = *pn, m = n * (n + 1) / 2, ik = 0, jk = 0, ij = 0;
    for (int k = 1; k <= m; k++) {
        xinv[VINDEX(k)] = x[VINDEX(k)];
    }
    for (int k = 1; k <= n; k++) {
        double piv = xinv[TINDEX(k, k, n)];
        for (int j = 1; j <= n; j++) {
            if (j == k) {
                continue;
            }
            jk = UINDEX(j, k, n);
            for (int i = j; i <= n; i++) {
                if (i == k) {
                    continue;
                }
                ik = UINDEX(i, k, n);
                ij = TINDEX(i, j, n);
                xinv[ij] = xinv[ij] - xinv[ik] * xinv[jk] / piv;
            }
        }
        for (int i = 1; i <= n; i++) {
            if (i == k) {
                continue;
            }
            ik = UINDEX(i, k, n);
            xinv[ik] = xinv[ik] / piv;
        }
        xinv[TINDEX(k, k, n)] = -1 / piv;
    }
    for (int k = 1; k <= m; k++) {
        xinv[VINDEX(k)] = -xinv[VINDEX(k)];
    }
    return;
}

void smacofMPInverseSDCLMatrix(const double *vmat, double *vinv,
                               const int *pn) {
    int n = *pn, nn = n * (n + 1) / 2;
    double add = 1.0 / ((double)n);
    double *vadd = malloc((sizeof *vadd) * nn);
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            vadd[TINDEX(i, j, n)] = vmat[TINDEX(i, j, n)] + add;
        }
    }
    (void)smacofInvertPDMatrix(vadd, vinv, pn);
    for (int k = 1; k <= nn; k++) {
        vinv[VINDEX(k)] -= add;
    }
    free(vadd);
    return;
}
*/

void smacofMultiplyAnyAnyMatrix(const int n, const int p, const int m,
                                double **a, double **x, double **y) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int s = 0; s < p; s++) {
                sum += a[i][s] * x[s][j];
            }
            y[i][j] = sum;
        }
    }
    return;
}

void smacofDistance(const int n, const int p, double **x, double **d) {
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int s = 0; s < p; s++) {
                sum += SQUARE(x[i][s] - x[j][s]);
            }
            d[i][j] = sqrt(fabs(sum));
        }
    }
    return;
}

/*
int main(void) {
    int n = 5;
    int p = 2;
    int width = 15;
    int precision = 10;
    double rr[15] = {4, -1, -1, -1, -1, 5, -1, -1, -1, 6, -1, -1, 7, -1, 8};
    double rx[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 1.0, 2.0, 4.0, 3.0};
    double(*cr)[n][n] = malloc(sizeof(*cr));
    assert(!(cr == NULL));
    double(*cx)[n][p] = malloc(sizeof(*cx));
    assert(!(cx == NULL));
    double(*cy)[n][p] = malloc(sizeof(*cy));
    assert(!(cy == NULL));
    double(*evec)[n][n] = malloc(sizeof(*evec));
    assert(!(evec == NULL));
    double(*cd)[p][p] = malloc(sizeof(*cd));
    assert(!(cd == NULL));
    double eval[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    (void)smacofFromSymmetricRtoC(n, rr, cr);
    (void)smacofPrintAnyMatrix(n, n, width, precision, cr);
    (void)smacofFromAnyRtoC(n, p, rx, cx);
    (void)smacofPrintAnyMatrix(n, p, width, precision, cx);
    (void)smacofCenter(n, p, cx);
    (void)smacofPrintAnyMatrix(n, p, width, precision, cx);
    (void)smacofGramSchmidt(n, p, cx, cd);
    (void)smacofPrintAnyMatrix(n, p, width, precision, cx);
    (void)smacofPrintAnyMatrix(p, p, width, precision, cd);
    (void)smacofMultiplyAnyAnyMatrix(n, p, cr, cx, cy);
    (void)smacofPrintAnyMatrix(n, p, width, precision, cy);
    (void)smacofJacobi(n, p, cr, evec, eval, 100, 10, true);
    (void)smacofPrintAnyMatrix(n, n, width, precision, evec);
    for (int i = 0; i < p; i++) {
        printf("%10.6f ", eval[i]);
    }
    printf("\n");
    return EXIT_SUCCESS;
}
*/
