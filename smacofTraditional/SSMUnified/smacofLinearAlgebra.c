#include "smacofCommon.h"

void smacofJacobi(const int n, const int m, double **a, double **evec,
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
    for (int i = 0; i < m; i++) {
        fold += SQUARE(a[i][i]);
    }
    while (true) {
        for (int j = 0; j < m; j++) {
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
                    oldi[k] = a[i][k];
                    oldj[k] = a[j][k];
                }
                for (int k = 0; k < n; k++) {
                    a[i][k] = u * oldi[k] - v * oldj[k];
                    a[j][k] = v * oldi[k] + u * oldj[k];
                    a[k][i] = a[i][k];
                    a[k][j] = a[j][k];
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
                a[j][i] = a[i][j];
            }
        }
        fnew = 0.0;
        for (int i = 0; i < m; i++) {
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
    for (int i = 0; i < m; i++) {
        eval[i] = a[i][i];
    }
    (void)smacofFreeAnyVector(oldi);
    (void)smacofFreeAnyVector(oldj);
    return;
}

/*
void smacofDoubleJacobi(const int n, double **a, double **b, double **evec,
                        double *eval, double **bvec, double *bval,
                        const int itmax, const int ieps, const bool verbose) {
    (void)smacofJacobi(n, n, b, evec, eval, itmax, ieps, verbose);
    (void)smacofCopyAnyMatrix(n, n, evec, bvec);
    (void)smacofCopyAnyVector(n, eval, bval);
    double **h = smacofMakeAnyMatrix(n, n);
    (void)smacofMultiplyAnyAnyMatrix(n, n, n, a, evec, h);
    (void)smacofCrossprodAnyAnyMatrix(n, n, n, evec, h, a);
    // X'AX with X=K\Lambda^{-1/2}
    double *div = smacofMakeAnyVector(n);
    for (int i = 0; i < n; i++) {
        div[i] = sqrt(fabs(bval[i]));
        if (div[i] < 1e-15) {
            div[i] = 0.0;
        } else {
            div[i] = 1.0 / div[i];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] *= (div[i] * div[j]);
        }
    }
    (void)smacofJacobi(n, n, a, evec, eval, itmax, ieps, false);
    (void)smacofPrintAnyVector(n, 15, 10, eval);
    (void)smacofFreeAnyMatrix(n, h);
    (void)smacofFreeAnyVector(div);
}
*/

void smacofMultipleAnySymmetricAnyMatrix(const int n, const int m, double **x,
                                         double **a, double **u) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    sum += x[k][i] * x[l][j] * a[k][l];
                }
            }
            u[i][j] = sum;
        }
    }
}

void smacofScaleMatrixColumns(const int n, const int m, const double p,
                              double **x, double *y, double **v) {
    for (int j = 0; j < m; j++) {
        if (y[j] < 1e-15) {
            y[j] = 0.0;
        } else {
            y[j] = pow(y[j], p);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            v[i][j] = x[i][j] * y[j];
        }
    }
}

void smacofDoubleJacobi(const int n, double **a, double **b, double **evec,
                        double *eval, const int itmax, const int ieps,
                        const bool verbose) {
    double **bvec = smacofMakeAnyMatrix(n, n);
    double *bval = smacofMakeAnyVector(n);
    double **v = smacofMakeAnyMatrix(n, n);
    double **u = smacofMakeAnyMatrix(n, n);
    (void)smacofJacobi(n, n, b, bvec, bval, itmax, ieps, verbose);
    (void)smacofScaleMatrixColumns(n, n, -0.5, bvec, bval, v);
    (void)smacofMultipleAnySymmetricAnyMatrix(n, n, v, a, u);
    (void)smacofJacobi(n, n, u, bvec, eval, itmax, ieps, verbose);
    (void)smacofMultiplyAnyAnyMatrix(n, n, n, v, bvec, evec);
    (void)smacofFreeAnyMatrix(n, bvec);
    (void)smacofFreeAnyMatrix(n, v);
    (void)smacofFreeAnyMatrix(n, u);
    (void)smacofFreeAnyVector(bval);
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

void smacofInvertPositiveDefiniteMatrix(const int n, double **x,
                                        double **xinv) {
    (void)smacofCopyAnyMatrix(n, n, x, xinv);
    for (int k = 0; k < n; k++) {
        double piv = xinv[k][k];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i != k) && (j != k)) {
                    xinv[i][j] -= xinv[i][k] * xinv[k][j] / piv;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            if (i != k) {
                xinv[k][i] /= piv;
                xinv[i][k] /= piv;
            }
        }
        xinv[k][k] = -1.0 / piv;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            xinv[i][j] = -xinv[i][j];
        }
    }
    return;
}

void smacofMultiplyAnyAnyMatrix(const int n, const int p, const int m,
                                double **x, double **y, double **z) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < p; k++) {
                sum += x[i][k] * y[k][j];
            }
            z[i][j] = sum;
        }
    }
    return;
}

void smacofCrossprodAnyAnyMatrix(const int n, const int p, const int m,
                                 double **x, double **y, double **z) {
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < p; k++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += x[i][k] * y[i][j];
            }
            z[k][j] = sum;
        }
    }
    return;
}

void smacofMultiplySymmetricAnyMatrix(const int n, const int p, double **a,
                                      double **x, double **y) {
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                if (k <= i) {
                    sum += a[i][k] * x[k][s];
                } else {
                    sum += a[k][i] * x[k][s];
                }
            }
            y[i][s] = sum;
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
