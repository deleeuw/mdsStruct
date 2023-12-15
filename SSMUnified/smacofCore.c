#include "smacof.h"

double smacofStress(const int n, const bool weights, double **delta, double **w,
                    double **dmat) {
    double sum = 0.0, mw = 2.0 / (double)(n * (n - 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double fac = SQUARE(delta[i][j] - dmat[i][j]);
            if (weights) {
                fac *= w[i][j];
            } else {
                fac *= mw;
            }
            sum += fac;
        }
    }
    return sum;
}

void smacofGuttmanTransform(const int n, const int p, const bool weights,
                            double **dhat, double **dmat, double **w,
                            double **vinv, double **xold, double **xnew) {
    double **xaux = smacofMakeAnyMatrix(n, p);
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                int ij = MAX(i, j);
                int ji = MIN(i, j);
                double rat = (dhat[ij][ji] / dmat[ij][ji]);
                if (weights) {
                    rat *= w[ij][ji];
                }
                sum += rat * (xold[i][s] - xold[j][s]);
            }
            xaux[i][s] = sum;
        }
    }
    if (weights) {
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        continue;
                    }
                    int ij = MAX(i, j);
                    int ji = MIN(i, j);
                    double rat = vinv[ij][ji];
                    sum += rat * (xaux[i][s] - xaux[j][s]);
                }
                xnew[i][s] = -sum;
            }
        }
    } else {
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < n; i++) {
                xnew[i][s] = xaux[i][s] / (double)n;
            }
        }
    }
    free(xaux);
    return;
}

void smacofMakeVMatrix(const int n, double **weights, double **vmat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            vmat[i][j] = -weights[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            sum += vmat[ij][ji];
        }
        vmat[i][i] = -sum;
    }
    return;
}

void smacofInverseVMatrix(const int n, double **vmat, double **vinv) {
    double add = 1.0 / ((double)n);
    double **vadd = smacofMakeAnyMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            vadd[i][j] = vmat[i][j] + add;
        }
    }
    (void)smacofInvertPositiveDefiniteMatrix(n, vadd, vinv);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            vinv[i][j] -= add;
        }
    }
    smacofFreeMatrix(n, vadd);
    return;
}
