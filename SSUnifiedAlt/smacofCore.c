#include "smacof.h"

/*
double smacofStress(const int n, const bool weights, double **delta, double **w,
                    double **dmat) {
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double fac1 = SQUARE(delta[i][j] - dmat[i][j]);
            double fac2 = SQUARE(dmat[i][j]);
            if (weights) {
                sum1 += w[i][j] * fac1;
                sum2 += w[i][j] * fac2;
            } else {
                sum1 += fac1;
                sum2 += fac2;
            }
        }
    }
    return sum1 / sum2;
}
*/

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

void smacofGuttmanTransformAlt(const int n, const int p, const bool weights,
                               double **dhat, double **dmat, double **w,
                               double **vinv, double **xold, double **xnew) {
    double **xaux = smacofMakeAnyMatrix(n, p);
    for (int s = 0; s < p; s++) {
        for (int i = 1; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                double rat = (dhat[i][j] / dmat[i][j]);
                if (weights) {
                    rat *= w[i][j];
                }
                sum += rat * (xold[i][s] - xold[j][s]);
            }
            xaux[i][s] = sum;
        }
    }
    if (weights) {
        for (int s = 0; s < p; s++) {
            for (int i = 1; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    double rat = vinv[i][j];
                    sum += rat * (xaux[i][s] - xaux[j][s]);
                }
                xnew[i][s] = sum;
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
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double wij = weights[i][j];
            vmat[i][j] = -wij;
            vmat[i][i] += wij;
            vmat[j][j] += wij;
        }
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
