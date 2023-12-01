
#include "smacofUnweighted.h"

void smacofUnweightedMakeBMatrix(const int n, double **delta, double **dmat,
                                 double **bmat) {
    double cell = 0.0, fac = 1.0 / (double)(n * (n - 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            double dfix = dmat[i][j];
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = fac * (delta[i][j] / dfix);
            }
            bmat[i][j] = -cell;
            bmat[j][i] = -cell;
        }
    }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sum -= bmat[i][j];
        }
        bmat[i][i] = sum;
    }
    return;
}

void smacofUnweightedGuttman(const int n, const int p, double **bmat,
                             double **xold, double **xnew) {
    (void)smacofMultiplyAnyAnyMatrix(n, n, p, bmat, xold, xnew);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            xnew[i][j] /= (double)n;
        }
    }
    return;
}

double smacofUnweightedMakeStress(const int n, double **delta, double **dist) {
    double sum = 0.0, stress = 0.0;
    int m = n * (n - 1) / 2;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += SQUARE(delta[i][j] - dist[i][j]);
        }
    }
    stress = sum / (2.0 * (double)m);
    return stress;
}
