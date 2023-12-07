#include "smacofWeighted.h"

void smacofWeightedMakeBMatrix(const int n, double **delta, double **weights,
                               double **dmat, double **bmat) {
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double dfix = dmat[i][j], cell = 0.0;
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = weights[i][j] * delta[i][j] / dfix;
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
            sum += bmat[i][j];
        }
        bmat[i][i] = -sum;
    }
    return;
}

void smacofWeightedMakeVMatrix(const int n, double **weights, double **vmat) {
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double cell = weights[i][j];
            vmat[i][j] -= cell;
            vmat[j][i] -= cell;
            vmat[i][i] += cell;
            vmat[j][j] += cell;
        }
    }
    return;
}

void smacofWeightedInverseVMatrix(const int n, double **vmat, double **vinv) {
    double add = 1.0 / ((double)n);
    double **vadd = smacofMakeAnyMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vadd[i][j] = vmat[i][j] + add;
        }
    }
    (void)smacofInvertPositiveDefiniteMatrix(n, vadd, vinv);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vinv[i][j] -= add;
        }
    }
    smacofFreeAnyMatrix(n, vadd);
    return;
}

void smacofWeightedGuttman(const int n, const int p, double **vinv,
                           double **bmat, double **xold, double **xnew) {
    double **ymat = smacofMakeAnyMatrix(n, p);
    (void)smacofMultiplyAnyAnyMatrix(n, n, p, bmat, xold, ymat);
    (void)smacofMultiplyAnyAnyMatrix(n, n, p, vinv, ymat, xnew);
    (void)smacofFreeAnyMatrix(n, ymat);
    return;
}

double smacofWeightedMakeStress(const int n, double **delta, double **weights,
                                double **dist) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += weights[i][j] * SQUARE(delta[i][j] - dist[i][j]);
        }
    }
    return sum;
}
