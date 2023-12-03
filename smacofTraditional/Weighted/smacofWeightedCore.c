#include "smacofWeighted.h"


void smacofWeightedMakeBMatrix(cont int n, double **delta,
                               double **weights, double **dmat, double **bmat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
        double dfix = dmat[i][j]], cell = 0.0;
        if (dfix < 1e-15) {
            cell = -0.0;
        } else {
            cell = weights[i][j] * delta[i][j] / dfix;
        }
        bmat[i][j] -= cell;
        bmat[j][i] -= cell;
        bmat[i][i] += cell;
        bmat[j][j] += cell;
        }
    }
    return;
}


void smacofWeightedMakeVMatrix(const int n, double **weights,
                               double **vmat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double cell = weights[i][j];
            vmat[i][j] -= cell;
            vmat[j][i] -= cell;
            vmat[i][i] += cell;
            vmat[j][j] += cell;
        }
    }
    return;
}

void smacofWeightedGuttman(const int n, const int p, double **vinv, double **bmat,
                           double **xold, double **xnew) {
    double **ymat = smacofMakeAnyMatrix(n, p);
    (void)smacofMultiplyAnyAnyMatrix(n, p, bmat, xold, ymat);
    (void)smacofMultiplyAnyAnyMatrix(n, p, vinv, ymat, xnew);
    (void)smacofFreeAnyMatrix(ymat);
    return;
}

double smacofWeightedMakeStress(const int n, double **delta, double **weights,
                              double **dist) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
        sum += weights[i][j] * SQUARE(delta[i][j] - dist[i][j]);
    }
    return sum / 2.0;
}
