#include "smacofWeighted.h"

void smacofWeightedNormWeights(const int n, double **weights) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += weights[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            weights[i][j] /= sum;
        }
    }
    return;
}

void smacofUnweighting(const int n, double wmax, double **delta, double **weights, 
    double **dmat, double **dhat) {
    for (int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            dhat[i][j] = dmat[i][j] + (weights[i][j] / wmax) * (delta[i][j] - dmat[i][j]);
        }
    }
    return;
}


double smacofWeightedMaxWeights(const int n, double **weights) {
    double wmax = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            wmax = MAX(wmax, weights[i][j]);
        }
    }
    return wmax;
}

void smacofWeightedNormDelta(const int n, double **delta, double **weights) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += weights[i][j] * SQUARE(delta[i][j]);
        }
    }
    sum = 1 / sqrt(sum);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            delta[i][j] *= sum;
        }
    }
    return;
}

void smacofWeightedScale(const int n, const int p, double **delta,
                         double **weights, double **dmat, double **xold) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            swde += weights[i][j] * dmat[i][j] * delta[i][j];
            swdd += weights[i][j] * SQUARE(dmat[i][j]);
        }
    }
    lbd = swde / swdd;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dmat[i][j] *= lbd;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            xold[i][j] *= lbd;
        }
    }
    return;
}
