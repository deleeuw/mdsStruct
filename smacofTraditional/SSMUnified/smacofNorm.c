#include "smacof.h"

void smacofNormDelta(const int n, const bool weights, double **delta,
                     double **w) {
    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < (n - 1); j++) {
            double fac = SQUARE(delta[i][j]);
            if (weights) {
                fac *= w[i][j];
            }
            sum += fac;
        }
    }
    sum = 1.0 / sqrt(sum);
    if (!weights) {
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < (n - 1); j++) {
                delta[i][j] *= sum * sqrt((double)(n * (n - 1)));
            }
        }
    }
    return;
}

void smacofScale(const int n, const int p, const bool weights, double **delta,
                 double **w, double **dold, double **xold) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            swde += dold[i][j] * delta[i][j];
            swdd += SQUARE(dold[i][j]);
        }
    }
    lbd = swde / swdd;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dold[i][j] *= lbd;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xold[i][s] *= lbd;
        }
    }
    return;
}

void smacofNormWeights(const int n, double **w) {
    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < (n - 1); j++) {
            sum += w[i][j];
        }
    }
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < (n - 1); j++) {
            w[i][j] /= sum;
        }
    }
    return;
}

void smacofUnweighting(const int n, double wmax, double **delta, double **w,
                       double **dmat, double **dhat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dhat[i][j] =
                dmat[i][j] + (w[i][j] / wmax) * (delta[i][j] - dmat[i][j]);
        }
    }
    return;
}

double smacofMaxWeights(const int n, double **w) {
    double wmax = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            wmax = MAX(wmax, w[i][j]);
        }
    }
    return wmax;
}
