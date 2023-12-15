#include "smacof.h"

void smacofNormDelta(const int n, const bool weights, double **delta,
                     double **w) {
    double sum = 0.0, mw = 2.0 / (double)(n * (n - 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double fac = SQUARE(delta[i][j]);
            if (weights) {
                fac *= w[i][j];
            } else {
                fac *= mw;
            }
            sum += fac;
        }
    }
    sum = 1.0 / sqrt(sum);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            delta[i][j] *= sum;
        }
    }
    return;
}

void smacofScale(const int n, const int p, const bool weights, double **delta,
                 double **w, double **dold, double **xold) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double fac1 = dold[i][j] * delta[i][j];
            double fac2 = SQUARE(dold[i][j]);
            if (weights) {
                fac1 *= w[i][j];
                fac2 *= w[i][j];
            }
            swde += fac1;
            swdd += fac2;
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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            sum += w[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
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
