#include "smacofUnweighted.h"

void smacofUnweightedNormDelta(const int n, double **delta) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += SQUARE(delta[i][j]);
        }
    }
    sum = 1.0 / sqrt(sum);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            delta[i][j] *= sum * sqrt((double)(n * (n - 1)));
        }
    }
    return;
}

void smacofUnweightedScale(const int n, const int p, double **delta,
                           double **dold, double **xold) {
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
