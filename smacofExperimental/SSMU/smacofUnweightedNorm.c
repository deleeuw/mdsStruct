#include "../Include/smacof.h"

void smacofUnweightedNormDelta(const unsigned n, double delta[n][n]) {
    double sum = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            sum += SQUARE(delta[i][j]);
        }
    }
    sum = 1.0 / sqrt(sum);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            delta[i][j] *= sum * sqrt((double)(n * (n - 1)));
        }
    }
    return;
}

void smacofUnweightedScale(const unsigned n, const unsigned p,
                           const double delta[n][n], double dold[n][n],
                           double xold[n][p]) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            swde += dold[i][j] * delta[i][j];
            swdd += SQUARE(dold[i][j]);
        }
    }
    lbd = swde / swdd;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            dold[i][j] *= lbd;
        }
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned s = 0; s < p; s++) {
            xold[i][s] *= lbd;
        }
    }
    return;
}
