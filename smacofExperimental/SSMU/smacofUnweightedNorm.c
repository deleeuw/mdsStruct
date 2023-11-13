#include "../Include/smacof.h"

void smacofUnweightedNormDelta(const size_t n, double (*delta)[n]) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            sum += SQUARE(delta[i][j]);
        }
    }
    sum = 1.0 / sqrt(sum);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            delta[i][j] *= sum * sqrt((double)(n * (n - 1)));
        }
    }
    return;
}

void smacofUnweightedScale(const size_t n, const size_t p,
                           const double (*delta)[n], double (*dold)[n],
                           double (*xold)[n]) {
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            swde += dold[i][j] * delta[i][j];
            swdd += SQUARE(dold[i][j]);
        }
    }
    lbd = swde / swdd;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            dold[i][j] *= lbd;
        }
    }
    for (size_t i = 0; i < n; i++) {
        for (size_t s = 0; s < p; s++) {
            xold[i][s] *= lbd;
        }
    }
    return;
}
