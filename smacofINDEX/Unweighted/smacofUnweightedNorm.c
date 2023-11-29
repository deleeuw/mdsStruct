#include "smacofUnweighted.h"

void smacofUnweightedNormDelta(const unsigned n, double *delta) {
    unsigned m = n * (n - 1) / 2;
    double sum = 0.0;
    for (unsigned i = 1; i <= m; i++) {
        sum += SQUARE(delta[VINDEX(i)]);
    }
    sum = 1 / sqrt(sum);
    for (unsigned i = 1; i <= m; i++) {
        delta[VINDEX(i)] *= sum * sqrt((double)m);
    }
    return;
}

void smacofUnweightedScale(const unsigned n, const unsigned p,
                           const double *delta, double *dold, double *xold) {
    unsigned m = n * (n - 1) / 2, np = n * p;
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (unsigned k = 1; k <= m; k++) {
        unsigned kv = VINDEX(k);
        swde += dold[kv] * delta[kv];
        swdd += SQUARE(dold[kv]);
    }
    lbd = swde / swdd;
    for (unsigned k = 1; k <= m; k++) {
        dold[VINDEX(k)] *= lbd;
    }
    for (unsigned i = 1; i <= np; i++) {
        xold[VINDEX(i)] *= lbd;
    }
    return;
}
