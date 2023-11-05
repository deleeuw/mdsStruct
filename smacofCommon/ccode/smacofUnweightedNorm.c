#include "../../smacofInclude/smacof.h"

void smacofUnweightedNormDelta(double *delta, const int *pm) {
    int m = *pm;
    double sum = 0.0;
    for (int i = 1; i <= m; i++) {
        sum += SQUARE(delta[VINDEX(i)]);
    }
    sum = 1 / sqrt(sum);
    for (int i = 1; i <= m; i++) {
        delta[VINDEX(i)] *= sum;
    }
    return;
}

void smacofUnweightedScale(const double *delta, double *dold,
                 double *xold, const int *pn, const int *pp, const int *pm) {
    int n = *pn, p = *pp, m = *pm;
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k);
        swde += dold[kv] * delta[kv];
        swdd += SQUARE(dold[kv]);
    }
    lbd = swde / swdd;
    for (int k = 1; k <= m; k++) {
        dold[VINDEX(k)] *= lbd;
    }
    for (int i = 1; i <= (n * p); i++) {
        xold[VINDEX(i)] *= lbd;
    }
    return;
}
