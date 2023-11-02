#include "smacof.h"

void smacofNormWeights(double *weights, const int *pm) {
    int m = *pm;
    double sum = 0.0;
    for (int i = 1; i <= m; i++) {
        sum += weights[VINDEX(i)];
    }
    for (int i = 1; i <= m; i++) {
        weights[VINDEX(i)] /= sum;
    }
    return;
}

void smacofNormDelta(double *delta, const double *weights, const int *pm) {
    int m = *pm;
    double sum = 0.0;
    for (int i = 1; i <= m; i++) {
        sum += weights[VINDEX(i)] * SQUARE(delta[VINDEX(i)]);
    }
    sum = 1 / sqrt(sum);
    for (int i = 1; i <= m; i++) {
        delta[VINDEX(i)] *= sum;
    }
    return;
}

void smacofCenter(double *x, const int *np, const int *pp) {
    int n = *np, p = *pp;
    for (int s = 1; s <= p; s++) {
        double sum = 0.0;
        for (int i = 1; i <= n; i++) {
            sum += x[MINDEX(i, s, n)];
        }
        sum /= (double)n;
        for (int i = 1; i <= n; i++) {
            x[MINDEX(i, s, n)] -= sum;
        }
    }
    return;
}

void smacofScale(const double *delta, const double *weights, double *dold,
                 double *xold, const int *pn, const int *pp, const int *pm) {
    int n = *pn, p = *pp, m = *pm;
    double swde = 0.0, swdd = 0.0, lbd = 0.0;
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k);
        swde += weights[kv] * dold[kv] * delta[kv];
        swdd += weights[kv] * SQUARE(dold[kv]);
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
