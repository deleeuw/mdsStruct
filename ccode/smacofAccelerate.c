#include "smacof.h"

void smacofRelax(const double *xold, double *xnew, const double *pechange,
                 const double *ppchange, const int *pnp, const int *pitel,
                 const bool *prelax, double *prate) {
    int np = *pnp, itel = *pitel;
    double echange = *pechange, pchange = *ppchange, rate = *prate;
    bool relax = *prelax;
    if (itel == 1) {
        rate = NaN;
    } else {
        rate = echange / pchange;
        if (relax) {
            double eopt = MAX(0, MIN(1, rate / (2 - rate)));
            for (int i = 1; i <= np; i++) {
                int iv = VINDEX(i);
                xnew[iv] = (1 + eopt) * xnew[iv] - eopt * xold[iv];
            }
        }
    }
    *prate = rate;
}

/*
 * int main() {
    int n = 4, p = 2, itmax = 100, eps = 10, width = 15, precision = 10,
        m = n * (n - 1) / 2;
    bool verbose = true;
    double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double d[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double weights[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double x[8] = {1.0, 1.0, -1.0, -1.0, -4.0, 4.0, -4.0, 4.0};
    (void)smacofNormWeights(weights, &m);
    (void)smacofNormDelta(delta, weights, &m);
    (void)smacofDistance(x, d, &n, &p);
    (void)smacofScale(delta, weights, d, x, &n, &p);
    (void)smacofDiagonalAdjust(delta, weights, x, &n, &p, &itmax, &eps,
                               &verbose);
}
*/
