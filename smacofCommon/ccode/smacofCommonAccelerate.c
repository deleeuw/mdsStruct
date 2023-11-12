#include "../../smacofInclude/smacof.h"

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

void smacofRate(void) {}
