#include "smacofCommon.h"

double smacofRelax(const double *xold, double *xnew, const double echange,
                   const double pchange, const unsigned np, const unsigned itel,
                   const bool relax) {
    double rate = 0.0;
    if (itel == 1) {
        rate = NaN;
    } else {
        rate = echange / pchange;
        if (relax) {
            double eopt = MAX(0, MIN(1, rate / (2 - rate)));
            for (unsigned i = 1; i <= np; i++) {
                xnew[i] = (1 + eopt) * xnew[i] - eopt * xold[i];
            }
        }
    }
    return rate;
}

/* to do: split rate and relax */
