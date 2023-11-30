#include "smacofCommon.h"

void smacofRelax(const unsigned n, const unsigned p, const double *xold,
                 double *xnew, double rate) {
    unsigned np = n * p;
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (unsigned i = 1; i <= np; i++) {
        xnew[VINDEX(i)] = (1 + eopt) * xnew[VINDEX(i)] - eopt * xold[VINDEX(i)];
    }
}

/* to do: split rate and relax */
