#include "../Include/smacof.h"

void smacofRelax(const size_t n, const size_t p, const double (*xold)[n],
                 double (*xnew)[n], const double echange, const double pchange,
                 const size_t itel, const bool relax) {
    double rate = smacofRate(echange, pchange, itel);
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (size_t i = 0; i < n; i++) {
        for (size_t s = 0; s < p; s++) {
            xnew[i][s] = (1 + eopt) * xnew[i][s] - eopt * xold[i][s];
        }
    }
}

double smacofRate(const double echange, const double pchange,
                  const size_t itel) {
    double rate = 0.0;
    if (itel == 1) {
        rate = NaN;
    } else {
        rate = echange / pchange;
    }
    return rate;
}
