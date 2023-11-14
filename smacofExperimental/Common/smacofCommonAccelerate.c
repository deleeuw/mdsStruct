#include "../Include/smacof.h"

void smacofRelax(const unsigned n, const unsigned p, const double (*xold)[n][p],
                 double (*xnew)[n][p], const double echange, const double pchange,
                 const unsigned itel, const bool relax) {
    double rate = smacofRate(echange, pchange, itel);
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (unsigned i = 0; i < n; i++) {
        for (unsigned s = 0; s < p; s++) {
            (*xnew)[i][s] = (1 + eopt) * (*xnew)[i][s] - eopt * (*xold)[i][s];
        }
    }
}

double smacofRate(const double echange, const double pchange,
                  const unsigned itel) {
    double rate = 0.0;
    if (itel == 1) {
        rate = NaN;
    } else {
        rate = echange / pchange;
    }
    return rate;
}
