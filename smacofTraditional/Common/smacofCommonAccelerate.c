#include "smacofCommon.h"

void smacofRelax(const unsigned n, const unsigned p,
                 const double **cxold, double **cxnew,
                 double rate) {
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (unsigned i = 0; i < n; i++) {
        for (unsigned s = 0; s < p; s++) {
            (*cxnew)[i][s] =
                (1 + eopt) * cxnew[i][s] - eopt * cxold[i][s];
        }
    }
}
