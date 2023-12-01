#include "smacofCommon.h"

void smacofRelax(const int n, const int p, const double rate, double **cxold,
                 double **cxnew) {
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            cxnew[i][s] = (1 + eopt) * cxnew[i][s] - eopt * cxold[i][s];
        }
    }
}
