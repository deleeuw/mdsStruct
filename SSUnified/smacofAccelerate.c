#include "smacof.h"

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew) {
    double eopt = MAX(0, MIN(1, rate / (2 - rate)));
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xnew[i][s] = (1 + eopt) * xnew[i][s] - eopt * xold[i][s];
        }
    }
}
