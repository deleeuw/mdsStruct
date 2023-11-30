#include "smacofUnweighted.h"

void smacofUnweightedMakeBMatrix(const unsigned n, const double **delta,
                                 const double **dmat, double **bmat) {
    double cell = 0.0, dn = (double)n, fac = 2.0 / (dn * (dn - 1.0));
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            double dfix = dmat[i][j];
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = fac * (delta[i][j] / dfix);
            }
            bmat[i][j] = -cell;
            bmat[j][i] = -cell;
        }
    }
    for (unsigned i = 0; i < n; i++) {
        double sum = 0.0;
        for (unsigned j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sum -= bmat[i][j];
        }
        bmat[i][i] = sum;
    }
    return;
}

void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double **bmat, const double **xold,
                             double **xnew) {
    (void)smacofMultiplyAnyAnyMatrix(n, n, p, bmat, xold, xnew);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            xnew[i][j] *= (double)((n - 1) / 2.0);
        }
    }
    return;
}

double smacofUnweightedMakeStress(const unsigned n, const double **delta,
                                  const double **dist) {
    double sum = 0.0, stress = 0.0;
    unsigned m = n * (n - 1) / 2;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            sum += SQUARE(delta[i][j] - dist[i][j]);
        }
    }
    stress = sum / (2.0 * (double)m);
    return stress;
}
