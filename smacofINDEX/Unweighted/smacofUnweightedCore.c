#include "smacofUnweighted.h"

void smacofUnweightedMakeBMatrix(const unsigned n, const double *delta,
                                 const double *dmat, double *bmat) {
    unsigned nn = n * (n + 1) / 2;
    for (unsigned k = 1; k <= nn; k++) {
        bmat[VINDEX(k)] = -0.0;
    }
    for (unsigned j = 1; j <= (n - 1); j++) {
        for (unsigned i = (j + 1); i <= n; i++) {
            unsigned ij = SINDEX(i, j, n);
            double dfix = dmat[ij], cell = 0.0;
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = delta[ij] / (dfix * (double)n);
            }
            bmat[TINDEX(i, j, n)] -= cell;
            bmat[TINDEX(i, i, n)] += cell;
            bmat[TINDEX(j, j, n)] += cell;
        }
    }
    return;
}

void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double *bmat, const double *xold,
                             double *xnew) {
    (void)smacofMultiplySymmetricMatrix(n, p, bmat, xold, xnew);
    return;
}

double smacofUnweightedMakeStress(const unsigned n, const double *delta,
                                  const double *dist) {
    unsigned m = n * (n - 1) / 2;
    double sum = 0.0;
    for (unsigned k = 1; k <= m; k++) {
        unsigned ik = VINDEX(k);
        sum += SQUARE(delta[ik] - dist[ik]);
    }
    return sum / (2.0 * (double)m);
}
