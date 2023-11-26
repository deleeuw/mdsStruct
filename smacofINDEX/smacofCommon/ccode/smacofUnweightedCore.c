#include "../../smacofInclude/smacof.h"

void smacofUnweightedMakeBMatrix(const double *delta, const double *dmat,
                                 double *bmat, const int *pn) {
    int n = *pn, nn = n * (n + 1) / 2;
    for (int k = 1; k <= nn; k++) {
        bmat[VINDEX(k)] = -0.0;
    }
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            int ij = SINDEX(i, j, n);
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

void smacofUnweightedGuttman(const double *bmat, const double *xold,
                             double *xnew, const int *pn, const int *pp) {
    (void)smacofMultiplySymmetricMatrix(bmat, xold, xnew, pn, pp);
    return;
}

void smacofUnweightedMakeStress(const double *delta, const double *dist,
                                double *stress, const int *pn) {
    int n = *pn, m = n * (n - 1) / 2;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int ik = VINDEX(k);
        sum += SQUARE(delta[ik] - dist[ik]);
    }
    *stress = sum / (2.0 * (double)m);
    return;
}
