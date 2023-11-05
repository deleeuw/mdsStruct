#include "../../smacofInclude/smacof.h"

void smacofUnweightedMakeBMatrix(const double *delta, const double *dmat,
                                 double *bmat, const int *irow, const int *icol,
                                 const int *pn, const int *pm) {
    int m = *pm, n = *pn, nn = n * (n - 1) / 2, mm = n * (n + 1) / 2;
    for (int k = 1; k <= mm; k++) {
        bmat[VINDEX(k)] = -0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double dfix = dmat[SINDEX(ik, jk, n)], cell = 0.0;
        if (dfix < 1e-15) {
            cell = -0.0;
        } else {
            cell = delta[kv] / dfix;
        }
        bmat[TINDEX(ik, jk, n)] -= cell;
        bmat[TINDEX(ik, ik, n)] += cell;
        bmat[TINDEX(jk, jk, n)] += cell;
    }
    return;
}

void smacofUnweightedGuttman(const double *bmat, const double *xold,
                             double *xnew, const int *pn, const int *pp) {
    int n = *pn, p = *pp;
    (void)smacofMultiplySymmetricMatrix(bmat, xold, xnew, pn, pp);
    return;
}

void smacofUnweightedMakeStress(const double *delta, const double *dist,
                                const int *pm, double *stress) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int ik = VINDEX(k);
        sum += SQUARE(delta[ik] - dist[ik]);
    }
    *stress = sum / 2.0;
    return;
}
