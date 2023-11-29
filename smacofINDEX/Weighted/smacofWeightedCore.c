#include "../Include/smacof.h"

// bmat always has to be a complete SDCL matrix

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dmat, double *bmat,
                               const int *irow, const int *icol, const int *pn,
                               const int *pm) {
    int m = *pm, n = *pn, mm = n * (n + 1) / 2;
    for (int k = 1; k <= mm; k++) {
        bmat[VINDEX(k)] = -0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double dfix = dmat[SINDEX(ik, jk, n)], cell = 0.0;
        if (dfix < 1e-15) {
            cell = -0.0;
        } else {
            cell = weights[kv] * delta[kv] / dfix;
        }
        bmat[TINDEX(ik, jk, n)] -= cell;
        bmat[TINDEX(ik, ik, n)] += cell;
        bmat[TINDEX(jk, jk, n)] += cell;
    }
    return;
}

// vmat always has to be a complete SDCL matrix

void smacofWeightedMakeVMatrix(const double *weights, double *vmat,
                               const int *irow, const int *icol, const int *pn,
                               const int *pm) {
    int m = *pm, n = *pn, nn = n * (n - 1) / 2;
    for (int k = 1; k <= nn; k++) {
        vmat[VINDEX(k)] = -0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double cell = weights[kv];
        vmat[TINDEX(ik, jk, n)] -= cell;
        vmat[TINDEX(ik, ik, n)] += cell;
        vmat[TINDEX(jk, jk, n)] += cell;
    }
    return;
}

void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const int *pn,
                           const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    double *ymat = (double *)calloc((size_t)np, (size_t)sizeof(double));
    (void)smacofMultiplySymmetricMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySymmetricMatrix(vinv, ymat, xnew, pn, pp);
    free(ymat);
    return;
}

void smacofWeightedMakeStress(const double *delta, const double *weights,
                              const double *dist, const int *irow,
                              const int *icol, const int *pn, const int *pm,
                              double *stress) {
    int m = *pm, n = *pn;
    ;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int kk = VINDEX(k), ik = irow[kk], jk = icol[kk];
        sum += weights[kk] * SQUARE(delta[kk] - dist[SINDEX(ik, jk, n)]);
    }
    *stress = sum / 2.0;
    return;
}
