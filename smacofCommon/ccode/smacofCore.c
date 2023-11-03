#include "smacof.h"

void smacofDistance(const double *x, double *d, const int *irow,
                    const int *icol, const int *pn, const int *pp,
                    const int *pm) {
    int n = *pn, p = *pp, m = *pm;
    for (int k = 1; k <= m; k++) {
        int i = irow[VINDEX(k)];
        int j = icol[VINDEX(k)];
        double sum = 0.0;
        for (int s = 1; s <= p; s++) {
            sum += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
        }
        d[VINDEX(k)] = sqrt(fabs(sum));
    }
    return;
}

// bmat always has to be a complete SDCL matrix

void smacofMakeBMatrix(const double *delta, const double *weights,
                       const double *dold, double *bmat, const int *irow,
                       const int *icol, const int *pn, const int *pm) {
    int m = *pm, n = *pn, nn = n * (n - 1) / 2;
    for (int k = 1; k <= nn; k++) {
        bmat[VINDEX(k)] = -0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double dfix = dold[kv], cell = 0.0;
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

void smacofMakeVMatrix(const double *weights, double *vmat, const int *irow,
                       const int *icol, const int *pn, const int *pm) {
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

void smacofGuttman(const double *vinv, const double *bmat, const double *xold,
                   double *xnew, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    int width = 15, precision = 10;
    double *ymat = (double *)calloc((size_t)np, (size_t)sizeof(double));
    (void)smacofMultiplySymmetricMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySymmetricMatrix(vinv, ymat, xnew, pn, pp);
    free(ymat);
    return;
}

void smacofStress(const double *delta, const double *weights,
                  const double *dist, const int *pm, double *stress) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int ik = VINDEX(k);
        sum += weights[ik] * SQUARE(delta[ik] - dist[ik]);
    }
    *stress = sum / 2.0;
    return;
}

/*
void smacofEtaSquare(const double *weights, const double *dist, const int *pm,
                     double *etasquare) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int ik = VINDEX(k);
        sum += weights[ik] * SQUARE(dist[ik]);
    }
    *etasquare = sum;
    return;
}

void smacofRho(const double *delta, const double *weights, const double *dist,
               const int *pm, double *rho) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int ik = VINDEX(k);
        sum += weights[ik] * delta[ik] * dist[ik];
    }
    *rho = sum;
    return;
}
*/