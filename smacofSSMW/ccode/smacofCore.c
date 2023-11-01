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

void smacofMakeBMatrix(const double *delta, const double *weights,
                       const double *dold, double *bmat, const int *pm) {
    int m = *pm;
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k);
        double dfix = dold[kv];
        if (dfix < 1e-15) {
            bmat[kv] = 0.0;
        } else {
            bmat[kv] = -weights[kv] * delta[kv] / dfix;
        }
    }
    return;
}

void smacofMakeVMatrix(const double *weights, double *v, const int *pm) {
    int m = *pm;
    for (int k = 1; k <= m; k++) {
        v[VINDEX(k)] = -weights[VINDEX(k)];
    }
    return;
}

void smacofGuttman(const double *vinv, const double *bmat, const double *xold,
                   double *xnew, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    double *ymat = (double *)calloc((size_t)np, (size_t)sizeof(double));
    (void)smacofMultiplySDCLMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySDCLMatrix(vinv, ymat, xnew, pn, pp);
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