#include "smacof.h"

void smacofDistance(const double *x, double *d, const int *irow, 
    const int *icol, const int *pn, const int *pp, const int *pm) {
    int n = *pn, p = *pp, m = *pm;
    for (int j = 1; j <= m; j++) {  
            int i = irow(VINDEX(j));
            int j = jrow(VINDEX(j));     
            double sum = 0.0;
            for (int s = 1; s <= p; s++) {
                sum += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            d[VINDEX(j)] = sqrt(fabs(sum));
    }
    return;
}

void smacofMakeBMatrix(const double *delta, const double *weights,
                       const double *dold, double *bmat, const int *pm) {
    int m = *pm;
    for (int j = 1; j <= m; j++) {
        int i = VINDEX(j);
        double dfix = dold[i];
        if (dfix < 1e-15) {
            bmat[i] = 0.0;
        } else {
            bmat[i] = -weights[i] * delta[i] / dfix;
        }
    }
    return;
}

void smacofMakeVMatrix(const double *weights, double *v, const int *pm) {
    int m = *pm;
    for (int j = 1; j <= m; j++) {
        v[VINDEX(j)] = -weights[VINDEX(jrow)];
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
        int i = VINDEX(k);
        sum += weights[i] * SQUARE(delta[i] - dist[i]);
    }
    *stress = sum / 2.0;
    return;
}

void smacofEtaSquare(const double *weights, const double *dist, const int *pm,
                     double *etasquare) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int i = VINDEX(k);
        sum += weights[i] * SQUARE(dist[i]);
    }
    *etasquare = sum;
    return;
}

void smacofRho(const double *delta, const double *weights, const double *dist,
               const int *pm, double *rho) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int i = VINDEX(k);
        sum += weights[i] * delta[i] * dist[i];
    }
    *rho = sum;
    return;
}