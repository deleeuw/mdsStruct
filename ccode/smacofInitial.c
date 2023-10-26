#include "smacof.h"

void smacofInitial(const double *delta, const double *weights, double *xini,
                   double *dini, const int *pinit, const int *pn,
                   const int *pp) {
    int n = *pn, p = *pp, np = n * p, init = *pinit;
    int itmax_j = 100, itmax_e = 100, eps_j = 10;
    int width = 15, precision = 10;
    double eps_e = pow(10, -10);
    bool verbose_j = false, verbose_e = false;
    switch (init) {
        case 1:
            (void)smacofTorgerson(delta, xini, pn, pp, &itmax_j, &eps_j,
                                  &verbose_j);
            break;
        case 2:
            //(void)smacofElegant(delta, xini, pn, pp, &itmax_e, &eps_e,
            //&verbose_e);
            break;
        case 3:
            (void)smacofInitSDCmatrix(delta, weights, xini, pn, pp, &itmax_j,
                                      &eps_j, &verbose_j);
            break;
        case 4:
            (void)smacofInitRandom(xini, pn, pp);
            break;
    }
    (void)smacofCenter(xini, pn, pp);
    (void)smacofDistance(xini, dini, pn, pp);
    (void)smacofScale(delta, weights, dini, xini, pn, pp);
    return;
    // maybe diagonal improvement
}

void smacofInitRandom(double *xini, const int *pn, const int *pp) {
    int np = *pn * *pp;
    for (int i = 1; i <= np; i++) {
        xini[VINDEX(i)] = drand48();
    }
    return;
}

void smacofInitSDCmatrix(const double *delta, const double *weights,
                         double *xini, const int *pn, const int *pp,
                         const int *pitmax_j, const int *peps_j,
                         const bool *pverbose_j) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2;
    int width = 15, precision = 10;
    double *a = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *b = (double *)calloc((size_t)(m + n), (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= m; i++) {
        int iv = VINDEX(i);
        a[iv] = -weights[iv] * SQUARE(delta[iv]);
    }
    (void)smacofAddSDCDiagonal(a, b, pn);
    (void)smacofJacobi(b, evec, eval, pn, pitmax_j, peps_j, pverbose_j);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    return;
}

void smacofDiagonalImprove() {}
