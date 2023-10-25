#include "smacof.h"

void smacofInitial(const double *delta, const double *weights, double *xini,
                   const int *init, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    int itmax_j = 100, itmax_e = 100;
    int eps_j = 10, eps_e = 10;
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
            (void)smacofInitSDCmatrix(delta, weights, pn, &itmax_j, &eps_j,
                                      &verbose_j);
            break;
        case 4:
            (void)smacofInitRandom(xini, &np) break;
    }
    // scale in here
    // maybe diagonal improvement
}

void smacofInitRandom(double *xini, int *np) {
    int np = *np;
    for (int i = 1; i <= np; i++) {
        xini[VINDEX(i)] = drand48();
    }
    return;
}

void smacofInitSDCmatrix(const double *delta, const double *weights,
                         double *xini, const int *pm) {
    int m = *pm;
    double *cross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    for (int i = 1; i <= m; i++) {
        int iv = VINDEX(i);
        cross[iv] = -weights[iv] * SQUARE(delta[iv]);
    }
    // give it a diagonal
    (void)smacofJacobi(cross, evec, eval, pn, pitmax, peps, pverbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    return;
}