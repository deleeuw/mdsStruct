#include "../../smacofInclude/smacof.h"

void smacofWeightedInitial(const double *delta, const double *weights,
                           const int *irow, const int *icol, double *xini,
                           const int *pinit, const int *pn, const int *pp,
                           const int *pm) {
    int init = *pinit, np = *pn * *pp;
    int itmax_j = 100, itmax_d = 100, eps_j = 10, eps_d = 10;
    int width = 15, precision = 10;
    bool verbose_j = false, verbose_d = false;
    double *dini = (double *)calloc((size_t)np, (size_t)sizeof(double));
    switch (init) {
        case 1:
            (void)smacofWeightedInitTorgerson(delta, weights, irow, icol, xini,
                                              pn, pp, pm);
            break;
        case 2:
            (void)smacofWeightedInitMaximumSum(delta, weights, irow, icol, xini,
                                               pn, pp, pm);
            break;
        case 3:
            (void)smacofInitRandom(xini, pn, pp);
            break;
        case 4:;
            break;
    }
    (void)smacofCenter(xini, pn, pp);
    (void)smacofDistance(xini, dini, pn, pp);
    (void)smacofWeightedScale(delta, weights, dini, xini, pn, pp, pm);
    free(dini);
    return;
}

void smacofWeightedInitTorgerson(const double *delta, const double *weights,
                                 const int *irow, const int *icol, double *xold,
                                 const int *pn, const int *pp, const int *pm) {
    int n = *pn, p = *pp, m = *pm, nn = n * (n - 1) / 2, itmax = 100, eps = 10;
    int width = 15, precision = 10, nrow = 1;
    bool verbose = false;
    double *dimp = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *cross = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        sum += weights[VINDEX(k)] * delta[VINDEX(k)];
    }
    for (int k = 1; k <= nn; k++) {
        dimp[VINDEX(k)] = sum;
    }
    for (int k = 1; k <= m; k++) {
        int ki = irow[VINDEX(k)];
        int kj = icol[VINDEX(k)];
        dimp[SINDEX(ki, kj, n)] = delta[VINDEX(k)];
    }
    (void)smacofDoubleCenter(dimp, cross, pn);
    (void)smacofJacobi(cross, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(fabs(eval[VINDEX(s)]));
        }
    }
    free(cross);
    free(evec);
    free(eval);
    free(dimp);
    return;
}

void smacofWeightedInitMaximumSum(const double *delta, const double *weights,
                                  const int *irow, const int *icol,
                                  double *xini, const int *pn, const int *pp,
                                  const int *pm) {
    int n = *pn, p = *pp, m = *pm, nn = n * (n - 1) / 2, itmax = 100, eps = 10;
    bool verbose = false;
    double *a = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *b = (double *)calloc((size_t)(m + n), (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int k = 1; k <= nn; k++) {
        b[VINDEX(k)] = 0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double cell = weights[kv] * SQUARE(delta[kv]);
        b[TINDEX(ik, jk, n)] -= cell;
        b[TINDEX(ik, ik, n)] += cell;
        b[TINDEX(jk, jk, n)] += cell;
    }
    (void)smacofJacobi(b, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(a);
    free(b);
    free(evec);
    free(eval);
    return;
}
