#include "smacofUnweighted.h"

void smacofUnweightedInitial(const unsigned n, const unsigned p,
                             const double *delta, double *xini,
                             const unsigned init) {
    unsigned m = n * (n - 1) / 2;
    double *dini = (double *)calloc((size_t)m, (size_t)sizeof(double));
    assert(dini != NULL);
    switch (init) {
        case 1:
            (void)smacofUnweightedInitTorgerson(n, p, delta, xini);
            break;
        case 2:
            (void)smacofUnweightedInitMaximumSum(n, p, delta, xini);
            break;
        case 3:
            (void)smacofInitRandom(n, p, xini);
            break;
        case 4:;
            break;
    }
    (void)smacofCenter(n, p, xini);
    (void)smacofDistance(n, p, xini, dini);
    (void)smacofUnweightedScale(n, p, delta, dini, xini);
    free(dini);
    return;
}

void smacofUnweightedInitTorgerson(const unsigned n, const unsigned p,
                                   const double *delta, double *xold) {
    unsigned m = n * (n - 1) / 2, nn = n * n, itmax = 100, eps = 10;
    bool verbose = false;
    double *dimp = (double *)calloc((size_t)m, (size_t)sizeof(double));
    assert(dimp != NULL);
    double *cross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    assert(cross != NULL);
    double *evec = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    assert(evec != NULL);
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    assert(eval != NULL);
    double sum = 0.0;
    for (unsigned k = 1; k <= m; k++) {
        sum += delta[VINDEX(k)];
    }
    for (unsigned k = 1; k <= nn; k++) {
        dimp[VINDEX(k)] = sum;
    }
    for (unsigned j = 1; j <= (n - 1); j++) {
        for (unsigned i = (j + 1); i <= n; i++) {
            unsigned ij = SINDEX(i, j, n);
            dimp[ij] = delta[ij];
        }
    }
    (void)smacofDoubleCenter(n, dimp, cross);
    (void)smacofJacobi(n, p, cross, evec, eval, itmax, eps, verbose);
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned s = 1; s <= p; s++) {
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

void smacofUnweightedInitMaximumSum(const unsigned n, const unsigned p,
                                    const double *delta, double *xini) {
    unsigned nn = n * (n + 1) / 2, itmax = 100, eps = 10;
    bool verbose = false;
    double *b = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    assert(b != NULL);
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    assert(evec != NULL);
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    assert(eval != NULL);
    for (unsigned k = 1; k <= nn; k++) {
        b[VINDEX(k)] = 0.0;
    }
    // double index
    for (unsigned j = 1; j <= (n - 1); j++) {
        for (unsigned i = (j + 1); i <= n; i++) {
            double cell = SQUARE(delta[SINDEX(i, j, n)]);
            b[TINDEX(i, j, n)] -= cell;
            b[TINDEX(i, i, n)] += cell;
            b[TINDEX(j, j, n)] += cell;
        }
    }
    (void)smacofJacobi(n, p, b, evec, eval, itmax, eps, verbose);
    for (unsigned i = 1; i <= n; i++) {
        for (unsigned s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(b);
    free(evec);
    free(eval);
    return;
}
