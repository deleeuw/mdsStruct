#include "../Include/smacof.h"

void smacofUnweightedInitial(const size_t n, const size_t p, const size_t init,
                             const double (*delta)[n], double (*xini)[n]) {
    size_t nn = SQUARE(n);
    double(*dini)[n] = malloc(sizeof(*dini) * n);
    assert(!(dini == NULL));
    switch (init) {
        case 1:
            //(void)smacofUnweightedInitTorgerson(n, p, delta, xini);
            break;
        case 2:
            //(void)smacofUnweightedInitMaximumSum(n, p, delta, xini);
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

/*

void smacofUnweightedInitTorgerson(const int n, const int p,
                                   const double delta[n][n],
                                   double xold[n][p]) {
    size_t m = n * (n - 1) / 2, nn = SQUARE(n), size_t itmax = 100,
             eps = 10;
    bool verbose = false;
    double(*dimp)[n][n] = calloc((size_t)nn, (size_t)sizeof(*dinp));
    double(*cross)[n][n] = calloc((size_t)nn, (size_t)sizeof(*cbold));
    double(*evec)[n][n] = calloc((size_t)nn, (size_t)sizeof(*cdini));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double sum = 0.0;
    for (size_t k = 1; k <= m; k++) {
        sum += delta[VINDEX(k)];
    }
    for (size_t k = 1; k <= nn; k++) {
        dimp[VINDEX(k)] = sum;
    }
    for (size_t j = 1; j <= (n - 1); j++) {
        for (size_t i = (j + 1); i <= n; i++) {
            size_t ij = SINDEX(i, j, n);
            dimp[ij] = delta[ij];
        }
    }
    (void)smacofDoubleCenter(dimp, cross, pn);
    (void)smacofJacobi(n, p, cross, *evec, *eval, itmax, eps, verbose);
    for (size_t i = 0; i < n; i++) {
        for (size_t s = 0; s < p; s++) {
            xold[i][s] =
                evec[i][s] * sqrt(fabs(eval[s]));
        }
    }
    free(cross);
    free(evec);
    free(eval);
    free(dimp);
    return;
}

void smacofUnweightedInitMaximumSum(const size_t n, const size_t p,
                                    const double delta[n][n],
                                    double xold[n][p]) {
    size_t nn = n * (n + 1) / 2, itmax = 100, eps = 10;
    bool verbose = false;
    double *b = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (size_t k = 1; k <= nn; k++) {
        b[VINDEX(k)] = 0.0;
    }
    // double index
    for (size_t j = 1; j <= (n - 1); j++) {
        for (size_t i = (j + 1); i <= n; i++) {
            double cell = SQUARE(delta[SINDEX(i, j, n)]);
            b[TINDEX(i, j, n)] -= cell;
            b[TINDEX(i, i, n)] += cell;
            b[TINDEX(j, j, n)] += cell;
        }
    }
    (void)smacofJacobi(b, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (size_t i = 1; i <= n; i++) {
        for (size_t s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(b);
    free(evec);
    free(eval);
    return;
}
*/
