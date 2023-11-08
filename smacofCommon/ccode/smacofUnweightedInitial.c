#include "../../smacofInclude/smacof.h"

void smacofUnweightedInitial(const double *delta, double *xini,
                             const int *pinit, const int *pn, const int *pp) {
    int init = *pinit, n = *pn, m = n * (n - 1) / 2;
    double *dini = (double *)calloc((size_t)m, (size_t)sizeof(double));
    switch (init) {
        case 1:
            (void)smacofUnweightedInitTorgerson(delta, xini, pn, pp);
            break;
        case 2:
            (void)smacofUnweightedInitMaximumSum(delta, xini, pn, pp);
            break;
        case 3:
            (void)smacofInitRandom(xini, pn, pp);
            break;
        case 4:;
            break;
    }
    (void)smacofCenter(xini, pn, pp);
    (void)smacofDistance(xini, dini, pn, pp);
    (void)smacofUnweightedScale(delta, dini, xini, pn, pp);
    free(dini);
    return;
}

void smacofUnweightedInitTorgerson(const double *delta, double *xold,
                                   const int *pn, const int *pp) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2, nn = n * n, itmax = 100,
        eps = 10;
    bool verbose = false;
    double *dimp = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *cross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        sum += delta[VINDEX(k)];
    }
    for (int k = 1; k <= nn; k++) {
        dimp[VINDEX(k)] = sum;
    }
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            int ij = SINDEX(i, j, n);
            dimp[ij] = delta[ij];
        }
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

void smacofUnweightedInitMaximumSum(const double *delta, double *xini,
                                    const int *pn, const int *pp) {
    int n = *pn, p = *pp, nn = n * (n + 1) / 2,
        itmax = 100, eps = 10;
    bool verbose = false;
    double *b = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int k = 1; k <= nn; k++) {
        b[VINDEX(k)] = 0.0;
    }
    // double index
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            double cell = SQUARE(delta[SINDEX(i, j, n)]);
            b[TINDEX(i, j, n)] -= cell;
            b[TINDEX(i, i, n)] += cell;
            b[TINDEX(j, j, n)] += cell;
        }
    }
    (void)smacofJacobi(b, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(b);
    free(evec);
    free(eval);
    return;
}
