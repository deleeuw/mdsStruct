#include "smacof.h"

void smacofDist(const double *x, double *d, const int *ii, const int *jj,
                const int *pm, const int *pn, const int *pp) {
    int n = *pn, p = *pp, m = *pm, width = 15, precision = 10;
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            double sum = 0.0;
            for (int s = 1; s <= p; s++) {
                sum += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            d[SINDEX(i, j, n)] = sqrt(fabs(sum));
        }
    }
    return;
}

void smacofGuttman(const double *delta, const double *weights,
                   const double *vinv, const double *dold, double *xold,
                   double *xnew, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p, m = n * (n - 1) / 2;
    double *bmat = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *ymat = (double *)calloc((size_t)np, (size_t)sizeof(double));
    for (int i = 1; i <= m; i++) {
        bmat[VINDEX(i)] =
            weights[VINDEX(i)] * delta[VINDEX(i)] / dold[VINDEX(i)];
    }
    (void)smacofMultiplySDCMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySDCMatrix(vinv, ymat, xnew, pn, pp);
    free(bmat);
    free(ymat);
    return;
}

void smacofLoss(const double *delta, const double *weights, const double *dist,
                const int *pm, double *loss) {
    int m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        sum += weights[VINDEX(k)] * SQUARE(delta[VINDEX(k)] - dist[VINDEX(k)]);
    }
    *loss = sum / 2.0;
}
