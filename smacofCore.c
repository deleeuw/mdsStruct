#include "smacof.h"

void smacofDist(const double *x, double *d, const int *ii, const int *jj,
                const int *pm, const int *pn, const int *pp) {
    int n = *pn, p = *pp, m = *pm;
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        int i = ii[VINDEX(k)], j = jj[VINDEX(k)];
        sum = 0.0;
        for (int s = 1; s <= p; s++) {
            sum += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
        }
        d[SINDEX(i, j, n)] = sqrt(fabs(sum));
    }
    return;
}

void smacofGuttman(const double *delta, const double *w, const double *vinv,
                   const double *d, const double *x, double *z, const int *pn,
                   const int *pp) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2;
    double *b = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *y = (double *)calloc((size_t)(n * p), (size_t)sizeof(double));
    for (int i = 1; i <= m; i++) {
        b[VINDEX(i)] = w[VINDEX(i)] * delta[VINDEX(i)] / d[VINDEX(i)];
    }
    (void)smacofMultiplySDCMatrix(b, x, y, pn, pp);
    (void)smacofMultiplySDCMatrix(vinv, y, z, pn, pp);
    free(b);
    free(y);
    return;
}

void smacofLoss(const double *delta, const double *w, const double *d,
                const int *m, double *loss) {}
