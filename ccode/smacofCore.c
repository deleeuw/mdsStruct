#include "smacof.h"

void smacofDist(const double *x, double *d, const int *ii, const int *jj,
                const int *pm, const int *pn, const int *pp) {
    int n = *pn, p = *pp;
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

void smacofMakeBmat(const double *delta, const double *weights,
                    const double *dold, double *bmat, const int *pm) {
    int m = *pm;
    for (int i = 1; i <= m; i++) {
        int iv = VINDEX(i);
        double dfix = dold[iv];
        if (dfix < 1e-15) {
            bmat[iv] = 0.0;
        } else {
            bmat[iv] = -weights[iv] * delta[iv] / dfix;
        }
    }
    return;
}

// void smacofInitial() {}

void smacofGuttman(const double *vinv, const double *bmat, const double *xold,
                   double *xnew, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    int width = 15, precision = 10;
    double *ymat = (double *)calloc((size_t)np, (size_t)sizeof(double));
    if (DEBUG) {
        printf("in guttman before matmult\n\n");
        printf("xold in %p\n\n", xold);
        (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
        printf("xnew in %p\n\n", xnew);
        (void)smacofPrintAnyMatrix(xnew, pn, pp, &width, &precision);
    }
    (void)smacofMultiplySDCMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySDCMatrix(vinv, ymat, xnew, pn, pp);
    if (DEBUG) {
        printf("in guttman after matmult\n\n");
        printf("xold in %p\n\n", xold);
        (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
        printf("xnew in %p\n\n", xnew);
        (void)smacofPrintAnyMatrix(xnew, pn, pp, &width, &precision);
    }
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
