#include "smacof.h"

void smacofRelax(const double *xold, double *xnew, const double *pechange,
                 const double *ppchange, const int *pnp, const int *pitel,
                 const bool *prelax, double *prate) {
    int np = *pnp, itel = *pitel;
    double echange = *pechange, pchange = *ppchange, rate = *prate;
    bool relax = *prelax;
    if (itel == 1) {
        rate = NaN;
    } else {
        rate = echange / pchange;
        if (relax) {
            double eopt = MAX(0, MIN(1, rate / (2 - rate)));
            for (int i = 1; i <= np; i++) {
                int iv = VINDEX(i);
                xnew[iv] = (1 + eopt) * xnew[iv] - eopt * xold[iv];
            }
        }
    }
    *prate = rate;
}

void smacofDiagonalAdjust(const double *delta, const double *weights, double *x,
                          const int *pn, const int *pp, const int *pitmax,
                          const int *peps, const bool *pverbose) {
    int n = *pn, p = *pp, np = n * p, m = n * (n - 1) / 2, itmax = *pitmax,
        itel = 1, width = 15, precision = 10;
    double sum = 0.0, ssum = 0.0, diff = 0.0, lwk = 0.0,
           eps = pow(10.0, -(double)*peps), stress = 0.0;
    bool verbose = *pverbose;
    double *lbd = (double *)calloc((size_t)p, (size_t)sizeof(double));
    double *dia = (double *)calloc((size_t)p, (size_t)sizeof(double));
    double *xwk = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *dwk = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *bwk = (double *)calloc((size_t)m, (size_t)sizeof(double));
    for (int s = 1; s <= p; s++) {
        lbd[s] = 1.0;
        ssum = 0.0;
        for (int i = 1; i <= n; i++) {
            sum = 0.0;
            for (int j = 1; j <= n; j++) {
                if (i == j) {
                    continue;
                }
                sum += weights[PINDEX(i, j, n)] *
                       (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            ssum += sum * x[MINDEX(i, s, n)];
        }
        dia[s] = ssum;
    }
    while (true) {
        for (int i = 1; i <= n; i++) {
            for (int s = 1; s <= p; s++) {
                xwk[MINDEX(i, s, n)] =
                    x[MINDEX(i, s, n)] * lbd[VINDEX(s)];
            }
        }
        (void)smacofDistance(xwk, dwk, pn, pp);
        (void)smacofMakeBMatrix(delta, weights, dwk, bwk, &m);
        (void)smacofPrintSDCMatrix(bwk, pn, &width, &precision);
        (void)smacofStress(delta, weights, dwk, &m, &stress);
        printf("stress %15.10f\n", stress);
        diff = 0.0;
        for (int s = 1; s <= p; s++) {
            sum = 0.0;
            for (int t = 1; t <= p; t++) {
                if (s == t) {
                    continue;
                }
                sum += bwk[PINDEX(s, t, p)] * (lbd[VINDEX(s)] - lbd[VINDEX(t)]);
            }
            lwk = fabs(sum / dia[s]);
            diff = MAX(diff, fabs(lbd[s] - lwk));
            lbd[s] = lwk;
        }
        if (verbose) {
            printf("%3d", itel);
            for (int s = 1; s <= p; s++) {
                printf("%15.10f", lbd[s]);
            }
            printf("\n");
        }
        if ((itel == itmax) || (diff < eps)) {
            break;
        }
        itel++;
    }
    (void)memcpy(x, xwk, (size_t)np * sizeof(double));
    free(lbd);
    free(dia);
    free(xwk);
    free(dwk);
    free(bwk);
    return;
}

/*
 * int main() {
    int n = 4, p = 2, itmax = 100, eps = 6, width = 15, precision = 10,
        m = n * (n - 1) / 2;
    bool verbose = true;
    double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double weights[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double x[8] = {1.0, 1.0, -1.0, -1.0, -4.0, 4.0, -4.0, 4.0};
    (void)smacofNormWeights(weights, &m);
    (void)smacofNormDelta(delta, weights, &m);
    (void)smacofDiagonalAdjust(delta, weights, x, &n, &p, &itmax, &eps,
                               &verbose);
    (void)smacofPrintAnyMatrix(x, &n, &p, &width, &precision);
}
*/
