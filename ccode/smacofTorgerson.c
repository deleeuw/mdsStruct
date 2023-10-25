#include "smacof.h"

void smacofTorgerson(const double *delta, double *xold, const int *pn,
                     const int *pp, const int *pitmax, const int *peps,
                     const bool *pverbose) {
    int n = *pn, p = *pp, m = n * (n + 1) / 2;
    double *cross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    (void)smacofDoubleCenter(delta, cross, pn);
    (void)smacofJacobi(cross, evec, eval, pn, pitmax, peps, pverbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(cross);
    free(evec);
    free(eval);
    return;
}

void smacofDoubleCenter(const double *delta, double *cross, const int *pn) {
    int n = *pn, ij = 0;
    double tsum = 0.0, cell = 0.0;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE(delta[PINDEX(i, j, n)]);
        }
        rsum[VINDEX(i)] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            ij = TINDEX(i, j, n);
            if (i == j) {
                cell = 0.0;
            } else {
                cell = SQUARE(delta[SINDEX(i, j, n)]);
            }
            cross[ij] =
                -0.5 * (cell - rsum[VINDEX(i)] - rsum[VINDEX(j)] + tsum);
        }
    }
    free(rsum);
    return;
}

/*
int main() {
    double delta[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    double xold[10]= {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    int n = 5, p = 2, width = 15, precision = 10;
    (void)smacofTorgerson(delta, xold, &n, &p);
    (void)smacofPrintAnyMatrix(xold, &n, &p, &width, &precision);
    return (EXIT_SUCCESS);
}
*/
