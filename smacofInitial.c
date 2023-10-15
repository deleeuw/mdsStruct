#include "smacof.h"

void smacofInitial(const double *delta, const int *pn, const int *pp,
                   double *x) {}

void smacofDoubleCenter(const double *delta, double *cross, const int *pn) {
    int n = *pn, ij = 0;
    double tsum = 0.0, sum = 0.0, cell = 0.0;
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
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            ij = SINDEX(i, j, n);
            cell = SQUARE(delta[ij]);
            cross[ij] =
                -0.5 * (cell - rsum[VINDEX(i)] - rsum[VINDEX(j)] + tsum);
        }
    }
    free(rsum);
    return;
}

/*
int main() {
    double delta[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double cross[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vinv[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double w[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    int n = 5;
    (void)smacofDoubleCenter(delta, cross, &n);
    (void)smacofVinverse(w, vinv, &n);
    (void)smacofPrintMM(cross, &n);
    (void)smacofPrintMM(vinv, &n);
    return (EXIT_SUCCESS);
}
 */