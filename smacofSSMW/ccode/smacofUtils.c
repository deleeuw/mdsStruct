#include "smacof.h"

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

void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      const int *pn, const int *pp,
                                      double *maxdiff) {
    *maxdiff = 0.0;
    int n = *pn, p = *pp;
    for (int s = 1; s <= p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (int i = 1; i <= n; i++) {
            posdiff =
                MAX(posdiff, fabs(x[MINDEX(i, s, n)] - y[MINDEX(i, s, n)]));
            negdiff =
                MAX(negdiff, fabs(x[MINDEX(i, s, n)] + y[MINDEX(i, s, n)]));
        }
        *maxdiff = MAX(*maxdiff, MIN(posdiff, negdiff));
    }
    return;
}

void smacofMaxDistanceDifference(const double *dold, const double *dnew,
                                 const int *pm, double *pdchange) {
    double dchange = 0.0;
    int m = *pm;
    for (int k = 1; k <= m; k++) {
        dchange = MAX(dchange, fabs(dold[VINDEX(k)] - dnew[VINDEX(k)]));
    }
    *pdchange = dchange;
    return;
}

void smacofRMSDifference(const double *x, double *y, const int *pn,
                         const int *pp, double *diff) {
    int n = *pn, p = *pp, np = n * p;
    double sum = 0.0;
    for (int s = 1; s <= p; s++) {
        double sump = 0.0, sumn = 0.0;
        for (int i = 1; i <= np; i++) {
            int k = VINDEX(i);
            sump += SQUARE(x[k] - y[k]);
            sumn += SQUARE(x[k] + y[k]);
        }
        if (sumn < sump) {
            for (int i = 1; i <= np; i++) {
                int k = VINDEX(i);
                y[k] = -y[k];
            }
            sum += sumn;
        } else {
            sum += sump;
        }
        *diff = sqrt(fabs(sum / ((double)np)));
        return;
    }
}

void smacofMakeIIandJJ(const int *pn, int *ii, int *jj) {
    int n = *pn;
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            int ij = SINDEX(i, j, n);
            ii[ij] = i;
            jj[ij] = j;
        }
    }
    return;
}

/*
 int main() {
    int n = 4, p = 2, width = 10, precision = 6;
    double v[6] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0};
    double x[8] = {1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0};
    double y[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    (void)smacofPrintSDCLMatrix(v, &n, &width, &precision);
    (void)smacofPrintAnyMatrix(x, &n, &p, &width, &precision);
    (void)smacofMultiplySDCLMatrix(v, x, y, &n, &p);
    (void)smacofPrintAnyMatrix(y, &n, &p, &width, &precision);
    double sum = 0.0;
    for (int i = 1; i <= 4; i++) {
      for (int s = 1; s <= 2; s++) {
        sum += x[MINDEX(i, s, 4)] * y[MINDEX(i, s, 4)];
      }
    }
    printf("%10.6f\n", sum);
}
*/