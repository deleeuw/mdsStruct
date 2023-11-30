#include "smacofCommon.h"

void smacofInitRandom(const int n, const int p,
                      double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            x[i][j] = drand48();
        }
    }
    return;
}

void smacofDoubleCenter(const int n, const double **x, double **y) {
    double tsum = 0.0, cell = 0.0;
    double *rsum = (double *)calloc((size_t) n, sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE(x[i][j]);
        }
        rsum[i] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (i == j) {
                cell = 0.0;
            } else {
                cell = SQUARE(x[i][j]);
            }
            cross[i][j] = -0.5 * (cell - (*rsum)[i] - (*rsum)[j] + tsum);
        }
    }
    free(rsum);
    return;
}

double smacofMaxConfigurationDifference(const int n, const int p,
                                        const double *x, const double *y) {
    double maxdiff = 0.0;
    for (int s = 0; s < p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (int i = 0; i < n; i++) {
            posdiff = MAX(posdiff, fabs(x[i][s] - y[i][s]));
            negdiff = MAX(negdiff, fabs(x[i][s] + y[i][s]));
        }
        maxdiff = MAX(maxdiff, MIN(posdiff, negdiff));
    }
    return maxdiff;
}

double smacofMaxDistanceDifference(const int n, const double **x,
                                   const double **y) {
    double dchange = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dchange = MAX(dchange, fabs(dold[i][j] - dnew[i][j]));
        }
    }
    return dchange;
}

double smacofRMSDifference(const int n, const int p,
                           const double **x, const double **y) {
    int np = n * p;
    double sum = 0.0;
    for (int s = 0; s < p; s++) {
        double sumxx = 0.0, sumyy = 0.0, sumxy = 0;
        for (int i = 0; i <= n; i++) {
            sumxx += SQUARE(x[i][s]);
            sumyy += SQUARE(x[i][s]);
            sumxy += x[i][s] * y[i][s];
        }
        sum += sumxx + sumyy - 2 * sumxy;
    }
    return sqrt(fabs(sum / ((double)np)));
}

void smacofZeroAnyMatrix(const int n, const int p,
                         double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            x[i][j] = 0.0;
        }
    }
    return;
}
