#include "smacof.h"

void smacofDoubleCenter(const int n, double **x, double **y) {
    double tsum = 0.0, cell = 0.0;
    double *rsum = smacofMakeAnyVector(n);
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
            y[i][j] = -0.5 * (cell - rsum[i] - rsum[j] + tsum);
        }
    }
    free(rsum);
    return;
}

double smacofMaxConfigurationDifference(const int n, const int p, double **x,
                                        double **y) {
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

double smacofMaxDistanceDifference(const int n, double **x, double **y) {
    double dchange = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dchange = MAX(dchange, fabs(x[i][j] - y[i][j]));
        }
    }
    return dchange;
}

double smacofRMSDifference(const int n, const int p, double **x, double **y) {
    int np = n * p;
    double sumxx = 0.0, sumyy = 0.0, sumxy = 0.0;
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            sumxx += SQUARE(x[i][s]);
            sumyy += SQUARE(y[i][s]);
            sumxy += x[i][s] * y[i][s];
        }
    }
    double rms = sqrt(fabs((sumxx + sumyy - 2.0 * fabs(sumxy))) / ((double)np));
    return rms;
}

void smacofZeroAnyMatrix(const int n, const int p, double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            x[i][s] = 0.0;
        }
    }
    return;
}
