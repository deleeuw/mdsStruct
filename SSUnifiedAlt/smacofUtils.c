#include "smacof.h"

void smacofDoubleCenter(const int n, double **x, double **y) {
    double tsum = 0.0, cell = 0.0;
    double *rsum = smacofMakeVector(n);
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
    double stot = 0.0;
    for (int s = 0; s < p; s++) {
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < n; i++) {
            s1 += SQUARE(x[i][s] - y[i][s]);
            s2 += SQUARE(x[i][s] + y[i][s]);
        }
        stot += MIN(s1, s2);
    }
    return sqrt(stot / ((double)np));
}

double smacofEtaSquareDifference(const int n, const int p, const bool weights,
                                 double **vmat, double **x, double **y) {
    double stot = 0.0, fac = 0.0;
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double dfi = x[i][s] - y[i][s];
            for (int j = 0; j <= i; j++) {
                if (weights) {
                    fac = vmat[i][j];
                } else {
                    if (i == j) {
                        fac = (double)n - 1.0;
                    } else {
                        fac = -2.0;
                    }
                }
                double dfj = x[j][s] - y[j][s];
                stot += fac * dfi * dfj;
            }
        }
    }
    return stot;
}

void smacofZeroAnyMatrix(const int n, const int p, double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            x[i][s] = 0.0;
        }
    }
    return;
}
