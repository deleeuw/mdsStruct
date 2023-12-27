#include "smacof.h"

void smacofDoubleCenter(const int n, double **x, double **y) {
    double tsum = 0.0;
    double *rsum = smacofMakeVector(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            sum += SQUARE(x[ij][ji]);
        }
        rsum[i] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            y[i][j] = -0.5 * (SQUARE(x[i][j]) - rsum[i] - rsum[j] + tsum);
        }
    }
    free(rsum);
    return;
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
