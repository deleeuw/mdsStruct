#include "../Include/smacof.h"

void smacofInitRandom(const size_t n, const size_t p, double (*xini)[n]) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < p; j++) {
            xini[i][j] = drand48();
        }
    }
    return;
}

void smacofDoubleCenter(const size_t n, const double delta[n][n],
                        double cross[n][n]) {
    double tsum = 0.0, cell = 0.0;
    double *rsum = malloc((sizeof *rsum) * n);
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE(delta[i][j]);
        }
        rsum[i] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            if (i == j) {
                cell = 0.0;
            } else {
                cell = SQUARE(delta[i][j]);
            }
            cross[i][j] = -0.5 * (cell - rsum[i] - rsum[j] + tsum);
        }
    }
    free(rsum);
    return;
}

void smacofMaxConfigurationDifference(const size_t n, const size_t p,
                                      const double x[n][p],
                                      const double y[n][p], double *maxdiff) {
    *maxdiff = 0.0;
    for (size_t s = 0; s < p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (size_t i = 0; i < n; i++) {
            posdiff = MAX(posdiff, fabs(x[i][s] - y[i][s]));
            negdiff = MAX(negdiff, fabs(x[i][s] + y[i][s]));
        }
        *maxdiff = MAX(*maxdiff, MIN(posdiff, negdiff));
    }
    return;
}

void smacofMaxDistanceDifference(const size_t n, const double dold[n][n],
                                 const double dnew[n][n], double *pdchange) {
    double dchange = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            dchange = MAX(dchange, fabs(dold[i][j] - dnew[i][j]));
        }
    }
    *pdchange = dchange;
    return;
}

void smacofRMSDifference(const size_t n, const size_t p, const double x[n][p],
                         double y[n][p], double *diff) {
    size_t np = n * p;
    double sum = 0.0;
    for (size_t s = 0; s < p; s++) {
        double sump = 0.0, sumn = 0.0;
        for (size_t i = 0; i <= n; i++) {
            sump += SQUARE(x[i][s] - y[i][s]);
            sumn += SQUARE(x[i][s] + y[i][s]);
        }
        if (sumn < sump) {
            for (size_t i = 0; i < n; i++) {
                y[i][s] = -y[i][s];
            }
            sum += sumn;
        } else {
            sum += sump;
        }
    }
    *diff = sqrt(fabs(sum / ((double)np)));
    return;
}
