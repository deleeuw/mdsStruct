#include "../Include/smacof.h"

void smacofInitRandom(const unsigned n, const unsigned p, double (*xini)[n][p]) {
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            (*xini)[i][j] = drand48();
        }
    }
    return;
}

void smacofDoubleCenter(const unsigned n, const double (*delta)[n][n],
                        double (*cross)[n][n]) {
    double tsum = 0.0, cell = 0.0;
    double (*rsum)[n] = malloc(sizeof *rsum);
    for (unsigned i = 0; i < n; i++) {
        double sum = 0.0;
        for (unsigned j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE((*delta)[i][j]);
        }
        (*rsum)[i] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (unsigned j = 0; j < n; j++) {
        for (unsigned i = 0; i < n; i++) {
            if (i == j) {
                cell = 0.0;
            } else {
                cell = SQUARE((*delta)[i][j]);
            }
            (*cross)[i][j] = -0.5 * (cell - (*rsum)[i] - (*rsum)[j] + tsum);
        }
    }
    free(rsum);
    return;
}

double smacofMaxConfigurationDifference(const unsigned n, const unsigned p,
                                      const double (*x)[n][p],
                                      const double (*y)[n][p]) {
    double maxdiff = 0.0;
    for (unsigned s = 0; s < p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (unsigned i = 0; i < n; i++) {
            posdiff = MAX(posdiff, fabs((*x)[i][s] - (*y)[i][s]));
            negdiff = MAX(negdiff, fabs((*x)[i][s] + (*y)[i][s]));
        }
        maxdiff = MAX(maxdiff, MIN(posdiff, negdiff));
    }
    return maxdiff;
}

double smacofMaxDistanceDifference(const unsigned n, const double (*dold)[n][n],
                                 const double (*dnew)[n][n]) {
    double dchange = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            dchange = MAX(dchange, fabs((*dold)[i][j] - (*dnew)[i][j]));
        }
    }
    return dchange;
}

double smacofRMSDifference(const unsigned n, const unsigned p, const double (*x)[n][p],
                         double (*y)[n][p]) {
    unsigned np = n * p;
    double sum = 0.0;
    for (unsigned s = 0; s < p; s++) {
        double sump = 0.0, sumn = 0.0;
        for (unsigned i = 0; i <= n; i++) {
            sump += SQUARE((*x)[i][s] - (*y)[i][s]);
            sumn += SQUARE((*x)[i][s] + (*y)[i][s]);
        }
        if (sumn < sump) {
            for (unsigned i = 0; i < n; i++) {
              (*y)[i][s] = -(*y)[i][s];
            }
            sum += sumn;
        } else {
            sum += sump;
        }
    }
    double diff = sqrt(fabs(sum / ((double)np)));
    return diff;
}
