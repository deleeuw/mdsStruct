#include "smacofCommon.h"

void smacofInitRandom(const unsigned n, const unsigned p, double *xini) {
    unsigned np = n * p;
    for (unsigned i = 1; i <= np; i++) {
        xini[VINDEX(i)] = drand48();
    }
    return;
}

void smacofDoubleCenter(const unsigned n, const double *delta, double *cross) {
    double tsum = 0.0, cell = 0.0;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    assert(!(rsum == NULL));
    for (unsigned i = 1; i <= n; i++) {
        double sum = 0.0;
        for (unsigned j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE(delta[PINDEX(i, j, n)]);
        }
        rsum[VINDEX(i)] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (unsigned j = 1; j <= n; j++) {
        for (unsigned i = j; i <= n; i++) {
            unsigned ij = TINDEX(i, j, n);
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

double smacofMaxConfigurationDifference(const unsigned n, const unsigned p,
                                        const double *x, const double *y) {
    double maxdiff = 0.0;
    for (unsigned s = 1; s <= p; s++) {
        double posdiff = 0.0, negdiff = 0.0;
        for (unsigned i = 1; i <= n; i++) {
            posdiff =
                MAX(posdiff, fabs(x[MINDEX(i, s, n)] - y[MINDEX(i, s, n)]));
            negdiff =
                MAX(negdiff, fabs(x[MINDEX(i, s, n)] + y[MINDEX(i, s, n)]));
        }
        maxdiff = MAX(maxdiff, MIN(posdiff, negdiff));
    }
    return maxdiff;
}

double smacofMaxDistanceDifference(const unsigned m, const double *dold,
                                   const double *dnew) {
    double dchange = 0.0;
    for (unsigned k = 1; k <= m; k++) {
        dchange = MAX(dchange, fabs(dold[VINDEX(k)] - dnew[VINDEX(k)]));
    }
    return dchange;
}

double smacofRMSDifference(const unsigned n, const unsigned p, const double *x,
                           double *y) {
    unsigned np = n * p;
    double sum = 0.0;
    for (unsigned s = 1; s <= p; s++) {
        double sump = 0.0, sumn = 0.0;
        for (unsigned i = 1; i <= np; i++) {
            unsigned k = VINDEX(i);
            sump += SQUARE(x[k] - y[k]);
            sumn += SQUARE(x[k] + y[k]);
        }
        if (sumn < sump) {
            for (unsigned i = 1; i <= np; i++) {
                unsigned k = VINDEX(i);
                y[k] = -y[k];
            }
            sum += sumn;
        } else {
            sum += sump;
        }
    }
    return sqrt(fabs(sum / ((double)np)));
}

void smacofMakeIIandJJ(const unsigned n, unsigned *ii, unsigned *jj) {
    for (unsigned j = 1; j <= (n - 1); j++) {
        for (unsigned i = (j + 1); i <= n; i++) {
            unsigned ij = SINDEX(i, j, n);
            ii[ij] = i;
            jj[ij] = j;
        }
    }
    return;
}
