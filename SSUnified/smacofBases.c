#include "smacof.h"

int smacofDoubleComparison(const void *a, const void *b) {
    return ((int)copysign(1.0, *(double *)a - *(double *)b));
}

bool smacofCheckIncreasing(double *innerknots, const double lowend,
                           const double highend, const int ninner) {
    if (ninner == 0) {
        return false;
    }
    if (lowend >= innerknots[VINDEX(1)]) {
        return true;
    }
    if (highend <= innerknots[VINDEX(ninner)]) {
        return true;
    }
    for (int i = 1; i < ninner; i++) {
        if (innerknots[i] <= innerknots[i - 1]) {
            return true;
        }
    }
    return false;
}

void smacofExtendPartition(const double *innerknots, const int *multiplicities,
                           const int *order, const int *ninner,
                           const double *lowend, const double *highend,
                           double *extended) {
    int k = 1;
    for (int i = 1; i <= *order; i++) {
        extended[VINDEX(k)] = *lowend;
        k++;
    }
    for (int j = 1; j <= *ninner; j++)
        for (int i = 1; i <= multiplicities[VINDEX(j)]; i++) {
            extended[VINDEX(k)] = innerknots[VINDEX(j)];
            k++;
        }
    for (int i = 1; i <= *order; i++) {
        extended[VINDEX(k)] = *highend;
        k++;
    }
}

void smacofBisect(const double *x, const double *knots, const int *lowindex,
                  const int *highindex, int *index) {
    int l = *lowindex, u = *highindex, mid = 0;
    while ((u - l) > 1) {
        mid = (int)floor((u + l) / 2);
        if (*x < knots[VINDEX(mid)])
            u = mid;
        else
            l = mid;
    }
    *index = l;
    return;
}

void smacofMakeInnerKnots(const int ninner, const bool percentiles, const int n,
                          const double lowend, const double highend,
                          double *dsort, double *innerknots) {
    if (ninner == 0) {
        return;
    }
    for (int i = 1; i <= ninner; i++) {
        innerknots[i - 1] = lowend + (i / (ninner + 1.0)) * (highend - lowend);
    }
    if (percentiles) {
        (void)qsort(dsort, (size_t)n, (size_t)sizeof(double),
                    smacofDoubleComparison);
        for (int i = 1; i <= ninner; i++) {
            double q = i / (ninner + 1.0);
            innerknots[i - 1] = smacofPercentile(n, q, dsort);
        }
    }
    return;
}

double smacofPercentile(const int n, const double q, double *dsort) {
    // R7 method of Hyndman and Fan
    if (q == 1.0) {
        return (dsort[n - 1]);
    }
    int n1 = (int)floor(q * (n - 1) + 1) - 1;
    int n2 = n1 + 1;
    double r = q * (n + 1.0) - floor(q * (n + 1.0));
    double x = dsort[n1] + r * (dsort[n2] - dsort[n1]);
    return x;
}

void smacofBsplines(const double *x, const double *knots, const int *order,
                    const int *nknots, int *index, double *q) {
    int lowindex = 1, highindex = *nknots, m = *order, j, jp1;
    double drr, dll, saved, term;
    double *dr = (double *)calloc((size_t)m, sizeof(double));
    double *dl = (double *)calloc((size_t)m, sizeof(double));
    (void)smacofBisect(x, knots, &lowindex, &highindex, index);
    int l = *index;
    for (j = 1; j <= m; j++) {
        q[VINDEX(j)] = 0.0;
    }
    if (*x == knots[VINDEX(*nknots)]) {
        q[VINDEX(m)] = 1.0;
        return;
    }
    q[VINDEX(1)] = 1.0;
    j = 1;
    if (j >= m) return;
    while (j < m) {
        dr[VINDEX(j)] = knots[VINDEX(l + j)] - *x;
        dl[VINDEX(j)] = *x - knots[VINDEX(l + 1 - j)];
        jp1 = j + 1;
        saved = 0.0;
        for (int r = 1; r <= j; r++) {
            drr = dr[VINDEX(r)];
            dll = dl[VINDEX(jp1 - r)];
            term = q[VINDEX(r)] / (drr + dll);
            q[VINDEX(r)] = saved + drr * term;
            saved = dll * term;
        }
        q[VINDEX(jp1)] = saved;
        j = jp1;
    }
    free(dr);
    free(dl);
    return;
}

void smacofBsplineBasis(const double *x, double *knots, const int *order,
                        const int *nknots, const int *nvalues, double *result) {
    int m = *order, l = 0;
    double *q = (double *)calloc((size_t)m + 1, sizeof(double));
    for (int i = 1; i <= *nvalues; i++) {
        (void)smacofBsplines(x + VINDEX(i), knots, order, nknots, &l, q);
        for (int j = 1; j <= m; j++) {
            // int r = MIN(l - m + j, *nknots - m);
            result[MINDEX(i, l - m + j, *nvalues)] = q[VINDEX(j)];
        }
    }
    free(q);
    return;
}

// Given a matrix X form the matrix of decreasing sums in the same space
// y_{ij} = sum_{k=j}^m x_{ik}

void smacofCumsumMatrix(const int n, const int m, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = m - 2; j >= 0; j--) {
            x[i][j] += x[i][j + 1];
        }
    }
    return;
}

void smacofBernsteinBasis(const int n, const int m, const double *y,
                          const bool ordinal, double **z) {
    double fac = 0.0, rat = 0.0;
    for (int i = 0; i < n; i++) {
        z[i][0] = pow(1.0 - y[i], (double)(m - 1));
        for (int j = 1; j < m; j++) {
            rat = ((double)(m - j)) / ((double)j);
            if (y[i] == 0.0) {
                z[i][j] = 0.0;
                continue;
            }
            if (y[i] == 1.0) {
                if (j == (m - 1)) {
                    z[i][j] = 1.0;
                } else {
                    z[i][j] = 0.0;
                }
                continue;
            }
            fac = y[i] / (1.0 - y[i]);
            z[i][j] = rat * fac * z[i][j - 1];
        }
    }
    if (ordinal) {
        (void)smacofCumsumMatrix(n, m, z);
    }
    return;
}
