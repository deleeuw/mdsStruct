#include "smacof.h"

bool smacofCheckIncreasing(const int ninner, const double *innerknots,
                           const double lowend, const double highend) {
    if (lowend >= innerknots[0]) {
        return true;
    }
    if (highend <= innerknots[ninner - 1]) {
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

void smacofBsplineBasis(const double *x, const double *knots, const int *order,
                        const int *nknots, const int *nvalues, double *result) {
    int m = *order, l = 0;
    double *q = (double *)calloc((size_t)m + 1, sizeof(double));
    for (int i = 1; i <= *nvalues; i++) {
        (void)smacofBsplines(x + VINDEX(i), knots, order, nknots, &l, q);
        for (int j = 1; j <= m; j++) {
            int r = MIN(l - m + j, *nknots - m);
            result[MINDEX(i, r, *nvalues)] = q[VINDEX(j)];
        }
    }
    free(q);
    return;
}

// Given a matrix X form the matrix of decreasing sums in the same space
// y_{ij} = sum_{k=j}^m x_{ik}

void smacofCumsumMatrix(const int n, const int m, double **x) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = j; k < m; k++) {
                sum += x[i][k];
            }
            x[i][j] = sum;
        }
    }
    return;
}

void smacofBernsteinBase(const int n, const int m, const double *y,
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
    (void)smacofPrintAnyMatrix(stdout, n, m, 15, 10, z);
    return;
}

/*
int n = 20;
int m = 4;
double x[20] = {0.0000000000,  3.0000000000,  2.0000000000,  3.0000000000,
3.0000000000,  0.0000000000,  1.0000000000,  2.0000000000,
2.0000000000,  1.0000000000,  0.0000000000,  3.0000000000,
3.0000000000,  2.0000000000,  3.0000000000,  0.0000000000};
int precision = 6;
bool ordinal = false;

int main() {
    double **z = smacofMakeAnyMatrix(n, m);
    assert(!(z == NULL));
    double *y = malloc(n * sizeof(double));
    assert(!(y == NULL));
    double max = -INFINITY, min = INFINITY;
    for (int i = 0; i < n; i++) {
        max = MAX(max, x[i]);
        min = MIN(min, x[i]);
    }
    for (int i = 0; i < n; i++) {
        y[i] = (x[i] - min) / (max - min);
    }
    (void)smacofBernsteinBase(n, m, y, ordinal, z);
    (void)smacofPrintAnyMatrix(n, m, 15, 10, z);
    free(z);
    free(y);
    return EXIT_SUCCESS;
}
*/
