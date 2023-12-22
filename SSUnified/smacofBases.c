#include "smacof.h"

void checkIncreasing(const double *innerknots, const double *lowend,
                     const double *highend, const int *ninner, bool *fail) {
    *fail = false;
    if (*lowend >= innerknots[VINDEX(1)]) {
        *fail = true;
        return;
    }
    if (*highend <= innerknots[VINDEX(*ninner)]) {
        *fail = true;
        return;
    }
    for (int i = 1; i < *ninner; i++) {
        if (innerknots[i] <= innerknots[i - 1]) {
            *fail = true;
            return;
        }
    }
}

void extendPartition(const double *innerknots, const int *multiplicities,
                     const int *order, const int *ninner, const double *lowend,
                     const double *highend, double *extended) {
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

void bisect(const double *x, const double *knots, const int *lowindex,
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

void bsplines(const double *x, const double *knots, const int *order,
              const int *nknots, int *index, double *q) {
    int lowindex = 1, highindex = *nknots, m = *order, j, jp1;
    double drr, dll, saved, term;
    double *dr = (double *)calloc((size_t)m, sizeof(double));
    double *dl = (double *)calloc((size_t)m, sizeof(double));
    (void)bisect(x, knots, &lowindex, &highindex, index);
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
        j  = jp1;
    }
    free(dr);
    free(dl);
    return;
}

void bsplineBasis(const double *x, const double *knots, const int *order,
                  const int *nknots, const int *nvalues, double *result) {
    int m = *order, l = 0;
    double *q = (double *)calloc((size_t)m + 1, sizeof(double));
    for (int i = 1; i <= *nvalues; i++) {
        (void)bsplines(x + VINDEX(i), knots, order, nknots, &l, q);
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
    //(void)smacofPrintAnyMatrix(stdout, n, m, 15, 10, z);
    return;
}

/*
double innerknots[3] = {.1, .3, .7};
int multiplicities[3] = {1, 1, 1};
const int order = 2;
const int ninner = 3;
const double lowend = 0;
const double highend = 1;
const int nknots = 7;
const int nvalues = 10;
double extended[7] = {0};
double x[10] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};
double result[50] = {0};

int main(void) {
    (void)smacofExtendPartition(innerknots, multiplicities, &order, &ninner,
                                &lowend, &highend, extended);
    (void)smacofPrintVector(stdout, nknots, 15, 10, extended);
    (void)smacofBsplineBasis(x, extended, &order, &nknots, &nvalues, result);
    (void)smacofPrintVector(stdout, 50, 15, 10, result);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%15.10f ", result[MINDEX(i, j, 10)]);
        }
        printf("\n");
    }
    return EXIT_SUCCESS;
}
*/