#include "smacof.h"

// Computes the Moore-Penrose inverse of an SDCL
// matrix of order n and rank n-1 from only the elements
// below the diagonal

void smacofMPInverseSDCLMatrix(const double *v, double *vinv, const int *pn) {
    int n = *pn, m = n * (n - 1) / 2, l = n * (n + 1) / 2;
    double add = 1.0 / ((double)n);
    double *vadd = (double *)calloc((size_t)l, (size_t)sizeof(double));
    double *vmid = (double *)calloc((size_t)l, (size_t)sizeof(double));
    (void)smacofAddSDCLDiagonal(v, vadd, pn);
    for (int k = 1; k <= l; k++) {
        vadd[VINDEX(k)] += add;
    }
    (void)smacofInvertPDMatrix(vadd, vmid, pn);
    for (int k = 1; k <= l; k++) {
        vmid[VINDEX(k)] -= add;
    }
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            vinv[SINDEX(i, j, n)] = vmid[TINDEX(i, j, n)];
        }
    }
    free(vadd);
    free(vmid);
    return;
}

void smacofMultiplySDCLMatrix(const double *a, const double *x, double *y,
                              const int *pn, const int *pp) {
    int n = *pn, p = *pp;
    for (int s = 1; s <= p; s++) {
        for (int i = 1; i <= n; i++) {
            double sum = 0.0;
            for (int j = 1; j <= n; j++) {
                if (j == i) {
                    continue;
                }
                int ij = PINDEX(i, j, n);
                sum += -a[ij] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
            y[MINDEX(i, s, n)] = sum;
        }
    }
    return;
}

// smacofAddSDCLDiagonal() takes the n(n -1)/2 elements of the strict lower
// triangular part of an SDCL matrix in column-major storage and returns the n(n
// + 1)/2 elements of the lower triangular part, including the diagonal, again
// in column-major storage.

void smacofAddSDCLDiagonal(const double *a, double *b, const int *pn) {
    int n = *pn;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            sum += a[PINDEX(i, j, n)];
        }
        rsum[VINDEX(i)] = sum;
    }
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            int ij = TINDEX(i, j, n);
            if (i == j) {
                b[ij] = -rsum[VINDEX(i)];
            } else {
                b[ij] = a[SINDEX(i, j, n)];
            }
        }
    }
    free(rsum);
    return;
}

/*
int main() {
    double v[3] = {-1, -1, -1};
    double vinv[3] = {0, 0, 0};
    int n = 3, width = 15, precision = 10;
    (void)smacofPrintSDCLMatrix(v, &n, &width, &precision);
    (void)smacofMPInverseSDCLMatrix(v, vinv, &n);
    (void)smacofPrintSDCLMatrix(vinv, &n, &width, &precision);
    return(EXIT_SUCCESS);
}
*/