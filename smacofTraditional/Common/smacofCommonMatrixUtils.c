#include "smacofCommon.h"

// make an n x p matrix

double **smacofMakeAnyMatrix(const int n, const int p) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc((size_t)p, sizeof(double));
    }
    return x;
}

double *smacofMakeAnyVector(const int n) {
    double *x = (double *)calloc((size_t)n, sizeof(double *));
    return x;
}

// make a lower digonal matrix of order n
// a strict lower diagonal matrix of order n is
// a lower diagonal matrix of order n - 1

double **smacofMakeLowerDiagonalMatrix(const int n) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc((i + 1), sizeof(double));
    }
    return x;
}

void smacofCopyAnyMatrix(const int n, const int p, double **x, double **y) {
    for (int i = 0; i < n; i++) {
        memcpy(y[i], x[i], p * sizeof(double));
    }
    return;
}

// free an n x p matrix

void smacofFreeAnyMatrix(const int n, double **x) {
    for (int i = 0; i < n; i++) {
        free(x[i]);
    }
    free(x);
    return;
}

void smacofFreeAnyVector(double *x) {
    free(x);
    return;
}
