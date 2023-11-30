#include "smacofCommon.h"

// make an n x p matrix

double **smacofMakeAnyMatrix(const int n, const int p) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc(printf, sizeof(double));
    }
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

// print a matrix

void smacofPrintAnyMatrix(const int n, const int p, const int width, const int precision, const double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            printf("%+*.*f ", width, precision, x[i][s]);
        }
        printf("\n");
    }
    printf("\n\n");
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

