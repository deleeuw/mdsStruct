#include "smacof.h"

double smacofStress(const int n, const bool weights, double **delta, double **w,
                    double **dmat) {
    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j <= (i - 1); j++) {
            double fac = SQUARE(delta[i][j] - dmat[i][j]);
            if (weights) {
                fac *= w[i][j];
            }
            sum += fac;
        }
    }
    return sum;
}

void smacofGuttmanTransform(const int n, const int p, const bool weights,
                            double **delta, double **dmat, double **w,
                            double **vinv, double **xold, double **xnew) {
    double **xaux = smacofMakeAnyMatrix(n, p);
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                int ij = MAX(i, j);
                int ji = MIN(i, j);
                double rat = (delta[ij][ji] / dmat[ij][ji]);
                if (weights) {
                    rat *= w[ij][ji];
                }
                sum += rat * (xold[i][s] - xold[j][s]);
            }
            xaux[i][s] = -sum;
        }
    }
    if (weights) {
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        continue;
                    }
                    int ij = MAX(i, j);
                    int ji = MIN(i, j);
                    double rat = vinv[ij][ji];
                    sum += rat * (xaux[i][s] - xaux[j][s]);
                }
                xnew[i][s] = -sum;
            }
        }
    } else {
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < n; i++) {
                xnew[i][s] = xaux[i][s] / (double)n;
            }
        }
    }
    free(xaux);
    return;
}

/*
int main(void) {
    int n = 4, p = 2;
    bool weights = false;
    double **delta = smacofMakeSymmetricMatrix(n);
    double **dmat = smacofMakeSymmetricMatrix(n);
    double **w = smacofMakeSymmetricMatrix(0);
    double **vinv = smacofMakeSymmetricMatrix(0);
    double **xold = smacofMakeAnyMatrix(n, p);
    double **xnew = smacofMakeAnyMatrix(n, p);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                // w[i][j] = 0.0;
                delta[i][j] = 0.0;
                dmat[i][j] = 0.0;
                // vinv[i][j] = 3.0/16.0;
            } else {
                // w[i][j] = 1.0;
                // vinv[i][j] = -1.0/16.0;
            }
        }
    }
    xold[0][0] = -1.0;
    xold[0][1] = -1.0;
    xold[1][0] = 1.0;
    xold[1][1] = -1.0;
    xold[2][0] = 1.0;
    xold[2][1] = 1.0;
    xold[3][0] = -1.0;
    xold[3][1] = 1.0;
    dmat[1][0] =
        sqrt(SQUARE(xold[1][0] - xold[0][0]) + SQUARE(xold[1][1] - xold[0][1]));
    dmat[2][0] =
        sqrt(SQUARE(xold[2][0] - xold[0][0]) + SQUARE(xold[2][1] - xold[0][1]));
    dmat[2][1] =
        sqrt(SQUARE(xold[2][0] - xold[1][0]) + SQUARE(xold[2][1] - xold[1][1]));
    dmat[3][0] =
        sqrt(SQUARE(xold[3][0] - xold[0][0]) + SQUARE(xold[3][1] - xold[0][1]));
    dmat[3][1] =
        sqrt(SQUARE(xold[3][0] - xold[1][0]) + SQUARE(xold[3][1] - xold[1][1]));
    dmat[3][2] =
        sqrt(SQUARE(xold[3][0] - xold[2][0]) + SQUARE(xold[3][1] - xold[2][1]));
    delta[1][0] = 1.0;
    delta[2][0] = 2.0;
    delta[2][1] = 3.0;
    delta[3][0] = 4.0;
    delta[3][1] = 5.0;
    delta[3][2] = 6.0;
    (void)smacofPrintSymmetricMatrix(n, 15, 10, delta);
    //(void)smacofPrintSymmetricMatrix(n, 15, 10, weights);
    (void)smacofPrintSymmetricMatrix(n, 15, 10, dmat);
    //(void)smacofPrintSymmetricMatrix(n, 15, 10, vinv);
    (void)smacofPrintAnyMatrix(n, p, 15, 10, xold);
    //(void)smacofPrintAnyMatrix(n, p, 15, 10, xnew);
    (void)smacofGuttmanTransform(n, p, weights, delta, dmat, weights, vinv,
                                 xold, xnew);
    (void)smacofPrintAnyMatrix(n, p, 15, 10, xnew);
}
*/
