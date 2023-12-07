#include "../Common/smacofCommon.h"

int main(void) {
    int n = 4, p = 2, width = 15, precision = 10;
    double **weights = smacofMakeAnyMatrix(n, n);
    double **delta = smacofMakeAnyMatrix(n, n);
    double **dmat = smacofMakeAnyMatrix(n, n);
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **vmat = smacofMakeAnyMatrix(n, n);
    double **xmat = smacofMakeAnyMatrix(n, p);
    double ****hessian = NULL;
    xmat[0][0] = xmat[1][0] = xmat[0][1] = xmat[3][1] = 1.0 / 2.0;
    xmat[2][0] = xmat[3][0] = xmat[1][1] = xmat[2][1] = -1.0 / 2.0;
    dmat[1][0] = dmat[0][1] = dmat[2][1] = dmat[1][2] = dmat[3][0] =
        dmat[0][3] = dmat[3][2] = dmat[2][3] = 1.0;
    dmat[2][0] = dmat[0][2] = dmat[3][1] = dmat[1][3] = sqrt(2.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                weights[i][i] = delta[i][i] = dmat[i][j] = 0.0;
            } else {
                weights[i][j] = 1;
                delta[i][j] = fabs((double)(i - j)) + 1;
                vmat[i][j] = -1;
                bmat[i][j] = -weights[i][j] * delta[i][j] / dmat[i][j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double sumv = 0.0, sumb = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            sumv += vmat[i][j];
            sumb += bmat[i][j];
        }
        vmat[i][i] = -sumv;
        bmat[i][i] = -sumb;
    }
    (void)smacofPrintAnyMatrix(n, n, width, precision, weights);
    (void)smacofPrintAnyMatrix(n, n, width, precision, delta);
    (void)smacofPrintAnyMatrix(n, n, width, precision, dmat);
    (void)smacofPrintAnyMatrix(n, n, width, precision, bmat);
    (void)smacofPrintAnyMatrix(n, p, width, precision, xmat);
    for (int s = 0; s < p; s++) {
        for (int t = 0; t < p; t++) {
            hessian[s][t] = smacofMakeAnyMatrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    hessian[s][t][i][j] = (double)(s + t + i + j);
                }
            }
        }
    }
}

void smacofWeightedHessian(const int n, const int p, double **delta,
                           double **weights, double **dmat, double **bmat,
                           double **vmat, double **x, double ****hessian) {
    for (int s = 0; s < p; s++) {
        for (int t = 0; t < p; t++) {
            double **work = smacofMakeAnyMatrix(n, n);
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    double cell =
                        weights[i][j] * delta[i][j] / THIRD(dmat[i][j]);
                    cell *= (x[i][s] - x[j][s]);
                    cell *= (x[i][t] - x[j][t]);
                    hessian[s][t][i][j] = -cell;
                    sum += cell;
                }
                hessian[s][t][i][i] = sum;
            }
        }
    }
    for (int s = 0; s < p; s++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                hessian[s][s][i][j] += (vmat[i][j] - bmat[i][j]);
            }
        }
    }
}