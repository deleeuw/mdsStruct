#include "Include/smacof.h"

/*
int n = 20;
int m = 4;
double x[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 3, 1, 1, 4, 5, 3, 2, 1, 2, 3};
int width = 10;
int precision = 6;
*/

void smacofBernsteinBase(const int n, const int m, const double *x,
                         double (*z)[n][m]);

/*
int main() {
    double(*z)[n][m] = malloc(sizeof *z);
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
    (void)smacofBernsteinBase(n, m, y, z);
    (void)smacofPrintAnyMatrix(n, m, width, precision, z);
    free(z);
    free(y);
}

*/

void smacofBernsteinBase(const int n, const int m, const double *y,
                         double (*z)[n][m]) {
    double fac = 0.0, rat = 0.0;
    for (int i = 0; i < n; i++) {
        (*z)[i][0] = pow(1.0 - y[i], (double)(m - 1));
        for (int j = 1; j < m; j++) {
            rat = ((double)(m - j)) / ((double)j);
            if (y[i] == 0.0) {
                (*z)[i][j] = 0.0;
                continue;
            }
            if (y[i] == 1.0) {
                if (j == (m - 1)) {
                    (*z)[i][j] = 1.0;
                } else {
                    (*z)[i][j] = 0.0;
                }
                continue;
            }
            fac = y[i] / (1.0 - y[i]);
            (*z)[i][j] = rat * fac * (*z)[i][j - 1];
        }
    }
}

void smacofBSplineBase(void) {

}

void smacofMonomialBase(void) {

}