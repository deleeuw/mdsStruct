#include "../Include/smacof.h"

void smacofUnweightedMakeBMatrix(const size_t n, const double delta[n][n],
                                 const double dmat[n][n], double bmat[n][n]) {
  double cell = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            bmat[i][j] = -0.0;
        }
    }
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double dfix = dmat[i][j];
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = delta[i][j] / (dfix * (double)n);
            }
            bmat[i][j] -= cell;
            bmat[j][i] -= cell;
            bmat[i][i] += cell;
            bmat[j][j] += cell;
        }
    }
    return;
}

void smacofUnweightedGuttman(const size_t n, const int p,
                             const double bmat[n][n], const double xold[n][p],
                             double xnew[n][p]) {
    (void)smacofMultiplySymmetricMatrix(n, p, bmat, xold, xnew);
    return;
}

void smacofUnweightedMakeStress(const size_t n, const double delta[n][n],
                                const double dist[n][n], double *stress) {
    double sum = 0.0, m = n * (n - 1);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            sum += SQUARE(delta[i][j] - dist[i][j]);
        }
    }
    *stress = sum / (2.0 * (double)m);
    return;
}
