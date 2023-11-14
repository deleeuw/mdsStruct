#include "../Include/smacof.h"

void smacofUnweightedMakeBMatrix(const unsigned n, const double (*delta)[n][n],
                                 const double (*dmat)[n][n], double (*bmat)[n][n]) {
    double cell = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            (*bmat)[i][j] = -0.0;
        }
    }
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            double dfix = (*dmat)[i][j];
            if (dfix < 1e-15) {
                cell = -0.0;
            } else {
                cell = (*delta)[i][j] / (dfix * (double)n);
            }
            (*bmat)[i][j] -= cell;
            (*bmat)[j][i] -= cell;
            (*bmat)[i][i] += cell;
            (*bmat)[j][j] += cell;
        }
    }
    return;
}

void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double (*bmat)[n][n], const double (*xold)[n][p],
                             double (*xnew)[n][p]) {
    (void)smacofMultiplySymmetricMatrix(n, p, bmat, xold, xnew);
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < p; j++) {
        (*xnew)[i][j] /= (double)n;
      }
    }
    return;
}

double smacofUnweightedMakeStress(const unsigned n, const double (*delta)[n][n],
                                const double (*dist)[n][n]) {
    double sum = 0.0, stress = 0.0;
    unsigned m = n * (n - 1) / 2;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            sum += SQUARE((*delta)[i][j] - (*dist)[i][j]);
        }
    }
    stress = sum / (2.0 * (double)m);
    return stress;
}
