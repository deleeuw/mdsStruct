#include "smacofWeighted.h"

void smacofWeightedInitial(const int n, const int p, const int init,
                           double **delta, double **weights, double **xini) {
    switch (init) {
        case 1:
            (void)smacofWeightedInitTorgerson(n, p, delta, xini);
            break;
        case 2:
            (void)smacofWeightedInitMaximumSum(n, p, delta, weights, xini);
            break;
        case 3:
            (void)smacofInitRandom(n, p, xini);
            break;
    }
    (void)smacofCenter(n, p, xini);
    return;
}

void smacofWeightedInitTorgerson(const int n, const int p, double **delta,
                                 double **xold) {
    int itmax = 100, eps = 15;
    bool verbose = false;
    double **cross = smacofMakeAnyMatrix(n, n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    (void)smacofDoubleCenter(n, delta, cross);
    (void)smacofJacobi(n, p, cross, evec, eval, itmax, eps, verbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            xold[i][s] = fac * evec[i][s];
        }
    }
    (void)smacofFreeAnyMatrix(n, cross);
    (void)smacofFreeAnyMatrix(n, evec);
    (void)smacofFreeAnyVector(eval);
    return;
}

void smacofWeightedInitMaximumSum(const int n, const int p, double **delta,
                                  double **weights, double **xold) {
    int itmax = 100, eps = 15;
    bool verbose = false;
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double cell = weights[i][j] * SQUARE(delta[i][j]);
            bmat[i][j] -= cell;
            bmat[j][i] -= cell;
            bmat[i][i] += cell;
            bmat[j][j] += cell;
        }
    }
    (void)smacofJacobi(n, p, bmat, evec, eval, itmax, eps, verbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            xold[i][s] = fac * evec[i][s];
        }
    }
    (void)smacofFreeAnyMatrix(n, bmat);
    (void)smacofFreeAnyMatrix(n, evec);
    (void)smacofFreeAnyVector(eval);
    return;
}
