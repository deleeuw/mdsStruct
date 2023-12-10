#include "smacof.h"


void smacofInitial(const int n, const int p, const int init, const bool weights,
                   double **delta, double **w, double **xini) {
    switch (init) {
        case 1:
            (void)smacofInitTorgerson(n, p, weights, delta, w, xini);
            break;
        case 2:
            (void)smacofInitMaximumSum(n, p, weights, delta, w, xini);
            break;
        case 3:
            (void)smacofInitRandom(n, p, xini);
            break;
    }
    (void)smacofCenter(n, p, xini);
    return;
}

void smacofInitTorgerson(const int n, const int p, const bool weights,
                         double **delta, double **w, double **xold) {
    int itmax = 100, eps = 15;
    bool verbose = false;
    double **cross = smacofMakeSymmetricMatrix(n);
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

void smacofInitMaximumSum(const int n, const int p, const bool weights,
                          double **delta, double **w, double **xold) {
    int itmax = 100, eps = 15;
    bool verbose = false;
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double cell = w[i][j] * SQUARE(delta[i][j]);
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

void smacofInitRandom(const int n, const int p, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            x[i][j] = drand48();
        }
    }
    return;
}