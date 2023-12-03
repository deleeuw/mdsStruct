#include "smacofUnweighted.h"

void smacofUnweightedInitial(const int n, const int p, const int init,
                             double **delta, double **xini) {
    double **dini = smacofMakeAnyMatrix(n, n);
    assert(!(dini == NULL));
    switch (init) {
        case 1:
            (void)smacofUnweightedInitTorgerson(n, p, delta, xini);
            break;
        case 2:
            (void)smacofUnweightedInitMaximumSum(n, p, delta, xini);
            break;
        case 3:
            (void)smacofInitRandom(n, p, xini);
            break;
    }
    (void)smacofCenter(n, p, xini);
    (void)smacofDistance(n, p, xini, dini);
    (void)smacofFreeAnyMatrix(n, dini);
    return;
}

void smacofUnweightedInitTorgerson(const int n, const int p, double **delta,
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

void smacofUnweightedInitMaximumSum(const int n, const int p, double **delta,
                                    double **xini) {
    int itmax = 100, eps = 10;
    bool verbose = false;
    double **b = smacofMakeAnyMatrix(n, n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = 0.0;
        }
    }
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double cell = SQUARE(delta[i][j]);
            b[i][j] -= cell;
            b[j][i] -= cell;
            b[i][i] += cell;
            b[j][j] += cell;
        }
    }
    (void)smacofPrintAnyMatrix(n, n, 15, 10, b);
    (void)smacofJacobi(n, p, b, evec, eval, itmax, eps, verbose);
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < n; s++) {
            xini[i][s] = evec[i][s] * sqrt(eval[s]);
        }
    }
    (void)smacofFreeAnyMatrix(n, b);
    (void)smacofFreeAnyMatrix(n, evec);
    (void)smacofFreeAnyVector(eval);
    return;
}
