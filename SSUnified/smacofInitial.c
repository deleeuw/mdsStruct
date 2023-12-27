#include "smacof.h"

void smacofInitial(const int n, const int p, const int typeinit,
                   const bool weights, double **delta, double **w,
                   double **xold) {
    switch (typeinit) {
        case 1:
            (void)smacofInitTorgerson(n, p, delta, xold);
            break;
        case 2:
            (void)smacofInitMaximumSum(n, p, weights, delta, w, xold);
            break;
        case 3:
            (void)smacofInitRandom(n, p, xold);
            break;
    }
    (void)smacofCenter(n, p, xold);
    return;
}

void smacofInitTorgerson(const int n, const int p, double **delta,
                         double **xold) {
    int jitmax = 100, jeps = 15;
    bool jverbose = false;
    double **cross = smacofMakeSymmetricMatrix(n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeVector(n);
    (void)smacofDoubleCenter(n, delta, cross);
    (void)smacofJacobi(n, p, cross, evec, eval, jitmax, jeps, jverbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            xold[i][s] = fac * evec[i][s];
        }
    }
    (void)smacofFreeMatrix(n, cross);
    (void)smacofFreeMatrix(n, evec);
    (void)smacofFreeVector(eval);
    return;
}

void smacofInitMaximumSum(const int n, const int p, const bool weights,
                          double **delta, double **w, double **xold) {
    int jitmax = 100, jeps = 15;
    bool jverbose = false;
    double **bmat = smacofMakeSymmetricMatrix(n);
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeVector(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int ij = MAX(i, j);
            int ji = MIN(i, j);
            double cell = SQUARE(delta[ij][ji]);
            if (weights) {
                cell *= w[ij][ji];
            }
            bmat[ij][ji] = -cell;
            sum += cell;
        }
        bmat[i][i] = sum;
    }
    (void)smacofJacobi(n, p, bmat, evec, eval, jitmax, jeps, jverbose);
    for (int s = 0; s < p; s++) {
        double fac = sqrt(fabs(eval[s]));
        for (int i = 0; i < n; i++) {
            xold[i][s] = fac * evec[i][s];
        }
    }
    (void)smacofFreeMatrix(n, bmat);
    (void)smacofFreeMatrix(n, evec);
    (void)smacofFreeVector(eval);
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
