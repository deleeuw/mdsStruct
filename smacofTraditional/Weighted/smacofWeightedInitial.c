#include "smacofWeighted.h"

void smacofWeightedInitial(const int n, const int p, const int init, double **delta, double **weights,
                           double **xini) {
    double *dini = (double *)calloc((size_t)m, (size_t)sizeof(double));
    switch (init) {
        case 1:
            (void)smacofWeightedInitTorgerson(delta, weights, irow, icol, xini,
                                              pn, pp, pm);
            break;
        case 2:
            (void)smacofWeightedInitMaximumSum(delta, weights, irow, icol, xini,
                                               pn, pp, pm);
            break;
        case 3:
            (void)smacofInitRandom(xini, pn, pp);
            break;
        case 4:;
            break;
    }
    (void)smacofCenter(xini, pn, pp);
    (void)smacofDistance(xini, dini, pn, pp);
    (void)smacofWeightedScale(delta, weights, dini, xini, pn, pp, pm);
    free(dini);
    return;
}

void smacofWeightedInitTorgerson(const int n, const int p, double **delta, double **weights, double **xold) {
    bool verbose = false;
    double *dimp = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *cross = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double sum = 0.0;
    for (int k = 0; k < m; k++) {
        sum += weights[k] * delta[k];
    }
    for (int k = 0; k < nn; k++) {
        dimp[k] = sum;
    }
    for (int k = 0; k < m; k++) {
        int ki = irow[k];
        int kj = icol[k];
        dimp[ki][kj] = delta[k];
    }
    (void)smacofDoubleCenter(dimp, cross, pn);
    (void)smacofJacobi(cross, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xold[i][s] = evec[i][s] * sqrt(fabs(eval[s]));
        }
    }
    free(cross);
    free(evec);
    free(eval);
    free(dimp);
    return;
}

void smacofWeightedInitMaximumSum(const int n, const int p,
                                  double **delta, double **weights,
                                  double **xini) {
    bool verbose = false;
    double *a = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *b = (double *)calloc((size_t)(m + n), (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int k = 0; k < nn; k++) {
        b[k] = 0.0;
    }
    for (int k = 0; k < m; k++) {
        int ik = irow[k], jk = icol[k];
        double cell = weights[k] * SQUARE(delta[k]);
        b[i][j] -= cell;
        b[i][i] += cell;
        b[j][j] += cell;
    }
    (void)smacofJacobi(b, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            xini[i][s] = evec[i][s] * sqrt(eval[s]);
        }
    }
    free(a);
    free(b);
    free(evec);
    free(eval);
    return;
}
