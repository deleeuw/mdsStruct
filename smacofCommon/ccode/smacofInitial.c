#include "smacof.h"

void smacofInitial(const double *delta, const double *weights, const int *irow,
                   const int *icol, double *xini, const int *pinit,
                   const int *pn, const int *pp, const int *pm) {
    int init = *pinit, np = *pn * *pp;
    int itmax_j = 100, itmax_d = 100, eps_j = 10, eps_d = 10;
    int width = 15, precision = 10;
    bool verbose_j = false, verbose_d = false;
    double *dini = (double *)calloc((size_t)np, (size_t)sizeof(double));
    switch (init) {
        case 1:
            (void)smacofInitTorgerson(delta, weights, irow, icol, xini, pn, pp,
                                      pm);
            break;
        case 2:
            (void)smacofInitMaximumSum(delta, weights, irow, icol, xini, pn, pp,
                                       pm, &itmax_j, &eps_j, &verbose_j);
            break;
        case 3:
            (void)smacofInitRandom(xini, pn, pp);
            break;
    }
    (void)smacofCenter(xini, pn, pp);
    (void)smacofDistance(xini, dini, irow, icol, pn, pp, pm);
    (void)smacofScale(delta, weights, dini, xini, pn, pp, pm);
    free(dini);
    return;
}

void smacofInitTorgerson(const double *delta, const double *weights,
                         const int *irow, const int *icol, double *xold,
                         const int *pn, const int *pp, const int *pm) {
    int n = *pn, p = *pp, m = *pm, nn = n * (n - 1) / 2, itmax = 100, eps = 10;
    int width = 15, precision = 10, nrow = 1;
    bool verbose = false;
    double *dimp = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *cross = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double sum = 0.0;
    for (int k = 1; k <= m; k++) {
        sum += weights[VINDEX(k)] * delta[VINDEX(k)];
    }
    for (int k = 1; k <= nn; k++) {
        dimp[VINDEX(k)] = sum;
    }
    for (int k = 1; k <= m; k++) {
        int ki = irow[VINDEX(k)];
        int kj = icol[VINDEX(k)];
        dimp[SINDEX(ki, kj, n)] = delta[VINDEX(k)];
    }
    (void)smacofDoubleCenter(dimp, cross, pn);
    (void)smacofJacobi(cross, evec, eval, pn, pp, &itmax, &eps, &verbose);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(fabs(eval[VINDEX(s)]));
        }
    }
    free(cross);
    free(evec);
    free(eval);
    free(dimp);
    return;
}

void smacofInitRandom(double *xini, const int *pn, const int *pp) {
    int n = *pn, p = *pp, np = n * p;
    for (int i = 1; i <= np; i++) {
        xini[VINDEX(i)] = drand48();
    }
    return;
}

void smacofInitMaximumSum(const double *delta, const double *weights,
                          const int *irow, const int *icol, double *xini,
                          const int *pn, const int *pp, const int *pm,
                          const int *pitmax_j, const int *peps_j,
                          const bool *pverbose_j) {
    int n = *pn, p = *pp, m = *pm, nn = n * (n - 1) / 2;
    double *a = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *b = (double *)calloc((size_t)(m + n), (size_t)sizeof(double));
    double *evec = (double *)calloc((size_t)SQUARE(n), (size_t)sizeof(double));
    double *eval = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int k = 1; k <= nn; k++) {
        b[VINDEX(k)] = 0.0;
    }
    for (int k = 1; k <= m; k++) {
        int kv = VINDEX(k), ik = irow[kv], jk = icol[kv];
        double cell = weights[kv] * SQUARE(delta[kv]);
        b[TINDEX(ik, jk, n)] -= cell;
        b[TINDEX(ik, ik, n)] += cell;
        b[TINDEX(jk, jk, n)] += cell;
    }
    (void)smacofJacobi(b, evec, eval, pn, pp, pitmax_j, peps_j, pverbose_j);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xini[MINDEX(i, s, n)] =
                evec[MINDEX(i, s, n)] * sqrt(eval[VINDEX(s)]);
        }
    }
    free(a);
    free(b);
    free(evec);
    free(eval);
    return;
}

void smacofDiagonalAdjust(const double *delta, const double *weights, double *x,
                          const int *pn, const int *pp, const int *pitmax,
                          const int *peps, const bool *pverbose) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2, itmax = *pitmax, itel = 1;
    double sum = 0.0, diff = 0.0, mul = 0.0,
           eps = pow((double)10.0, -(double)*peps), stress = 0.0;
    bool verbose = *pverbose;
    double *lbd = (double *)calloc((size_t)p, (size_t)sizeof(double));
    double *dia = (double *)calloc((size_t)p, (size_t)sizeof(double));
    double *tia = (double *)calloc((size_t)p, (size_t)sizeof(double));
    double *dwk = (double *)calloc((size_t)m, (size_t)sizeof(double));
    for (int s = 1; s <= p; s++) {
        lbd[VINDEX(s)] = 1.0;
        sum = 0.0;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (i == j) {
                    continue;
                }
                sum += weights[PINDEX(i, j, n)] *
                       SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
            }
        }
        dia[VINDEX(s)] = sum;
    }
    while (true) {
        for (int j = 1; j <= (n - 1); j++) {
            for (int i = (j + 1); i <= n; i++) {
                sum = 0.0;
                for (int s = 1; s <= p; s++) {
                    sum += SQUARE(lbd[VINDEX(s)] *
                                  (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]));
                }
                dwk[SINDEX(i, j, n)] = sqrt(sum);
            }
        }
        (void)smacofStress(delta, weights, dwk, &m, &stress);
        diff = 0.0;
        for (int s = 1; s <= p; s++) {
            sum = 0.0;
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= n; j++) {
                    if (i == j) {
                        continue;
                    }
                    sum += weights[PINDEX(i, j, n)] * delta[PINDEX(i, j, n)] *
                           SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) /
                           dwk[PINDEX(i, j, n)];
                }
            }
            tia[VINDEX(s)] = sum;
            mul = sum / dia[VINDEX(s)];
            diff = MAX(diff, fabs(1 - mul));
            lbd[VINDEX(s)] *= mul;
        }
        if (verbose) {
            printf("%3d %15.10f %15.10f %15.10f %15.10f\n", itel, stress, diff,
                   lbd[VINDEX(1)], lbd[VINDEX(2)]);
        }
        if ((itel == itmax) || (diff < eps)) {
            break;
        }
        itel++;
    }
    for (int s = 1; s <= p; s++) {
        for (int i = 1; i <= n; i++) {
            x[MINDEX(i, s, n)] *= lbd[VINDEX(s)];
        }
    }
    free(lbd);
    free(dia);
    free(tia);
    free(dwk);
    return;
}
