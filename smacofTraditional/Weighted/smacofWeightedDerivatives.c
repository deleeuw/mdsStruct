#include "../Common/smacofCommon.h"

/*
void smacofWeightedGradient(const int *nconst double *delta, const double
*weights, const double *vinv, const double *dold, const double *xold, double
*xnew, double *gradient, ) { int n = *pn, p = *pp, np = n * p, m = n * (n - 1) /
2; double *bmat = (double *)calloc((size_t)m, (size_t)sizeof(double)); double
*ymat = (double *)calloc((size_t)np, (size_t)sizeof(double)); for (int i = 1; i
<= m; i++) { double dfix = dold[VINDEX(i)]; int iv = VINDEX(i); if (dfix <
1e-15) { bmat[iv] = 0.0; } else { bmat[iv] = weights[iv] * delta[iv] / dfix;
        }
    }
    (void)smacofMultiplySymmetricMatrix(bmat, xold, ymat, pn, pp);
    (void)smacofMultiplySymmetricMatrix(vinv, ymat, xnew, pn, pp);
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            gradient[MINDEX(i, s, n)] =
                xold[MINDEX(i, s, n)] - xnew[MINDEX(i, s, n)];
        }
    }
    free(bmat);
    free(ymat);
    return;
}
*/

void smacofWeightedHessian(const int n, const int p, double **delta,
                           double **weights, double **dmat, double **bmat,
                           double **vmat, double **x, double ****hessian) {
    for (int s = 0; s < p; s++) {
        for (int t = 0; t < p; t++) {
            hessian[s][t] = smacofMakeAnyMatrix(n, n);
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

// The Hessian is np by np, and it is the SDCL matrix kronecker(I_p, V - B) + H
// First make H. The diagonal blocks are SDCL, the off-diagonal blocks regular
// matrices. Then add V and B to diagonal blocks
// It makes sense to use an array of p(p+1)/2 pointers, each pointing to a block

/*
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const int *pn,
                           const int *pp) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2, q = p * (p - 1) / 2;
    double **diagH = (double **)calloc((size_t)p, (size_t)sizeof(double *));
    double **nonDiagH = (double **)calloc((size_t)q, (size_t)sizeof(double *));
    for (int s = 1; s <= p; s++) {
        diagH[VINDEX(s)] = (double *)calloc((size_t)m, (size_t)sizeof(double));
    }
    // make nonDiagH[SINDEX(i, j, p)]
    for (int t = 1; t <= (p - 1); t++) {
        for (int s = (t + 1); s <= p; s++) {
            int st = SINDEX(s, t, p);
            nonDiagH[st] =
                (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
            for (int i = 1; i <= n; i++) {
                int is = MINDEX(i, s, n);
                int it = MINDEX(i, t, n);
                for (int j = 1; j <= n; j++) {
                    if (i == j) {
                        continue;
                    }
                    int ij = PINDEX(i, j, n);
                    int js = MINDEX(j, s, n);
                    int jt = MINDEX(j, t, n);
                    double f = weights[ij] * delta[ij] / THIRD(dmat[ij]);
                    nonDiagH[st][ij] =
                        f * (xconf[is] - xconf[js]) * (xconf[it] - xconf[jt]);
                }
            }
        }
    }
    // make diagH[VINDEX(i)]
    for (int s = 1; s <= p; s++) {
        int ss = VINDEX(s);
        diagH[VINDEX(s)] = (double *)calloc((size_t)m, (size_t)sizeof(double));
        for (int j = 1; j <= (n - 1); j++) {
            int js = MINDEX(j, s, n);
            for (int i = (j + 1); i <= n; i++) {
                int ij = SINDEX(i, j, n);
                int is = MINDEX(i, s, n);
                double f = weights[ij] * delta[ij] / THIRD(dmat[ij]);
                diagH[ss][ij] = f * SQUARE(xconf[is] - xconf[js]);
                diagH[ss][ij] += (v[ij] - bmat[ij]);
            }
        }
    }
    for (int s = 1; s <= p; s++) {
        free(diagH[VINDEX(s)]);
    }
    free(diagH);
    for (int t = 1; t <= (p - 1); t++) {
        for (int s = (t + 1); s <= p; s++) {
            int st = SINDEX(s, t, p);
            free(nonDiagH[st]);
        }
    }
    free(nonDiagH);
    return;
}
*/
