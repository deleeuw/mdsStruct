#include "smacof.h"

void smacofInitElegant(const double *delta, const double *weights,
                       const int *pn, const int *pm, const double *plbd,
                       double *proot, const int *pitmax, const double *peps,
                       const bool *pverbose) {
    int n = *pn, m = *pm, itmax = *pitmax, mm = m * (m + 1) / 2,
        nn = n * (n - 1) / 2;
    int ik = 0, jk = 0, il = 0, jl = 0;
    int width = 6, precision = 4;
    bool verbose = *pverbose;
    double eps = *peps, cell = 0.0, lbd = *plbd, root = *proot;
    int *ii = (int *)calloc((size_t)nn, (size_t)sizeof(double));
    int *jj = (int *)calloc((size_t)nn, (size_t)sizeof(double));
    double *u = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *cross = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    (void)smacofMakeIIandJJ(pn, ii, jj);
    (void)smacofDoubleCenter(delta, cross, pn);
    for (int t = 1; t <= m; t++) {
        int i = ii[VINDEX(t)];
        int j = jj[VINDEX(t)];
        for (int s = t; s <= m; s++) {
            int k = ii[VINDEX(s)];
            int l = jj[VINDEX(s)];
            cell = SQUARE(KDELTA(i, k) - KDELTA(j, k) - KDELTA(i, l) +
                          KDELTA(j, l));
            u[TINDEX(s, t, m)] =
                cell * sqrt(weights[VINDEX(s)] * weights[VINDEX(t)]);
        }
    }
    (void)smacofPrintLTMatrix(u, &m, &width, &precision);
    (void)smacofPerronRoot(u, pm, plbd, proot, pitmax, peps, pverbose);
    printf("%15.10f\n", *proot);
    free(u);
    free(ii);
    free(jj);
}

void smacofPerronRoot(double *a, const int *pn, const double *plbd,
                      double *proot, const int *pitmax, const double *peps,
                      const bool *verbose) {
    int n = *pn, itel = 1, itmax = *pitmax;
    double lbd = *plbd, eps = *peps, root = *proot;
    double *r = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double *b = (double *)calloc((size_t)(n * n), (size_t)sizeof(double));
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            double h = a[TINDEX(i, j, n)];
            if (i == j) {
                b[MINDEX(i, j, n)] = h + lbd;
            } else {
                b[MINDEX(i, j, n)] = h;
                b[MINDEX(j, i, n)] = h;
            }
        }
    }
    while (true) {
        double rmin = INFINITY, rmax = 0.0;
        for (int i = 1; i <= n; i++) {
            double sum = 0.0;
            for (int j = 1; j <= n; j++) {
                sum += b[MINDEX(i, j, n)];
            }
            r[VINDEX(i)] = sum;
            rmin = MIN(rmin, sum);
            rmax = MAX(rmax, sum);
        }
        root = ((rmin + rmax) / 2.0) - lbd;
        if (verbose) {
            printf("itel %3d rmin %15.10f rmax %15.10f root %15.10f\n", itel,
                   rmin, rmax, root);
        }
        if ((itel == itmax) || ((rmax - rmin) < eps)) {
            break;
        }
        itel++;
        for (int j = 1; j <= n; j++) {
            for (int i = 1; i <= n; i++) {
                b[MINDEX(i, j, n)] *= r[VINDEX(j)] / r[VINDEX(i)];
            }
        }
    }
    *proot = root;
    free(r);
    free(b);
}

int main() {
    double delta[10] = {1, 2, 3, 4, 5, 6, 6, 6, 6, 6};
    double weights[10] = {1, 1, 1, 1, 1, 1, 7, 8, 9, 10};
    int n = 5, m = 10, itmax = 0;
    double eps = 1e-10, lbd = 1.0, root = 0.0;
    bool verbose = true;
    (void)smacofElegant(delta, weights, &n, &m, &lbd, &root, &itmax, &eps,
                        &verbose);
}
