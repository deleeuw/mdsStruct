#include "smacof.h"

void smacofInitial(const double *delta, double *xold, const int *pn,
                   const int *pp) {
    int n = *pn, p = *pp, m = n * (n - 1) / 2, itmax = 100;
    double eps = 1e-10;
    int width = 10, precision = 6;
    double *cross = (double *)calloc((size_t)m, (size_t)sizeof(double));
    (void)smacofDoubleCenter(delta, cross, &n);
    if (DEBUG) {
        printf("cross\n");
        (void)smacofPrintSDCMatrix(cross, pn, &width, &precision);
    }
    (void)smacofSimultaneousIteration(cross, xold, pn, pp, &itmax, &eps);
}

void smacofDoubleCenter(const double *delta, double *cross, const int *pn) {
    int n = *pn, ij = 0;
    double tsum = 0.0, sum = 0.0, cell = 0.0;
    double *rsum = (double *)calloc((size_t)n, (size_t)sizeof(double));
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            sum += SQUARE(delta[PINDEX(i, j, n)]);
        }
        rsum[VINDEX(i)] = sum / ((double)n);
        tsum += sum;
    }
    tsum /= SQUARE((double)n);
    for (int j = 1; j <= (n - 1); j++) {
        for (int i = (j + 1); i <= n; i++) {
            ij = SINDEX(i, j, n);
            cell = SQUARE(delta[ij]);
            cross[ij] =
                -0.5 * (cell - rsum[VINDEX(i)] - rsum[VINDEX(j)] + tsum);
        }
    }
    free(rsum);
    return;
}

void smacofSimultaneousIteration(double *cross, double *xold, const int *pn,
                                 const int *pp, const int *itmax,
                                 const double *eps) {
    int n = *pn, p = *pp, np = n * p, itel = 1;
    double oldsum = 0.0, newsum = 0.0, maxdiff = 0.0;
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] = drand48();
        }
    }
    double *xnew = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *r = (double *)calloc((size_t)p, (size_t)sizeof(double));
    (void)smacofGramSchmidt(xold, r, &n, &p);
    oldsum = 0.0;
    while (true) {
        (void)smacofMultiplySDCMatrix(cross, xold, xnew, &n, &p);
        (void)smacofGramSchmidt(xnew, r, &n, &p);
        (void)smacofMaxConfDifference(xold, xnew, &maxdiff, &n, &p);
        newsum = 0.0;
        for (int i = 1; i <= p; i++) {
            newsum += r[VINDEX(i)];
        }
        if ((itel == *itmax) || ((newsum - oldsum) < *eps)) {
            break;
        }
        itel++;
        oldsum = newsum;
        (void)memcpy(xold, xnew, (size_t)(np * sizeof(double)));
    }
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] = xnew[MINDEX(i, s, n)] * sqrt(r[VINDEX(s)]);
        }
    }
    free(xnew);
    free(r);
    return;
}

void smacofGramSchmidt(double *x, double *r, int *pn, int *pp) {
    int ip, n = *pn, p = *pp, s = 1;
    while (s <= p) {
        for (int t = 1; t < s; t++) {
            double sum = 0.0;
            for (int i = 1; i <= n; i++) {
                sum += x[MINDEX(i, t, n)] * x[MINDEX(i, s, n)];
            }
            for (int i = 1; i <= n; i++) {
                x[MINDEX(i, s, n)] -= sum * x[MINDEX(i, t, n)];
            }
        }
        double sum = 0.0;
        for (int i = 1; i <= n; i++) {
            sum += x[MINDEX(i, s, n)] * x[MINDEX(i, s, n)];
        }
        sum = sqrt(sum);
        r[VINDEX(s)] = sum;
        for (int i = 1; i <= n; i++) {
            x[MINDEX(i, s, n)] /= sum;
        }
        s++;
    }
    return;
}

/*
int main() {
    double delta[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double cross[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double vinv[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double w[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    int n = 5;
    (void)smacofDoubleCenter(delta, cross, &n);
    (void)smacofVinverse(w, vinv, &n);
    (void)smacofPrintMM(cross, &n);
    (void)smacofPrintMM(vinv, &n);
    return (EXIT_SUCCESS);
}
 */