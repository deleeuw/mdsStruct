#include "smacof.h"

void smacofCCD(const int n, const int ncol, double *dmatvec, double *wvec,
               double *bcoef, double *dhatvec, double **basis, const int ditmax,
               const int ieps3, const bool dverbose, const bool weights,
               const bool ordinal) {
    int itel = 1, m = n * (n - 1) / 2;
    double *s = (double *)calloc((size_t)ncol, (size_t)sizeof(double));
    double *r = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double deps = pow(10.0, -(double)ieps3), sold = 0.0, snew = 0.0;
    for (int j = 0; j < ncol; j++) {
        double sum = 0.0;
        for (int i = 0; i < m; i++) {
            if (weights) {
                sum += wvec[i] * SQUARE(basis[i][j]);
            } else {
                sum += SQUARE(basis[i][j]);
            }
        }
        if (sum == 0.0) {
            s[j] = 1.0;
        } else {
            s[j] = sum;
        }
    }
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < m; i++) {
        r[i] = dhatvec[i] - dmatvec[i];
        if (weights) {
            sum1 += wvec[i] * SQUARE(dhatvec[i] - dmatvec[i]);
            sum2 += wvec[i] * SQUARE(dmatvec[i]);
        } else {
            sum1 += SQUARE(dhatvec[i] - dmatvec[i]);
            sum2 += SQUARE(dmatvec[i]);
        }
    }
    sold = sum1 / sum2;
    snew = sold;
    while (true) {
        double meps = 0.0, chng = 0.0;
        for (int j = 0; j < ncol; j++) {
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                if (weights) {
                    sum += wvec[i] * basis[i][j] * r[i];
                } else {
                    sum += basis[i][j] * r[i];
                }
            }
            chng = -sum / s[j];
            if (ordinal) {
                chng = MAX(-bcoef[j], chng);
            }
            meps = MAX(meps, fabs(chng));
            bcoef[j] += chng;
            snew -= s[j] * chng * chng;
            for (int i = 0; i < m; i++) {
                dhatvec[i] += chng * basis[i][j];
                r[i] += chng * basis[i][j];
            }
        }
        // end of CCD cycle
        snew /= sum2;
        if (dverbose) {
            printf("CCD itel %3d sold %15.10f snew %15.10f meps %15.10f\n",
                   itel, sold, snew, meps);
        }
        if ((itel == ditmax) || (meps < deps)) {
            break;
        }
        sold = snew;
        itel = itel + 1;
    }  // end of repeat loop
    free(s);
    free(r);
    return;
}

/*
int  main(void) {
    int  n = 10, m = 3, itmax = 1000, eps = 10, width = 15, precision = 10;
    double d[10] = {0.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    double b[3] = {1.0, 2.0, 3.0};
    double dhat[10] = {0};
    double delta[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    double(*x)[n][m] = malloc(sizeof(*x));
    bool verbose = true, nonnegative = true, ordinal = false;
    double maxd = -INFINITY, mind = INFINITY;
    for (int  i = 0; i < n; i++) {
        maxd = MAX(maxd, delta[i]);
        mind = MIN(mind, delta[i]);
    }
    for (int  i = 0; i < n; i++) {
        delta[i] = (delta[i] - mind) / (maxd - mind);
    }
    (void)smacofBernsteinBase(n, m, delta, ordinal, x);
    for (int  i = 0; i < n; i++) {
        double sum = 0.0;
        for (int  j = 0; j < m; j++) {
            sum += (*x)[i][j] * b[j];
        }
        dhat[i] = sum;
    }
    (void)smacofUnweightedCCD(n, m, d, b, dhat, x, itmax, eps, verbose,
                              nonnegative);
    for (int  i = 0; i < n; i++) {
        printf("%15.10f ", dhat[i]);
    }
    printf("\n");
    free(x);
    return EXIT_SUCCESS;
}
*/
