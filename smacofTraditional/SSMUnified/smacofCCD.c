#include "smacof.h"

void smacofCCD(const int n, const int m, double *y, double *b, double *dhat,
               double **x, const int itmax, const int eps, const bool verbose,
               const bool nonnegative) {
    double *s = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *r = (double *)calloc((size_t)n, (size_t)sizeof(double));
    double deps = pow(10.0, -(double)eps);
    int itel = 1;
    double sold = 0.0;
    for (int j = 0; j < m; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += SQUARE(x[i][j]);
        }
        s[j] = sum;
    }
    for (int i = 0; i < n; i++) {
        r[i] = dhat[i] - y[i];
        sold += SQUARE(r[i]);
    }
    double snew = sold;
    while (true) {
        double meps = 0.0, chng = 0.0;
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += x[i][j] * r[i];
            }
            chng = -sum / s[j];
            if (nonnegative) {
                chng = MAX(-b[j], chng);
            }
            meps = MAX(meps, fabs(chng));
            b[j] += chng;
            snew -= s[j] * chng * chng;
            for (int i = 0; i < n; i++) {
                dhat[i] += chng * x[i][j];
                r[i] += chng * x[i][j];
            }
        }  // end of CCD cycle
        if (verbose) {
            printf("itel %3d sold %15.10f snew %15.10f meps %15.10f\n", itel,
                   sold, snew, meps);
        }
        if ((itel == itmax) || (meps < deps)) {
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
