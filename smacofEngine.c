#include "smacof.h"

// to be called from R

void smacofEngine(double *delta, double *weights, double *xold, double *xnew,
                  double *dnew, double *psnew, const int *ii, const int *jj,
                  const int *pm, const int *pn, const int *pp,
                  const int *pitmax, const int *peps1, const int *peps2,
                  const bool *verbose) {
    int m = *pm, n = *pn, p = *pp, nn = n * (n - 1) / 2, np = n * p, itel = 1;
    int width = 10, precision = 6;
    double sold = 0.0, snew = *psnew, change = 0.0, itmax = *pitmax;
    double eps1 = pow(10, -*peps1), eps2 = pow(10, -*peps2);
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vinv = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    (void)smacofNormWeights(weights, pm);
    if (DEBUG) {
        printf("normalized weights\n");
        (void)smacofPrintSymmetricHollowMatrix(weights, pn, &width, &precision);
    }
    (void)smacofNormDelta(delta, weights, pm);
    if (DEBUG) {
        printf("normalized delta\n");
        (void)smacofPrintSymmetricHollowMatrix(delta, pn, &width, &precision);
    }
    (void)smacofMPInverseSDCMatrix(weights, vinv, pn);
    if (DEBUG) {
        printf("MP Inverse V\n");
        (void)smacofPrintSDCMatrix(vinv, pn, &width, &precision);
    }
    (void)smacofInitial(delta, xold, pn, pp);
    (void)smacofCenter(xold, pn, pp);
    if (DEBUG) {
        printf("xold\n");
        (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
    }
    (void)smacofDist(xold, dold, ii, jj, pm, pn, pp);
    if (DEBUG) {
        printf("dold\n");
        (void)smacofPrintSymmetricHollowMatrix(dold, pn, &width, &precision);
    }
    (void)smacofScale(delta, weights, dold, xold, pm, pn, pp);
    if (DEBUG) {
        printf("scaled xold\n");
        (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
    }
    if (DEBUG) {
        printf("scaled dold\n");
        (void)smacofPrintSymmetricHollowMatrix(dold, pn, &width, &precision);
    }
    (void)smacofLoss(delta, weights, dold, pm, &sold);
    if (DEBUG) {
        printf("sold %10.6f\n\n", sold);
    }
    while (true) {
        (void)smacofGuttman(delta, weights, vinv, dold, xold, xnew, pn, pp);
        if (DEBUG) {
            printf("Guttman Transform\n");
            (void)smacofPrintAnyMatrix(xnew, pn, pp, &width, &precision);
        }
        (void)smacofDist(xnew, dnew, ii, jj, pm, pn, pp);
        (void)smacofLoss(delta, weights, dnew, pm, &snew);
        (void)smacofMaxDifference(xold, xnew, &change, pn, pp);
        if (*verbose) {
            printf("itel %4d sold %15.10f sdif %15.10f chng %15.10f\n", itel,
                   sold, sold - snew, change);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (change < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        (void)memcpy(xold, xnew, (size_t)(np * sizeof(double)));
        (void)memcpy(dold, dnew, (size_t)(m * sizeof(double)));
    }
    free(dold);
    free(vinv);
    return;
}

int main() {
    double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double weights[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double xold[8] = {1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0};
    double xnew[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double dnew[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double snew = 0.0;
    int ii[6] = {2, 3, 4, 3, 4, 4};
    int jj[6] = {1, 1, 1, 2, 2, 3};
    int m = 6, n = 4, p = 2, itmax = 1000, peps1 = 10, peps2 = 6;
    bool verbose = true;
    (void)smacofEngine(delta, weights, xold, xnew, dnew, &snew, ii, jj, &m, &n,
                       &p, &itmax, &peps1, &peps2, &verbose);
}
