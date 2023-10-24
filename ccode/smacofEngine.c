#include "smacof.h"

// to be called from R

void smacofEngine(double *delta, double *weights, double *xold, double *xnew,
                  double *dnew, double *psnew, const int *ii, const int *jj,
                  const int *pm, const int *pn, const int *pp, int *pitel,
                  const int *pitmax, const int *peps1, const int *peps2,
                  const bool *pverbose) {
    int m = *pm, itel = *pitel, itmax = *pitmax,
        itmax_j = 100;
    bool verbose = *pverbose, verbose_j = false;
    double sold = 0.0, snew = *psnew, cchange = 0.0, dchange = 0.0;
    double eps1 = pow(10, -*peps1), eps2 = pow(10, -*peps2),
           eps_j = pow(10, -15);
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *v = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vinv = (double *)calloc((size_t)m, (size_t)sizeof(double));
    (void)smacofNormWeights(weights, pm);
    (void)smacofNormDelta(delta, weights, pm);
    (void)smacofMakeVfromW(weights, v, pn);
    (void)smacofMPInverseSDCMatrix(weights, vinv, pn);
    (void)smacofTorgerson(delta, xold, pn, pp, &itmax_j, &eps_j, &verbose_j);
    (void)smacofCenter(xold, pn, pp);
    (void)smacofDist(xold, dold, ii, jj, pm, pn, pp);
    (void)smacofScale(delta, weights, dold, xold, pm, pn, pp);
    (void)smacofLoss(delta, weights, dold, pm, &sold);
    while (true) {
        (void)smacofGuttman(delta, weights, vinv, dold, xold, xnew, pn, pp);
        (void)smacofDist(xnew, dnew, ii, jj, pm, pn, pp);
        (void)smacofLoss(delta, weights, dnew, pm, &snew);
        (void)smacofMaxConfDifference(xold, xnew, &cchange, pn, pp);
        (void)smacofMaxDistDifference(dold, dnew, &dchange, pm);
        if (verbose) {
            printf(
                "itel %4d sold %15.10f sdif %15.10f cchange %15.10f dchange "
                "%15.10f\n",
                itel, sold, sold - snew, cchange, dchange);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (cchange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        (void)COPY(&xold, &xnew);
        (void)COPY(&dold, &dnew);
    }
    *psnew = snew;
    *pitel = itel;
    free(vinv);
    free(v);
    return;
}

int main() {
    double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double weights[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double xold[8] = {1, 1, -2, -2, 1, 1, -1, -1};
    double xnew[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double dnew[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double snew = 0.0;
    int ii[6] = {2, 3, 4, 3, 4, 4};
    int jj[6] = {1, 1, 1, 2, 2, 3};
    int m = 6, n = 4, p = 2, itel = 1, itmax = 100, peps1 = 15, peps2 = 10;
    bool verbose = true;
    (void)smacofEngine(delta, weights, xold, xnew, dnew, &snew, ii, jj, &m, &n,
                       &p, &itel, &itmax, &peps1, &peps2, &verbose);
}
