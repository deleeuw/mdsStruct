#include "smacof.h"

// to be called from R

void smacofEngine(double *delta, double *weights, double *xini, double *xnew,
                  double *dnew, double *bmat, double *psnew, const int *ii,
                  const int *jj, const int *pm, const int *pn, const int *pp,
                  int *pitel, const int *pitmax, const int *peps1,
                  const int *peps2, const int *pinit, const bool *pverbose) {
    int m = *pm, n = *pn, p = *pp, np = p * n, itel = *pitel, itmax = *pitmax,
        itmax_j = 100, init = *pinit;
    int width = 15, precision = 10;
    bool verbose = *pverbose, verbose_j = false, verbose_e = false;
    double sold = 0.0, snew = *psnew, cchange = 0.0, dchange = 0.0;
    double eps1 = pow(10, -*peps1), eps2 = pow(10, -*peps2),
           eps_j = pow(10, -15), eps_e = pow(10, -15);
    double *xold = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vmat = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vinv = (double *)calloc((size_t)m, (size_t)sizeof(double));
    (void)smacofNormWeights(weights, pm);
    (void)smacofNormDelta(delta, weights, pm);
    (void)smacofMakeVmat(weights, vmat, pn);
    (void)smacofMPInverseSDCMatrix(weights, vinv, pn);
    (void)smacofMakeVmat(weights, vmat, pn);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrintSHMatrix(delta, pn, &width, &precision);
        printf("weights\n\n");
        (void)smacofPrintSHMatrix(weights, pn, &width, &precision);
        printf("vmat\n\n");
        (void)smacofPrintSDCMatrix(vmat, pn, &width, &precision);
    }
    //(void)smacofInitial();
    switch (init) {
        case 1:
            (void)smacofTorgerson(delta, xini, pn, pp, &itmax_j, &eps_j,
                                  &verbose_j);
            break;
        case 2:
            //(void)smacofElegant(delta, xini, pn, pp, &itmax_e, &eps_e,
            //&verbose_e);
            break;
        case 3:
            //(void)smacofInitSDCmat(delta, weights, pn, &itmax_j, &eps_j,
            //                     &verbose_j);
            break;
        case 4:
            for (int i = 1; i <= n; i++) {
                for (int s = 1; s <= p; s++) {
                    xini[MINDEX(i, s, n)] = drand48();
                }
            }
            break;
    }
    if (DEBUG) {
        printf("xini in %p\n\n", xini);
        (void)smacofPrintAnyMatrix(xini, pn, pp, &width, &precision);
    }
    (void)memcpy(xold, xini, (size_t)np * sizeof(double));
    if (DEBUG) {
        printf("xold in %p\n\n", xold);
        (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
    }
    (void)smacofCenter(xold, pn, pp);
    (void)smacofDist(xold, dold, ii, jj, pm, pn, pp);
    (void)smacofScale(delta, weights, dold, xold, pm, pn, pp);
    (void)smacofLoss(delta, weights, dold, pm, &sold);
    while (true) {
        (void)smacofMakeBmat(delta, weights, dold, bmat, pm);
        if (DEBUG) {
            printf("+++++++++++++++++++++++++++++++\n\n");
            printf("beginning iteration %3d\n\n", itel);
            printf("bmat in %p\n\n", bmat);
            (void)smacofPrintSDCMatrix(bmat, pn, &width, &precision);
            printf("xold in %p\n\n", xold);
            (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
        }
        (void)smacofGuttman(vinv, bmat, xold, xnew, pn, pp);
        if (DEBUG) {
            printf("after guttman\n\n");
            printf("xold in %p\n\n", xold);
            (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
        }
        (void)smacofDist(xnew, dnew, ii, jj, pm, pn, pp);
        (void)smacofLoss(delta, weights, dnew, pm, &snew);
        (void)smacofMaxConfDifference(xold, xnew, &cchange, pn, pp);
        (void)smacofMaxDistDifference(dold, dnew, &dchange, pm);
        if (verbose) {
            printf(
                "itel %4d sold %15.10f sdif %15.10f cchange %15.10f "
                "dchange "
                "%15.10f\n",
                itel, sold, sold - snew, cchange, dchange);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (cchange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        if (DEBUG) {
            printf("after swap\n\n");
            printf("xold in %p\n\n", xold);
            (void)smacofPrintAnyMatrix(xold, pn, pp, &width, &precision);
            printf("dold in %p\n\n", dold);
            (void)smacofPrintSHMatrix(dold, pn, &width, &precision);
        }
    }
    *psnew = snew;
    *pitel = itel;
    free(vinv);
    free(vmat);
    free(dold);
    free(xold);
    return;
}

int main() {
    double delta[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double weights[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double bmat[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xini[8] = {1.0, 2.0, 3.0, 4.0, 1.0, 4.0, 9.0, 16.0};
    double xnew[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double dnew[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double snew = 0.0;
    int ii[6] = {2, 3, 4, 3, 4, 4};
    int jj[6] = {1, 1, 1, 2, 2, 3};
    int m = 6, n = 4, p = 2, itel = 1, itmax = 100, init = 4, peps1 = 15,
        peps2 = 10;
    bool verbose = true;
    (void)smacofEngine(delta, weights, xini, xnew, dnew, bmat, &snew, ii, jj,
                       &m, &n, &p, &itel, &itmax, &init, &peps1, &peps2,
                       &verbose);
}
