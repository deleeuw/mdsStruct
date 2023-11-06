#include <smacof.h>

// to be called from R

void smacofSSMWEngine(double *delta, double *weights, const int *irow,
                      const int *icol, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax) {
    int n = *pn, p = *pp, np = p * n, m = *pm, nn = n * (n + 1) / 2,
        itel = *pitel, itmax = *pitmax;
    int width = 15, precision = 10;
    bool verbose = *pverbose, verbose_j = false, verbose_e = false;
    double sold = 0.0, snew = *psnew, cchange = 0.0, dchange = 0.0,
           pchange = 1.0, echange = 1.0, rate = 1.0;
    double eps1 = pow(10.0, -(double)*peps1), eps2 = pow(10.0, -(double)*peps2),
           eps_e = pow(10.0, -15.0);
    double *xold = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *dini = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *bold = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *vmat = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    double *vinv = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    (void)smacofNormWeights(weights, pm);
    (void)smacofNormDelta(delta, weights, pm);
    (void)smacofMakeVMatrix(weights, vmat, irow, icol, pn, pm);
    (void)smacofMPInverseSDCLMatrix(vmat, vinv, pn);
    (void)smacofInitial(delta, weights, irow, icol, xini, pinit, pn, pp, pm);
    (void)smacofDistance(xini, dini, irow, icol, pn, pp, pm);
    (void)memcpy(xold, xini, (size_t)np * sizeof(double));
    (void)memcpy(dold, dini, (size_t)m * sizeof(double));
    (void)smacofMakeBMatrix(delta, weights, dold, bold, irow, icol, pn, pm);
    (void)smacofStress(delta, weights, dold, pm, &sold);
    while (true) {
        (void)smacofGuttman(vinv, bold, xold, xnew, pn, pp);
        (void)smacofRMSDifference(xold, xnew, pn, pp, &echange);
        /// this is confusing, skip when relax == FALSE
        (void)smacofRelax(xold, xnew, &echange, &pchange, &np, &itel, prelax,
                          &rate);
        (void)smacofDistance(xnew, dnew, irow, icol, pn, pp, pm);
        (void)smacofMakeBMatrix(delta, weights, dnew, bnew, irow, icol, pn, pm);
        (void)smacofStress(delta, weights, dnew, pm, &snew);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f rate "
                "%12.10f\n",
                itel, sold, snew, sold - snew, echange, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (echange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        pchange = echange;
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        (void)memcpy(bold, bnew, (size_t)nn * sizeof(double));
    }
    *psnew = snew;
    *pitel = itel;
    free(vinv);
    free(vmat);
    free(dold);
    free(dini);
    free(xold);
    free(bold);
    return;
}
