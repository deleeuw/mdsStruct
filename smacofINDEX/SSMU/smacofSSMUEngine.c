#include "../Unweighted/smacofUnweighted.h"

// to be called from R

void smacofSSMUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const unsigned init, const unsigned itmax, const unsigned peps1,
                      const unsigned peps2, const bool verbose, const bool relax,
                      unsigned *pitel, double *psnew) {
    unsigned np = p * n, m = n * (n - 1) / 2, nn = n * (n + 1) / 2, itel = *pitel;
    unsigned width = 15, precision = 10;
    double sold = 0.0, snew = *psnew, pchange = 1.0, echange = 1.0, rate = 1.0;
    double eps1 = pow(10.0, -(double)peps1), eps2 = pow(10.0, -(double)peps2);
    double *xold = (double *)calloc((unsigned)np, (unsigned)sizeof(double));
    double *dold = (double *)calloc((unsigned)m, (unsigned)sizeof(double));
    double *dini = (double *)calloc((unsigned)m, (unsigned)sizeof(double));
    double *bold = (double *)calloc((unsigned)nn, (unsigned)sizeof(double));
    (void)smacofUnweightedNormDelta(n, delta);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrintSHMatrix(n, delta, width, precision);
    }
    (void)smacofUnweightedInitial(n, p, delta, xini, init);
    if (DEBUG) {
        printf("xini\n\n");
        (void)smacofPrintAnyMatrix(n, p, xini, width, precision);
    }
    (void)smacofDistance(n, p, xini, dini);
    if (DEBUG) {
        printf("dini\n\n");
        (void)smacofPrintSHMatrix(n, dini, width, precision);
    }
    (void)memcpy(xold, xini, (unsigned)np * sizeof(double));
    (void)memcpy(dold, dini, (unsigned)m * sizeof(double));
    (void)smacofUnweightedMakeBMatrix(n, delta, dold, bold);
    if (DEBUG) {
        printf("bold\n\n");
        (void)smacofPrintSymmetricMatrix(n, bold, width, precision);
    }
sold = smacofUnweightedMakeStress(n, delta, dold);
    if (DEBUG) {
        printf("sold %15.10f\n\n", sold);
    }
    while (true) {
        (void)smacofUnweightedGuttman(n, p, bold, xold, xnew);
        echange = smacofRMSDifference(n, p, xold, xnew);
        /// this is confusing, skip when relax == FALSE
        // (void)smacofRelax(xold, xnew, &echange, pchange, &np, &itel, relax, &rate);
        (void)smacofDistance(n, p, xnew, dnew);
        (void)smacofUnweightedMakeBMatrix(n, delta, dnew, bnew);
        snew = smacofUnweightedMakeStress(n, delta, dnew);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate "
                "%12.10f\n",
                itel, sold, snew, sold - snew, echange, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (echange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        pchange = echange;
        (void)memcpy(xold, xnew, (unsigned)np * sizeof(double));
        (void)memcpy(dold, dnew, (unsigned)m * sizeof(double));
        (void)memcpy(bold, bnew, (unsigned)nn * sizeof(double));
    }
    *psnew = snew;
    *pitel = itel;
    free(dold);
    free(dini);
    free(xold);
    free(bold);
    return;
}
