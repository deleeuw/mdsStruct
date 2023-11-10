#include "../smacofInclude/smacof.h"

// to be called from R

void smacofSSMUEngine(double *delta, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax) {
    int n = *pn, p = *pp, np = p * n, m = n * (n - 1) / 2, nn = n * (n + 1) / 2,
        itel = *pitel, itmax = *pitmax;
    int width = 15, precision = 10;
    bool verbose = *pverbose;
    double sold = 0.0, snew = *psnew, pchange = 1.0, echange = 1.0, rate = 1.0;
    double eps1 = pow(10.0, -(double)*peps1), eps2 = pow(10.0, -(double)*peps2);
    double *xold = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *dini = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *bold = (double *)calloc((size_t)nn, (size_t)sizeof(double));
    (void)smacofUnweightedNormDelta(delta, pn);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrintSHMatrix(delta, &n, &width, &precision);
    }
    (void)smacofUnweightedInitial(delta, xini, pinit, pn, pp);
    if (DEBUG) {
        printf("xini\n\n");
        (void)smacofPrintAnyMatrix(xini, &n, &p, &width, &precision);
    }
    (void)smacofDistance(xini, dini, pn, pp);
    if (DEBUG) {
        printf("dini\n\n");
        (void)smacofPrintSHMatrix(dini, &n, &width, &precision);
    }
    (void)memcpy(xold, xini, (size_t)np * sizeof(double));
    (void)memcpy(dold, dini, (size_t)m * sizeof(double));
    (void)smacofUnweightedMakeBMatrix(delta, dold, bold, pn);
    if (DEBUG) {
        printf("bold\n\n");
        (void)smacofPrintSymmetricMatrix(bold, &n, &width, &precision);
    }
    (void)smacofUnweightedMakeStress(delta, dold, &sold, pn);
    if (DEBUG) {
        printf("sold %15.10f\n\n", sold);
    }
    while (true) {
        (void)smacofUnweightedGuttman(bold, xold, xnew, pn, pp);
        (void)smacofRMSDifference(xold, xnew, pn, pp, &echange);
        /// this is confusing, skip when relax == FALSE
        (void)smacofRelax(xold, xnew, &echange, &pchange, &np, &itel, prelax,
                          &rate);
        (void)smacofDistance(xnew, dnew, pn, pp);
        (void)smacofUnweightedMakeBMatrix(delta, dnew, bnew, pn);
        (void)smacofUnweightedMakeStress(delta, dnew, &snew, pn);
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
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        (void)memcpy(bold, bnew, (size_t)nn * sizeof(double));
    }
    *psnew = snew;
    *pitel = itel;
    free(dold);
    free(dini);
    free(xold);
    free(bold);
    return;
}
