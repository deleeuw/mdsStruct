#include "../Unweighted/smacofUnweighted.h"

// to be called from R

void smacofSSMUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const unsigned init, const unsigned itmax,
                      const unsigned ieps1, const unsigned ieps2,
                      const bool verbose, const bool relax, unsigned *pitel,
                      double *psnew) {
    unsigned np = p * n, m = n * (n - 1) / 2, nn = n * (n + 1) / 2,
             itel = *pitel;
    unsigned width = 15, precision = 10;
    double sold = 0.0, snew = *psnew, chold = INFINITY, chnew = 1.0, rate = 1.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double *xold = (double *)calloc((unsigned)np, (unsigned)sizeof(double));
    assert(xold != NULL);
    double *dold = (double *)calloc((unsigned)m, (unsigned)sizeof(double));
    assert(dold != NULL);
    double *dini = (double *)calloc((unsigned)m, (unsigned)sizeof(double));
    assert(dini != NULL);
    double *bold = (double *)calloc((unsigned)nn, (unsigned)sizeof(double));
        assert(bold != NULL);
    (void)smacofUnweightedNormDelta(n, delta);
    (void)smacofUnweightedInitial(n, p, delta, xini, init);
    (void)smacofDistance(n, p, xini, dini);
    (void)memcpy(xold, xini, (unsigned)np * sizeof(double));
    (void)memcpy(dold, dini, (unsigned)m * sizeof(double));
    (void)smacofUnweightedMakeBMatrix(n, delta, dold, bold);
    sold = smacofUnweightedMakeStress(n, delta, dold);
    while (true) {
        (void)smacofUnweightedGuttman(n, p, bold, xold, xnew);
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, xold, xnew, rate);
        }
        (void)smacofDistance(n, p, xnew, dnew);
        (void)smacofUnweightedMakeBMatrix(n, delta, dnew, bnew);
        snew = smacofUnweightedMakeStress(n, delta, dnew);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate "
                "%12.10f\n",
                itel, sold, snew, sold - snew, chnew, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (chnew < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        chold = chnew;
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
