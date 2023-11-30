#include "../Unweighted/smacofUnweighted.h"

void smacofSSMUEngine(const unsigned n, const unsigned p, double **delta,
                      double **xini, double **xnew, double **dnew, double **bnew,
                      const unsigned init, const unsigned itmax,
                      const unsigned ieps1, const unsigned ieps2,
                      const bool verbose, const bool relax, unsigned *pitel,
                      double *psnew) {
    unsigned itel = *pitel;
    double sold = 0.0, snew = *psnew, ddiff = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 1.0, chold = INFINITY, rate = 1.0;
    (void)smacofUnweightedNormDelta(n, delta);
    (void)smacofUnweightedInitial(n, p, init, delta, xini);
    (void)smacofDistance(n, p, xini, dini);
    (void)memcpy(xold, xini, sizeof *xold);
    (void)memcpy(dold, dini, sizeof *dold);
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
        ddiff = smacofMaxDistanceDifference(n, dold, dnew);
        (void)smacofUnweightedMakeBMatrix(n, delta, dnew, bnew);
        snew = smacofUnweightedMakeStress(n, delta, dnew);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate %12.10f ddiff %12.10f\n",
                itel, sold, snew, sold - snew, chnew, rate, ddiff);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (ddiff < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        chold = chnew;
        (void)memcpy(xold, xnew, sizeof *xold);
        (void)memcpy(dold, dnew, sizeof *dold);
        (void)memcpy(bold, bnew, sizeof *bold);
    }
    return;
}
