#include "../Unweighted/smacofUnweighted.h"

void smacofSSMUEngine(const int n, const int p, double **delta, double **xini,
                      double **xnew, double **dnew, double **bnew,
                      const int init, const int itmax, const int ieps1,
                      const int ieps2, const bool verbose, const bool relax,
                      int *pitel, double *psnew) {
    int itel = *pitel;
    double sold = 0.0, snew = *psnew, ddiff = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    double **dini = smacofMakeAnyMatrix(n, n);
    double **dold = smacofMakeAnyMatrix(n, n);
    double **xold = smacofMakeAnyMatrix(n, p);
    double **bold = smacofMakeAnyMatrix(n, n);
    (void)smacofUnweightedNormDelta(n, delta);
    (void)smacofUnweightedInitial(n, p, init, delta, xini);
    (void)smacofDistance(n, p, xini, dini);
    (void)smacofCopyAnyMatrix(n, p, xini, xold);
    (void)smacofCopyAnyMatrix(n, n, dini, dold);
    (void)smacofUnweightedMakeBMatrix(n, delta, dold, bold);
    sold = smacofUnweightedMakeStress(n, delta, dold);
    while (true) {
        (void)smacofUnweightedGuttman(n, p, bold, xold, xnew);
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
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
        (void)smacofCopyAnyMatrix(n, p, xnew, xold);
        (void)smacofCopyAnyMatrix(n, p, dnew, dold);
        (void)smacofCopyAnyMatrix(n, p, bnew, bold);
    }
    return;
}
