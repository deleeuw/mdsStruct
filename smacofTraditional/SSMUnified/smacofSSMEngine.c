#include "smacof.h"

void smacofSSMEngine(const int n, const int p, double **delta, double **w,
                     double **xold, double **xnew, double **dmat, double **dhat,
                     double **bmat, const int init, const int itmax,
                     const int ieps1, const int ieps2, const bool verbose,
                     const bool relax, const bool weights, char *iterstring) {
    int itel = 1;
    double sold = 0.0, snew = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    (void)smacofCopyAnyMatrix(n, n, delta, dhat);
    (void)smacofNormDelta(n, dhat);
    (void)smacofInitial(n, p, init, dhat, xold);
    (void)smacofDistance(n, p, xold, dmat);
    (void)smacofScale(n, p, dhat, dmat, xold);
    sold = smacofStress(n, weights, dhat, w, dmat);
    // the fun starts here
    while (true) {
        (void)smacofGuttmanTransform(n, p, weights, delta, dmat, w, vinv, xold,
                                     xnew);
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
        }
        (void)smacofDistance(n, p, xnew, dmat);
        snew = smacofStress(n, weights, dhat, w, dmat);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate %12.10f\n",
                itel, sold, snew, sold - snew, chnew, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (chnew < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        chold = chnew;
        (void)smacofCopyAnyMatrix(n, p, xnew, xold);
    }
    sprintf(iterstring,
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    return;
}
