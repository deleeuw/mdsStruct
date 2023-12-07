#include "../Weighted/smacofWeighted.h"

void smacofSSMWEngine(const int n, const int p, double **delta,
                      double **weights, double **xold, double **xnew,
                      double **dmat, double **bmat, double **vmat,
                      double **vinv, const int init, const int itmax,
                      const int ieps1, const int ieps2, const bool verbose,
                      const bool relax, char *iterstring) {
    int itel = 1;
    double sold = 0.0, snew = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    (void)smacofWeightedNormWeights(n, weights);
    (void)smacofWeightedNormDelta(n, delta, weights);
    (void)smacofWeightedMakeVMatrix(n, weights, vmat);
    (void)smacofWeightedInverseVMatrix(n, vmat, vinv);
    (void)smacofWeightedInitial(n, p, init, delta, weights, xold);
    (void)smacofDistance(n, p, xold, dmat);
    (void)smacofWeightedScale(n, p, delta, weights, dmat, xold);
    sold = smacofWeightedMakeStress(n, delta, weights, dmat);
    // the fun starts here
    while (true) {
        (void)smacofWeightedMakeBMatrix(n, delta, weights, dmat, bmat);
        (void)smacofWeightedGuttman(n, p, vinv, bmat, xold, xnew);
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
        }
        (void)smacofDistance(n, p, xnew, dmat);
        (void)smacofWeightedMakeBMatrix(n, delta, weights, dmat, bmat);
        snew = smacofWeightedMakeStress(n, delta, weights, dmat);
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
