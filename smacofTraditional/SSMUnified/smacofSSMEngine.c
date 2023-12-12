#include "smacof.h"

void smacofSSMEngine(const int n, const int p, double **delta, double **w,
                     double **xold, double **xnew, double **dmat, double **dhat,
                     double **vmat, double **vinv, const int init,
                     const int itmax, const int ieps1, const int ieps2,
                     const bool verbose, const bool relax, const bool weights,
                     char *iterstring) {
    int itel = 1;
    double sold = 0.0, snew = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    if (weights) {
        (void)smacofNormWeights(n, w);
        (void)smacofMakeVMatrix(n, w, vmat);
        (void)smacofInverseVMatrix(n, vmat, vinv);
    }
    (void)smacofCopyAnyMatrix(n, n, delta, dhat);
    (void)smacofNormDelta(n, weights, dhat, w);
    (void)smacofInitial(n, p, init, weights, dhat, w, xold);
    (void)smacofDistance(n, p, xold, dmat);
    (void)smacofScale(n, p, weights, dhat, w, dmat, xold);
    sold = smacofStress(n, weights, dhat, w, dmat);
    //(void)smacofPrintSymmetricMatrix(stdout, n, 15, 10, w);
    //(void)smacofPrintSymmetricMatrix(stdout, n, 15, 10, vmat);
    //(void)smacofPrintSymmetricMatrix(stdout, n, 15, 10, vinv);
    //(void)smacofPrintSymmetricMatrix(stdout, n, 15, 10, dhat);
    //(void)smacofPrintSymmetricMatrix(stdout, n, 15, 10, dmat);
    //(void)smacofPrintAnyMatrix(stdout, n, p, 15, 10, xold);
    // the fun starts here
    while (true) {
        (void)smacofGuttmanTransform(n, p, weights, dhat, dmat, w, vinv, xold,
                                     xnew);
        chnew = smacofEtaSquareDifference(n, p, weights, vmat, xold, xnew);
        rate = sqrt(chnew / chold);
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
        }
        (void)smacofDistance(n, p, xnew, dmat);
        snew = smacofStress(n, weights, dhat, w, dmat);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f etad "
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
