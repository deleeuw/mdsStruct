#include "../Unweighted/smacofUnweighted.h"

void smacofSSMUEngine(const int n, const int p, double **delta, double **xold,
                      double **xnew, double **dmat, double **bmat,
                      const int init, const int itmax, const int ieps1,
                      const int ieps2, const bool verbose, const bool relax,
                      int *pitel, double *psnew) {
    int itel = *pitel;
    double sold = 0.0, snew = *psnew;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    int width = 15, precision = 10;
    if (DEBUG) {
        printf("delta\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, delta);
    }
    (void)smacofUnweightedNormDelta(n, delta);
    if (DEBUG) {
        printf("delta normalized\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, delta);
    }
    (void)smacofUnweightedInitial(n, p, init, delta, xold);
    if (DEBUG) {
        printf("xini\n");
        (void)smacofPrintAnyMatrix(n, p, width, precision, xold);
    }
    (void)smacofDistance(n, p, xold, dmat);
    if (DEBUG) {
        printf("dmat\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, dmat);
    }
    (void)smacofUnweightedScale(n, p, delta, dmat, xold);
    if (DEBUG) {
        printf("xold\n");
        (void)smacofPrintAnyMatrix(n, p, width, precision, xold);
    }
    if (DEBUG) {
        printf("dmat\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, dmat);
    }
    (void)smacofUnweightedMakeBMatrix(n, delta, dmat, bmat);
    if (DEBUG) {
        printf("bmat\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, bmat);
    }
    sold = smacofUnweightedMakeStress(n, delta, dmat);
    if (DEBUG) {
        printf("** sold %15.10f\n\n", sold);
    }
    // the fun starts here
    while (true) {
        (void)smacofUnweightedGuttman(n, p, bmat, xold, xnew);
        if (DEBUG) {
            printf("iteration %3d xnew\n", itel);
            (void)smacofPrintAnyMatrix(n, p, width, precision, xnew);
        }
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
        }
        (void)smacofDistance(n, p, xnew, dmat);
        if (DEBUG) {
            printf("iteration %3d dmat\n", itel);
            (void)smacofPrintAnyMatrix(n, n, width, precision, dmat);
        }
        (void)smacofUnweightedMakeBMatrix(n, delta, dmat, bmat);
        if (DEBUG) {
            printf("iteration %3d bmat\n", itel);
            (void)smacofPrintAnyMatrix(n, n, width, precision, bmat);
        }
        snew = smacofUnweightedMakeStress(n, delta, dmat);
        if (DEBUG) {
            printf("iteration %3d ** snew %15.10f\n\n", itel, snew);
        }
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
    return;
}
