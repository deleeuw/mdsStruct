#include "../Unweighted/smacofUnweighted.h"

void smacofSSPUEngine(const int n, const int p, double **delta, 
                      double **basis, double *coef, double **xold,
                      double **xnew, double **dmat, double **dhat,
                      double **bmat, const int init, const int itmax,
                      const int ieps1, const int ieps2, const bool verbose,
                      const bool relax, const int degree, const int ordinal, 
                      char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2, citmax = 100, ceps = 10;
    bool cverbose = false;
    double sold = 0.0, smid = 0.0, snew = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 1.0, chold = INFINITY, rate = 1.0;
    double *dvec = smacofMakeAnyVector(m);
    double *vhat = smacofMakeAnyVector(m);
    (void)smacofCopyAnyMatrix(n, n, delta, dhat);
    (void)smacofUnweightedNormDelta(n, dhat);
    (void)smacofUnweightedInitial(n, p, init, dhat, xold);
    (void)smacofDistance(n, p, xold, dmat);
    sold = smacofUnweightedMakeStress(n, dhat, dmat);
    /*
     * the fun starts here
     */
    while (true) {
        (void)smacofUnweightedMakeBMatrix(n, dhat, dmat, bmat);
        (void)smacofUnweightedGuttman(n, p, bmat, xold, xnew);
        chnew = smacofRMSDifference(n, p, xold, xnew);
        rate = chnew / chold;
        if (relax) {
            (void)smacofRelax(n, p, rate, xold, xnew);
        }
        (void)smacofDistance(n, p, xnew, dmat);
        smid = smacofUnweightedMakeStress(n, dhat, dmat);
        (void)smacofFromSymmetricCtoR(n, dmat, dvec);
        (void)smacofUnweightedCCD(n, degree, dvec, coef, vhat, basis, citmax, ceps,
                                 cverbose, ordinal);
        (void)smacofFromSymmetricHollowRtoC(n, vhat, dhat);
        (void)smacofUnweightedNormDelta(n, dhat);
        snew = smacofUnweightedMakeStress(n, dhat, dmat);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif %+12.10f "
                "rmsd "
                "%+12.10f rate %12.10f\n",
                itel, sold, smid, snew, sold - snew, chnew, rate);
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
            "itel %3d sold %12.10f smid %12.10f, snew %12.10f sdif %+12.10f "
            "rmsd %+12.10f rate %12.10f\n",
            itel, sold, smid, snew, sold - snew, chnew, rate);
    smacofFreeAnyVector(dvec);
    smacofFreeAnyVector(vhat);
    return;
}
