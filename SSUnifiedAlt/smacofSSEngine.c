#include "smacof.h"

void smacofSSEngine(const int n, const int p, double **delta, double **w,
                    double **vmat, double **vinv, double **xold, double **xnew,
                    double **dmat, double **dhat, double **basis,
                    const int haveinit, const int itmax, const int ieps1,
                    const int ieps2, const bool verbose, const bool relax,
                    const int ditmax, const int ieps3, const bool dverbose,
                    const bool weights, const int degree, const int ordinal,
                    const int transform, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    double *dmatvec = NULL, *dhatvec = NULL, *wvec = NULL, *bcoef = NULL;
    dmatvec = smacofMakeVector(m);
    dhatvec = smacofMakeVector(m);
    bcoef = smacofMakeVector(degree);
    if (weights) {
        wvec = smacofMakeVector(m);
        (void)smacofSymmetricCtoR(n, w, wvec);
    }
    (void)smacofCopyAnyMatrix(n, n, delta, dhat);
    (void)smacofInitial(n, p, haveinit, weights, dhat, w, xold);
    (void)smacofDistance(n, p, xold, dmat);
    sold = smacofStress(n, p, weights, dhat, dmat, xold, w);
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
        smid = smacofStress(n, p, weights, dhat, dmat, xnew, w);
        if (transform) {
            (void)smacofSymmetricCtoR(n, dmat, dmatvec);
            (void)smacofSymmetricCtoR(n, dhat, dhatvec);
            (void)smacofCCD(n, degree, dmatvec, wvec, bcoef, dhatvec, basis,
                            ditmax, ieps3, dverbose, weights, ordinal);
            (void)smacofSymmetricRtoC(n, dhatvec, dhat);
        }
        snew = smacofStress(n, p, weights, dhat, dmat, xnew, w);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif %+12.10f "
                "etad %+12.10f rate %12.10f\n",
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
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    (void)smacofFreeVector(dhatvec);
    (void)smacofFreeVector(dmatvec);
    if (weights) {
        (void)smacofFreeVector(wvec);
        (void)smacofFreeMatrix(n, vmat);
        (void)smacofFreeMatrix(n, vinv);
    }
    return;
}
