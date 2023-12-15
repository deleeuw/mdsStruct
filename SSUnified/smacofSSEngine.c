#include "smacof.h"

void smacofSSMEngine(const int n, const int p, double **delta, double **w,
                     double **xold, double **xnew, double **dmat, double **dhat,
                     double **basis, const int init, const int itmax,
                     const int ieps1, const int ieps2, const bool verbose,
                     const bool relax, const bool weights, const int transform,
                     const int degree, const int ordinal, char *iterstring) {
    int itel = 1, m = n * (n - 1) / 2;
    double sold = 0.0, snew = 0.0, smid = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 0.0, chold = INFINITY, rate = 0.0;
    double **vmat = NULL, **vinv = NULL;
    double *dmatvec = NULL, *dhatvec = NULL, *wvec = NULL, *bcoef = NULL;
    if (weights) {
        vmat = smacofMakeSymmetricMatrix(n);
        vinv = smacofMakeSymmetricMatrix(n);
        (void)smacofNormWeights(n, w);
        (void)smacofMakeVMatrix(n, w, vmat);
        (void)smacofInverseVMatrix(n, vmat, vinv);
    } else {
        vmat = smacofMakeSymmetricMatrix(0);
        vinv = smacofMakeSymmetricMatrix(0);
    }
    if ((transform == POLYNOMIAL) || (transform == SPLINICAL)) {
        dmatvec = smacofMakeVector(m);
        dhatvec = smacofMakeVector(m);
        bcoef = smacofMakeVector(degree);
        if (weights) {
            wvec = smacofMakeVector(m);
            (void)smacofSymmetricCtoR(n, w, wvec);
        }
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
        smid = smacofStress(n, weights, dhat, w, dmat);
        if (transform == INTERVAL) {
            (void)smacofInterval(n, delta, dmat, dhat);
            (void)smacofNormDelta(n, weights, dhat, w);
        }
        if (transform == POLYNOMIAL) {
            (void)smacofSymmetricCtoR(n, dmat, dmatvec);
            (void)smacofCCD(n, degree, dmatvec, wvec, bcoef, dhatvec, basis, 3,
                            10, false, weights, ordinal);
            (void)smacofSymmetricRtoC(n, dhatvec, dhat);
            (void)smacofNormDelta(n, weights, dhat, w);
        }
        snew = smacofStress(n, weights, dhat, w, dmat);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f smid %12.10f snew %12.10f sdif %+12.10f "
                "etad "
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
            "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd %+12.10f "
            "rate %12.10f\n",
            itel, sold, snew, sold - snew, chnew, rate);
    if ((transform == POLYNOMIAL) || (transform == SPLINICAL)) {
        (void)smacofFreeVector(dhatvec);
        (void)smacofFreeVector(dmatvec);
        if (weights) {
            (void)smacofFreeVector(wvec);
        }
    }
    if (weights) {
        (void)smacofFreeMatrix(n, vmat);
        (void)smacofFreeMatrix(n, vinv);
    }
    return;
}
