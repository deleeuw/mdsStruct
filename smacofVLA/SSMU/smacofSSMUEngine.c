#include "../Include/smacof.h"

// to be called from R

double delta[6] = {1.0, 1.4142, 1.0, 1.0, 1.4142, 1.0};
double xini[8] = {0.0, -2.0, -2.0, 4.0, 1.0, -1.0, -2.0, 2.0};
double xnew[8] = {0};
double dnew[6] = {0};
double bnew[10] = {0};
double snew = 0.0;
unsigned init = 4, n = 4, p = 2, itel = 1, itmax = 100;
unsigned ieps1 = 15, ieps2 = 10;
bool verbose = true, relax = false;

int main(void) {
    (void)smacofSSMUEngine(n, p, delta, xini, xnew, dnew, bnew, init, itmax,
                           ieps1, ieps2, verbose, relax, &itel, &snew);
    return EXIT_SUCCESS;
}

void smacofSSMUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const unsigned init, const unsigned itmax,
                      const unsigned ieps1, const unsigned ieps2,
                      const bool verbose, const bool relax, unsigned *pitel,
                      double *psnew) {
    unsigned itel = *pitel;
    unsigned width = 15, precision = 10;
    double sold = 0.0, snew = *psnew, ddiff = 0.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    double chnew = 1.0, chold = 1.0, rate = 1.0;
    double *pchnew = &chnew, *prate = &rate;
    /*
     * Convert delta from an R vector to a C matrix
     */
    double(*cdelta)[n][n] = malloc(sizeof *cdelta);
    assert(!(cdelta == NULL));
    (void)smacofFromSymmetricHollowRtoC(n, delta, cdelta);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdelta);
    }
    /*
     * Convert xini from an R vector to a C matrix
     */
    double(*cxini)[n][p] = malloc(sizeof *cxini);
    assert(!(cxini == NULL));
    (void)smacofFromAnyRtoC(n, p, xini, cxini);
    if (DEBUG) {
        printf("xini\n\n");
        (void)smacofPrintAnyMatrix(n, p, width, precision, cxini);
    }
    /*
     * Convert xnew from an R vector to a C matrix
     */
    double(*cxnew)[n][p] = malloc(sizeof *cxnew);
    assert(!(cxnew == NULL));
    double(*cdnew)[n][n] = malloc(sizeof *cdnew);
    assert(!(cdnew == NULL));
    double(*cbnew)[n][n] = malloc(sizeof *cbnew);
    assert(!(cbnew == NULL));
    double(*cxold)[n][p] = malloc(sizeof *cxold);
    assert(!(cxold == NULL));
    double(*cdold)[n][n] = malloc(sizeof *cdold);
    assert(!(cdold == NULL));
    double(*cdini)[n][n] = malloc(sizeof *cdini);
    assert(!(cdini == NULL));
    double(*cbold)[n][n] = malloc(sizeof *cbold);
    assert(!(cbold == NULL));
    (void)smacofUnweightedNormDelta(n, cdelta);
    if (DEBUG) {
        printf("delta normalized\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdelta);
    }
    (void)smacofUnweightedInitial(n, p, init, cdelta, cxini);
    (void)smacofDistance(n, p, cxini, cdini);
    (void)memcpy(cxold, cxini, sizeof *cxold);
    (void)memcpy(cdold, cdini, sizeof *cdold);
    if (DEBUG) {
        printf("xold = xini\n\n");
        (void)smacofPrintAnyMatrix(n, p, width, precision, cxold);
    }
    if (DEBUG) {
        printf("dold = dini\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdold);
    }
    (void)smacofUnweightedMakeBMatrix(n, cdelta, cdold, cbold);
    if (DEBUG) {
        printf("bold\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cbold);
    }
    sold = smacofUnweightedMakeStress(n, cdelta, cdold);
    if (DEBUG) {
        printf("sold %+*.*f\n\n", width, precision, sold);
    }
    /*
     * This is where the fun begins
     */
    while (true) {
        (void)smacofUnweightedGuttman(n, p, cbold, cxold, cxnew);
        if (DEBUG) {
            printf("xnew\n\n");
            (void)smacofPrintAnyMatrix(n, p, width, precision, cxnew);
        }
        /*
        (void)smacofRelax(n, p, cxold, cxnew, chold, &chnew, &rate, itel,
        relax); if (DEBUG) { printf("xrelaxed\n\n");
            (void)smacofPrintAnyMatrix(n, p, width, precision, cxnew);
        }
        */
        (void)smacofDistance(n, p, cxnew, cdnew);
        if (DEBUG) {
            printf("dnew\n\n");
            (void)smacofPrintAnyMatrix(n, n, width, precision, cdnew);
        }
        ddiff = smacofMaxDistanceDifference(n, cdold, cdnew);
        (void)smacofUnweightedMakeBMatrix(n, cdelta, cdnew, cbnew);
        if (DEBUG) {
            printf("** dnew\n\n");
            (void)smacofPrintAnyMatrix(n, n, width, precision, cdnew);
        }
        if (DEBUG) {
            printf("bnew\n\n");
            (void)smacofPrintAnyMatrix(n, n, width, precision, cbnew);
        }   
        snew = smacofUnweightedMakeStress(n, cdelta, cdnew);
        if (DEBUG) {
            printf("snew %15.10f\n\n", snew);
        }
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate %12.10f ddiff %12.10f\n",
                itel, sold, snew, sold - snew, *pchnew, *prate, ddiff);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (ddiff < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        chold = *pchnew;
        (void)memcpy(cxold, cxnew, sizeof *cxold);
        (void)memcpy(cdold, cdnew, sizeof *cdold);
        (void)memcpy(cbold, cbnew, sizeof *cbold);
    }
    /*
    * We hand everything back to R
    * delta, dnew, xnew, bnew
    * how about snew, itel
    *
    * Consider passing a block of storage from R large enough
    * to contain the matrices and pass the block back to R
    * appropriately filled
    */
    free(cdelta);
    free(cdold);
    free(cdnew);
    free(cdini);
    free(cxold);
    free(cxnew);
    free(cxini);
    free(cbold);
    free(cbnew);
    return;
}
