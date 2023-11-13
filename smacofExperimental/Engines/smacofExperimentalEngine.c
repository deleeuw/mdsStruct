#include "../Include/smacof.h"

// to be called from R

double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
double xini[8];
double xnew[8];
double dnew[6];
double bnew[10];
double snew = 0.0;
size_t init = 3, n = 4, p = 2, itel = 1, itmax = 1;
int ieps1 = 15, ieps2 = 10;
bool verbose = true, relax = true;

int main(void) {
    (void)smacofSSMUEngine(n, p, delta, xini, xnew, dnew, bnew, init, itmax,
                           ieps1, ieps2, verbose, relax, &itel, &snew);
    return EXIT_SUCCESS;
}

void smacofSSMUEngine(const size_t n, const size_t p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const size_t init, const size_t itmax, const size_t ieps1,
                      const size_t ieps2, const bool verbose, const bool relax,
                      size_t *pitel, double *psnew) {
    size_t np = p * n, nn = SQUARE(n), itel = *pitel;
    int width = 15, precision = 10;
    double sold = 0.0, snew = *psnew, pchange = 1.0, echange = 1.0, rate = 1.0;
    double eps1 = pow(10.0, -(double)ieps1), eps2 = pow(10.0, -(double)ieps2);
    /*
     * Convert delta from an R vector to a C matrix
     */
    double(*cdelta)[n] = malloc((sizeof *cdelta) * n);
    assert(!(cdelta == NULL));
    (void)smacofFromSymmetricHollowRtoC(delta, n, cdelta);
    if (DEBUG) {
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdelta);
    }
    double(*cxini)[n] = malloc((sizeof *cxini) * p);
    assert(!(cxini == NULL));
    double(*cxold)[n] = malloc((sizeof *cxold) * p);
    assert(!(cxold == NULL));
    double(*cdold)[n] = malloc((sizeof *cdold) * n);
    assert(!(cdold == NULL));
    double(*cdini)[n] = malloc((sizeof *cdini) * n);
    assert(!(cdini == NULL));
    double(*cbold)[n] = malloc((sizeof *cbold) * n);
    assert(!(cbold == NULL));
    double(*cdnew)[n] = malloc((sizeof *cdnew) * n);
    assert(!(cdnew == NULL));
    double(*cbnew)[n] = malloc((sizeof *cbnew) * n);
    assert(!(cbnew == NULL));
    double(*cxnew)[n] = malloc((sizeof *cxnew) * p);
    assert(!(cxnew == NULL));
    (void)smacofUnweightedNormDelta(n, cdelta);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdelta);
    }
    (void)smacofUnweightedInitial(n, p, init, cdelta, cxini);
    if (DEBUG) {
        printf("xini\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cxini);
    }
    (void)smacofDistance(n, p, cxini, cdini);
    if (DEBUG) {
        printf("dini\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdini);
    }
    cxold = cxini;
    cdold = cdini;
    /*
    (void)memcpy(cxold, cxini, (sizeof *cxold) * n);
    (void)memcpy(cdold, cdini, (sizeof *cdold) * n);
    */
    if (DEBUG) {
        printf("xold\n\n");
        (void)smacofPrintAnyMatrix(n, p, width, precision, cxold);
    }
    if (DEBUG) {
        printf("dold\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cdold);
    }
    (void)smacofUnweightedMakeBMatrix(n, cdelta, cdold, cbold);
    if (DEBUG) {
        printf("bold\n\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, cbold);
    }
    (void)smacofUnweightedMakeStress(n, cdelta, cdold, &sold);
    if (DEBUG) {
        printf("sold %15.10f\n\n", sold);
    }
    while (true) {
        (void)smacofUnweightedGuttman(n, p, cbold, cxold, cxnew);
        if (DEBUG) {
            printf("xnew\n\n");
            (void)smacofPrintAnyMatrix(n, p, width, precision, cxnew);
        }
        (void)smacofRMSDifference(cxold, cxnew, n, p, &echange);
        /// this is confusing, skip when relax == FALSE
        (void)smacofRelax(n, p, cxold, cxnew, echange, pchange, &itel, relax);
        if (DEBUG) {
            printf("xrelaced\n\n");
            (void)smacofPrintAnyMatrix(n, p, width, precision, cxnew);
        }
        (void)smacofDistance(n, p, cxnew, cdnew);
        if (DEBUG) {
            printf("dnew\n\n");
            (void)smacofPrintAnyMatrix(n, n, width, precision, cdnew);
        }
        (void)smacofUnweightedMakeBMatrix(n, cdelta, cdnew, cbnew);
        if (DEBUG) {
            printf("bnew\n\n");
            (void)smacofPrintAnyMatrix(n, n, width, precision, cbnew);
        }
        (void)smacofUnweightedMakeStress(n, cdelta, cdnew, &snew);
        if (DEBUG) {
            printf("snew %15.10f\n\n", snew);
        }
        if (verbose) {
            printf(
                "itel %3d sold %12.10f snew %12.10f sdif %+12.10f rmsd "
                "%+12.10f rate "
                "%12.10f\n",
                itel, sold, snew, sold - snew, echange, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (echange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        pchange = echange;
        cxold = cxnew;
        cdold = cdnew;
        cbold = cbnew;
        /*
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        (void)memcpy(bold, bnew, (size_t)nn * sizeof(double));
        */
    }
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
