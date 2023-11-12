#include "smacof.h"

// to be called from R

double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
double xini[8];
double xnew[8];
double dnew[6];
double bnew[10];
double snew = 0.0;
rsize_t init = 1, n = 4, p = 2, itel = 1, itmax = 5
int eps1 = 15, eps2 = 10;
bool verbose = true, relax = true;

int main(void) {
    (void)smacofSSMUEngine(
        const size_t n, const size_t p, double *delta, double *xini,
        double *xnew, double *dnew, double *bnew, const size_t init,
        const size_t itmax, const int eps1, const int eps2,
        const bool verbose, const bool relax, size_t itel, double *snew);
    return EXIT_SUCCESS;
}

void smacofSSMUEngine(const size_t n, const size_t p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const size_t init, const size_t itmax,
                      const size_t eps1, const size_t eps2,
                      const bool verbose, const bool relax, size_t itel,
                      double *psnew) {
    size_t np = p * n, nn = SQUARE(n);
    size_t width = 15, precision = 10;
    double sold = 0.0, snew = *psnew, pchange = 1.0, echange = 1.0, rate = 1.0;
    double eps1 = EPS(eps1), eps2 = EPS(eps2);
    /*
     * Convert delta from an R vector to a C matrix
     */
    double(*cdelta)[n][n] = malloc((sizeof *cdelta) * nn);
    (void)smacofFromSymmetricHollowRtoC(delta, n, *cdelta);
    /*
     * Make some more C matrices

    double(*cdnew)[n][n] = malloc((sizeof *cdnew) * nn);
    double(*cbnew)[n][n] = malloc((sizeof *cbnew) * nn);
    double(*cxini)[n][p] = malloc((sizeof *cxini) * np);
    double(*cxnew)[n][p] = malloc((sizeof *cxnew) * np);
    double(*cxold)[n][p] = malloc((sizeof *cxold) * np);
    double(*cdold)[n][n] = malloc((sizeof *cdold) * nn);
    double(*cdini)[n][n] = malloc((sizeof *cdini) * nn);
    double(*cbold)[n][n] = malloc((sizeof *cbold) * nn);
    (void)smacofUnweightedNormDelta(n, *cdelta);
    if (DEBUG) {
        printf("delta\n\n");
        (void)smacofPrsize_tAnyMatrix(n, n, width, precision, *cdelta);
    }
    (void)smacofUnweightedInitial(n, p, init, *delta, *xini);
    if (DEBUG) {
        printf("xini\n\n");
        (void)smacofPrsize_tAnyMatrix(n, p, width, precision, *xini);
    }
    (void)smacofDistance(xini, dini, pn, pp);
    if (DEBUG) {
        printf("dini\n\n");
        (void)smacofPrsize_tSHMatrix(dini, &n, &width, &precision);
    }
    (void)memcpy(xold, xini, (size_t)np * sizeof(double));
    (void)memcpy(dold, dini, (size_t)m * sizeof(double));
    (void)smacofUnweightedMakeBMatrix(delta, dold, bold, pn);
    if (DEBUG) {
        printf("bold\n\n");
        (void)smacofPrsize_tSymmetricMatrix(bold, &n, &width, &precision);
    }
    (void)smacofUnweightedMakeStress(delta, dold, &sold, pn);
    if (DEBUG) {
        printf("sold %15.10f\n\n", sold);
    }
    while (true) {
        (void)smacofUnweightedGuttman(bold, xold, xnew, pn, pp);
        (void)smacofRMSDifference(xold, xnew, pn, pp, &echange);
        /// this is confusing, skip when relax == FALSE
        (void)smacofRelax(xold, xnew, &echange, &pchange, &np, &itel, prelax,
                          &rate);
        (void)smacofDistance(xnew, dnew, pn, pp);
        (void)smacofUnweightedMakeBMatrix(delta, dnew, bnew, pn);
        (void)smacofUnweightedMakeStress(delta, dnew, &snew, pn);
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
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        (void)memcpy(bold, bnew, (size_t)nn * sizeof(double));
    }
    *psnew = snew;
    *pitel = itel;
    free(dold);
    free(dini);
    free(xold);
    free(bold);
    return;
    }
    */
}
