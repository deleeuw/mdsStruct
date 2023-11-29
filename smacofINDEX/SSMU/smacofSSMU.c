#include "../Unweighted/smacofUnweighted.h"

int main(void) {
    double delta[15] = {1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0, 8.0,
                        9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
    double xini[12] = {0.0};
    double xnew[12] = {0.0};
    double dnew[15] = {0.0};
    double bnew[21] = {0.0};
    double snew = 0.0;
    unsigned init = 1;
    unsigned n = 6;
    unsigned p = 2;
    unsigned itel = 1;
    unsigned itmax = 1000;
    unsigned eps1 = 15;
    unsigned eps2 = 10;
    bool verbose = true;
    bool relax = false;
    unsigned width = 15;
    unsigned precision = 10;
    (void)smacofSSMUEngine(n, p, delta, xini, xnew, dnew, bnew,
                           init, itmax, eps1,
                           eps2, verbose, relax, &itel, &snew);
    printf("delta\n\n");
    (void)smacofPrintSHMatrix(n, delta, width, precision);
    printf("itel %4d snew %15.10f\n\n", itel, snew);
    printf("configuration\n\n");
    (void)smacofPrintAnyMatrix(n, p, xnew, width, precision);
    printf("distances\n\n");
    (void)smacofPrintSHMatrix(n, dnew, width, precision);
    return (EXIT_SUCCESS);
}