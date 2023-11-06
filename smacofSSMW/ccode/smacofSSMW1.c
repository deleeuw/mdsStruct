#include "../../smacofInclude/smacof.h"

int main() {
    double delta[15] = {1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0, 8.0,
                        9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
    double weights[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int irow[15] = {2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6};
    int icol[15] = {1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5};
    double xini[12] = {0.0};
    double xnew[12] = {0.0};
    double dnew[15] = {0.0};
    double bnew[21] = {0.0};
    double snew = 0.0;
    int init = 1;
    int n = 6;
    int p = 2;
    int m = 15;
    int itel = 1;
    int itmax = 1000;
    int eps1 = 15;
    int eps2 = 10;
    bool verbose = true;
    bool relax = true;
    int width = 15;
    int precision = 10;
    (void)smacofSSMWEngine(delta, weights, irow, icol, xini, xnew, dnew, bnew,
                           &snew, &init, &n, &p, &m, &itel, &itmax, &eps1,
                           &eps2, &verbose, &relax);
    printf("delta\n\n");
    (void)smacofPrintSHMatrix(delta, &n, &width, &precision);
    printf("weights\n\n");
    (void)smacofPrintSHMatrix(weights, &n, &width, &precision);
    printf("itel %4d snew %15.10f\n\n", itel, snew);
    printf("configuration\n\n");
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    printf("distances\n\n");
    (void)smacofPrintSHMatrix(dnew, &n, &width, &precision);
    return (EXIT_SUCCESS);
}