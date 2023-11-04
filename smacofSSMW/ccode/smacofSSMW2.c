#include <smacof.h>

int main() {
    double delta[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 0.0};
    double weights[6] = {0.0, 1.0, 1.0, 1.0, 1.0, 0.0};
    int irow[6] = {2, 3, 4, 3, 4, 4};
    int icol[6] = {1, 1, 1, 2, 2, 3};
    double xini[8] = {0.0};
    double xnew[8] = {0.0};
    double dnew[6] = {0.0};
    double bnew[10] = {0.0};
    double snew = 0.0;
    int init = 1;
    int n = 4;
    int p = 2;
    int m = 6;
    int itel = 1;
    int itmax = 10000;
    int eps1 = 15;
    int eps2 = 10;
    bool verbose = false;
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
    (void)printf("itel %4d snew %15.10f\n\n", itel, snew);
    printf("configuration\n\n");
    (void)smacofPrintAnyMatrix(xnew, &n, &p, &width, &precision);
    printf("distances\n\n");
    (void)smacofPrintSHMatrix(dnew, &n, &width, &precision);
    return (EXIT_SUCCESS);
}