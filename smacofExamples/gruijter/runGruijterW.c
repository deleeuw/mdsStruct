#include "../../smacofInclude/smacof.h"

// global variables

double delta[36];
double bnew[45];
double dnew[36];
double xini[18];
double xnew[18];
double weights[36];
int irow[36];
int icol[36];
char fname[] = "gruijterW.txt";
int n = 9, p = 2, m = 36, itel = 1, itmax = 1000, init = 1, eps1 = 15,
    eps2 = 10, width = 6, precision = 4;
double snew = 0.0;
bool verbose = false, relax = true;

int main(void) {
    (void)smacofWeightedReadInputFile(fname, irow, icol, delta, weights);
    (void)smacofSSMWEngine(delta, weights, irow, icol, xini, xnew, dnew, bnew,
                           &snew, &init, &n, &p, &m, &itel, &itmax, &eps1,
                           &eps2, &verbose, &relax);
    return EXIT_SUCCESS;
}
