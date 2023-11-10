#include "../../smacofInclude/smacof.h"

// global variables

double delta[105];
double bnew[120];
double dnew[105];
double xini[30];
double xnew[30];
char fname[] = "eqdistU.txt";
int n = 15, p = 2, m = 105, itel = 1, itmax = 10000, init = 3, eps1 = 15,
    eps2 = 10, width = 6, precision = 4;
double snew = 0.0;
bool verbose = false, relax = true;
int main(void) {
    (void)smacofUnweightedReadInputFile(fname, delta);
    (void)smacofSSMUEngine(delta, xini, xnew, dnew, bnew, &snew, &init, &n, &p,
                           &m, &itel, &itmax, &eps1, &eps2, &verbose, &relax);
    return EXIT_SUCCESS;
}
