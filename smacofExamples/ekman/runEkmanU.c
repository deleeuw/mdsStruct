#include "../../smacofInclude/smacof.h"

// global variables

double delta[91];
double bnew[105];
double dnew[91];
double xini[28];
double xnew[28];
int irow[91];
int icol[91];
int ntimes = 10000;
char fname[] = "ekmanU.txt";
int n = 14, p = 2, m = 91, itel = 1, itmax = 100, init = 1, eps1 = 15,
    eps2 = 10, width = 6, precision = 4;
double snew;
bool verbose = false, relax = true, timer = true;

int main() {
    (void)smacofUnweightedReadInputFile(fname, irow, icol, delta);
    for (int i = 1; i <= ntimes; i++) {
        (void)smacofSSMUEngine(delta, irow, icol, xini, xnew, dnew, bnew, &snew,
                               &init, &n, &p, &m, &itel, &itmax, &eps1, &eps2,
                               &verbose, &relax);
    }
    return EXIT_SUCCESS;
}
