#include "../../smacofInclude/smacof.h"

// global variables

double delta[630];
double weights[630];
int irow[630];
int icol[630];
double bnew[666];
double dnew[630];
double xini[72];
double xnew[72];
char fname[] = "morseW.txt";
int n = 36, p = 2, m = 630, itel = 1, itmax = 1000, init = 1, eps1 = 15,
  eps2 = 10, width = 2, precision = 2;
double snew = 0.0;
bool verbose = false, relax = true;

int main(void) {
    (void)smacofWeightedReadInputFile(fname, irow, icol, delta, weights);
    (void)smacofSSMWEngine(delta, weights, irow, icol, xini, xnew, dnew, bnew,
                           &snew, &init, &n, &p, &m, &itel, &itmax, &eps1,
                           &eps2, &verbose, &relax);
    return EXIT_SUCCESS;
}
