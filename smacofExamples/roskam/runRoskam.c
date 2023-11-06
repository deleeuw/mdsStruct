#include "../../smacofCommon/ccode/smacof.h"

int main() {
  double bnew[1176] = {0};
  double dnew[1128] = {0};
  double xini[96] = {0};
  double xnew[96] = {0};
  int irow[351] = {0};
  int icol[351] = {0};
  double delta[351] = {0};
  double weights[351] = {0};
  char fname[] = "roskam.txt";
  int n = 48, p = 2, itel = 1, itmax = 2, init = 2, eps1 = 15, eps2 = 10, m = 351;
  double snew = 0.0;
  bool verbose = true, relax = true;
  (void)smacofReadInputFile(fname, irow, icol, delta, weights);
  (void)smacofSSMWEngine(delta, weights, irow,
                  icol, xini, xnew,
                  dnew, bnew, &snew, &init,
                  &n, &p, &m, &itel,
                  &itmax, &eps1, &eps2,
                  &verbose, &relax);
  return(EXIT_SUCCESS);
}
