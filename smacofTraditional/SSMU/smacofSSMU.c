#include "../Unweighted/smacofUnweighted.h"

int main(void) {
    (void)smacofSSMUEngine(n, p, delta, xini, xnew, dnew, bnew, init, itmax,
                           eps1, eps2, verbose, relax, &itel, &snew);
    return EXIT_SUCCESS;
}
