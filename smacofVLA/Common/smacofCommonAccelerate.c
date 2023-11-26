#include "../Include/smacof.h"

void smacofRelax(const unsigned n, const unsigned p,
                 const double (*cxold)[n][p], double (*cxnew)[n][p],
                 const double chold, double *pchnew, double *prate,
                 const unsigned itel, const bool relax) {
    double chnew = *pchnew, rate = *prate, eopt = 0.0;
    chnew = smacofRMSDifference(n, p, cxold, cxnew);
    if (itel == 1) {
        rate = 1.0, eopt = 1.0;
    } else {
        rate = chnew / chold;
        eopt = MAX(0, MIN(1, rate / (2 - rate)));
    }
    if ((itel > 1) && relax) {
        for (unsigned i = 0; i < n; i++) {
            for (unsigned s = 0; s < p; s++) {
                (*cxnew)[i][s] =
                    (1 + eopt) * (*cxnew)[i][s] - eopt * (*cxold)[i][s];
            }
        }
    }
    *pchnew = chnew;
    *prate = rate;
    return;
}
