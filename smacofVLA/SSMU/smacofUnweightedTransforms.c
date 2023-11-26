#include "../Include/smacof.h"

void smacofUnweightedInterval(const unsigned n, const double (*delta)[n][n],
                              const double (*dmat)[n][n],
                              double (*dhat)[n][n]) {
    double deltamin = INFINITY;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            deltamin = MIN(deltamin, (*delta)[i][j]);
        }
    }
    double sed = 0.0, see = 0.0, se = 0.0, sd = 0.0, sdd = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            double eij = (*delta)[i][j] - deltamin;
            double dij = (*dmat)[i][j];
            se += eij;
            sd += dij;
            sed += eij * dij;
            see += SQUARE(eij);
            sdd += SQUARE(dij);
        }
    }
    double alpha = (sed - se * sd) / (see - SQUARE(se));
    double gamma = sd - alpha * se;
    if ((alpha < 0.0) || (gamma < 0.0)) {
        double s1 = sdd - (SQUARE(sed) / see);
        double s2 = sdd - SQUARE(sd);
        if (s1 <= s2) {
            alpha = 0.0;
            gamma = sd;
        } else {
            alpha = sed / see;
            gamma = 0.0;
        }
    }
    double beta = gamma - alpha * deltamin;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            (*dhat)[i][j] = alpha * (*delta)[i][j] + beta;
        }
    }
    return;
}
