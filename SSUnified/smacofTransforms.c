#include "smacof.h"

void smacofInterval(const int n, double **delta, double **dmat, double **dhat) {
    double deltamin = INFINITY, alpha = 0.0, beta = 0.0, gamma = 0.0, s1 = 0.0,
           s2 = 0.0;
    double sed = 0.0, see = 0.0, se = 0.0, sd = 0.0, sdd = 0.0,
           dm = (double)(n * (n - 1) / 2);
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            deltamin = MIN(deltamin, delta[i][j]);
        }
    }
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            double eij = delta[i][j] - deltamin;
            double dij = dmat[i][j];
            se += eij;
            sd += dij;
            sed += eij * dij;
            see += SQUARE(eij);
            sdd += SQUARE(dij);
        }
    }
    se /= dm;
    sd /= dm;
    sed /= dm;
    see /= dm;
    sdd /= dm;
    alpha = (sed - se * sd) / (see - SQUARE(se));
    gamma = (sd - alpha * se);
    if ((alpha < 0.0) || (gamma < 0.0)) {
        s1 = sdd - SQUARE(sd);         // case alpha = 0
        s2 = sdd - SQUARE(sed) / see;  // case gamma = 0
        if (s1 < s2) {
            alpha = 0.0;
            gamma = sd;
        } else {
            alpha = sed / see;
            gamma = 0.0;
        }
    }
    beta = gamma - alpha * deltamin;
    for (int j = 0; j < (n - 1); j++) {
        for (int i = (j + 1); i < n; i++) {
            dhat[i][j] = alpha * delta[i][j] + beta;
            dhat[j][i] = dhat[i][j];
        }
    }
    return;
}
