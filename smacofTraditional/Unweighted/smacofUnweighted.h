#ifndef SMACOF_UNWEIGHTED_H
#define SMACOF_UNWEIGHTED_H

#include "../Common/smacofCommon.h"

void smacofSSMUEngine(const int n, const int p, double **delta, double **xini,
                      double **xnew, double **dnew, double **bnew,
                      const int init, const int itmax, const int ieps1,
                      const int ieps2, const bool verbose, const bool relax,
                      int *pitel, double *psnew);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const int n, const int p, double **delta,
                           double **dold, double **xold);
void smacofUnweightedNormDelta(const int n, double **delta);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const int n, double **delta, double **dmat,
                                 double **bmat);
void smacofUnweightedGuttman(const int n, const int p, double **bmat,
                             double **xold, double **xnew);
double smacofUnweightedMakeStress(const int n, double **delta, double **dist);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const int n, const int p, double **delta,
                                   double **xold);
void smacofUnweightedInitMaximumSum(const int n, const int p, double **delta,
                                    double **xold);
void smacofUnweightedInitial(const int n, const int p, const int init,
                             double **delta, double **xold);

#endif /* SMACOF_UNWEIGHTED_H */
