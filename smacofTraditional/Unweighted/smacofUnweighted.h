#ifndef SMACOF_UNWEIGHTED_H
#define SMACOF_UNWEIGHTED_H

#include "../Common/smacofCommon.h"

// smacofUnweightedNorm.c

void smacofUnweightedScale(const unsigned n, const unsigned p,
                           const double **delta, double **dold, double **xold);
void smacofUnweightedNormDelta(const unsigned n, double **delta);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const unsigned n, const double (*delta)[n][n],
                                 const double **smacofUnweightedMakeStress,
                                 double **bmat);
void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double **bmat, const double **xold, **xnew);
double smacofUnweightedMakeStress(const unsigned n, const double **delta,
                                  const double **dist);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const unsigned n, const unsigned p,
                                   const double (**delta, double **xold);
void smacofUnweightedInitMaximumSum(const unsigned n, const unsigned p,
                                    const double **delta, double **xold);
void smacofUnweightedInitial(const unsigned n, const unsigned p,
                             const unsigned init, const double **delta, double **xold);

#endif /* SMACOF_UNWEIGHTED_H */
