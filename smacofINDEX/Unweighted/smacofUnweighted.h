#ifndef SMACOF_UNWEIGHTED_H
#define SMACOF_UNWEIGHTED_H

#include "../Common/smacofCommon.h"

void smacofSSMUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const unsigned init, const unsigned itmax,
                      const unsigned ieps1, const unsigned ieps2,
                      const bool verbose, const bool relax, unsigned *pitel,
                      double *psnew);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const unsigned n, const double *delta,
                                 const double *dmat, double *bmat);
void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double *bmat, const double *xold,
                             double *xnew);
double smacofUnweightedMakeStress(const unsigned n, const double *delta,
                                  const double *dist);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const unsigned n, const unsigned p,
                           const double *delta, double *dold, double *xold);
void smacofUnweightedNormDelta(const unsigned n, double *delta);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const unsigned n, const unsigned p,
                                   const double *delta, double *xold);
void smacofUnweightedInitMaximumSum(const unsigned n, const unsigned p,
                                    const double *delta, double *xini);
void smacofUnweightedInitial(const unsigned n, const unsigned p,
                             const double *delta, double *xini,
                             const unsigned init);

// smacofUnweightedRead.c

void smacofUnweightedReadInputFile(char *fname, double *delta);

#endif /* SMACOF_UNWEIGHTED_H */
