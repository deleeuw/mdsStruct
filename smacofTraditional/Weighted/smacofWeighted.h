#ifndef SMACOF_WEIGHTED_H
#define SMACOF_WEIGHTED_H

#include "../Common/smacofCommon.h"

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dold, double *bmat,
                               const unsigned *irow, const unsigned *icol,
                               const unsigned *pn, const unsigned *pm);
void smacofWeightedMakeVMatrix(const double *weights, double *vmat,
                               const unsigned *irow, const unsigned *icol,
                               const unsigned *pn, const unsigned *pm);
void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const unsigned *pn,
                           const unsigned *pp);
void smacofWeightedMakeStress(const double *delta, const double *weights,
                              const double *dist, const unsigned *irow,
                              const unsigned *icol, const unsigned *pn,
                              const unsigned *pm, double *stress);

// smacofWeightedNorm.c

void smacofWeightedScale(const int n, const int p, const double **delta,
                         const double **weights, double **dmat, double **xold);
void smacofWeightedNormDelta(const int n, double **delta, const double **weights);
void smacofWeightedNormWeights(const int n, double **weights);

// smacofWeightedInitial.c

void smacofWeightedInitTorgerson(const double *delta, const double *weights,
                                 const unsigned *irow, const unsigned *icol,
                                 double *xold, const unsigned *pn,
                                 const unsigned *pp, const unsigned *pm);
void smacofWeightedInitMaximumSum(const double *delta, const double *weights,
                                  const unsigned *irow, const unsigned *icol,
                                  double *xini, const unsigned *pn,
                                  const unsigned *pp, const unsigned *pm);
void smacofWeightedInitial(const double *delta, const double *weights,
                           const unsigned *irow, const unsigned *icol,
                           double *xini, const unsigned *pinit,
                           const unsigned *pn, const unsigned *pp,
                           const unsigned *pm);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const unsigned *pn, const unsigned *pp);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v,
                           const unsigned *pn, const unsigned *pp);

#endif /* SMACOF_WEIGHTED_H */
