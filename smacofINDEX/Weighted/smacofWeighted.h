#ifndef SMACOF_WEIGHTED_H
#define SMACOF_WEIGHTED_H

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dold, double *bmat,
                               const unsigned *irow, const unsigned *icol, const unsigned n,
                               const unsigned m);
void smacofWeightedMakeVMatrix(const double *weights, double *vmat,
                               const unsigned *irow, const unsigned *icol, const unsigned n,
                               const unsigned m);
void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const unsigned n,
                           const unsigned p);
void smacofWeightedMakeStress(const double *delta, const double *weights,
                              const double *dist, const unsigned *irow,
                              const unsigned *icol, const unsigned n, const unsigned m,
                              double *stress);



// smacofWeightedNorm.c

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const unsigned n,
                         const unsigned p, const unsigned m);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const unsigned m);
void smacofWeightedNormWeights(double *weights, const unsigned m);

// smacofWeightedInitial.c

void smacofWeightedInitTorgerson(const double *delta, const double *weights,
                                 const unsigned *irow, const unsigned *icol, double *xold,
                                 const unsigned n, const unsigned p, const unsigned m);
void smacofWeightedInitMaximumSum(const double *delta, const double *weights,
                                  const unsigned *irow, const unsigned *icol,
                                  double *xini, const unsigned n, const unsigned p,
                                  const unsigned m);
void smacofWeightedInitial(const double *delta, const double *weights,
                           const unsigned *irow, const unsigned *icol, double *xini,
                           const unsigned *pinit, const unsigned n, const unsigned p,
                           const unsigned m);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const unsigned n, const unsigned p);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const unsigned n,
                           const unsigned p);
//

void smacofWeightedReadInputFile(char *fname, unsigned *irow, unsigned *icol,
                                 double *delta, double *weights);





#endif /* SMACOF_WEIGHTED_H */