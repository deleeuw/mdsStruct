#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG false
#define NaN (nan(""))
#define PI (2.0 * asin(1.0))
#define EPS(i) ((double)pow(10.0, -(i)))

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofEngine.c

void smacofSSMWEngine(double *delta, double *weights, const int *irow,
                      const int *icol, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax);

void smacofSSMUEngine(double *delta, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax);

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dold, double *bmat,
                               const int *irow, const int *icol, const int *pn,
                               const int *pm);
void smacofWeightedMakeVMatrix(const double *weights, double *vmat,
                               const int *irow, const int *icol, const int *pn,
                               const int *pm);
void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const int *pn,
                           const int *pp);
void smacofWeightedMakeStress(const double *delta, const double *weights,
                              const double *dist, const int *irow,
                              const int *icol, const int *pn, const int *pm,
                              double *stress);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const double *delta, const double *dmat,
                                 double *bmat, const int *pn);
void smacofUnweightedGuttman(const double *bmat, const double *xold,
                             double *xnew, const int *pn, const int *pp);
void smacofUnweightedMakeStress(const double *delta, const double *dist,
                                double *stress, const int *pn);

// smacofWeightedNorm.c

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const int *pn,
                         const int *pp, const int *pm);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const int *pm);
void smacofWeightedNormWeights(double *weights, const int *pm);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const double *delta, double *dold, double *xold,
                           const int *pn, const int *pp);
void smacofUnweightedNormDelta(double *delta, const int *pn);

// smacofCommonUtils.c

void smacofInitRandom(double *xini, const int *pn, const int *pp);
void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      const int *pn, const int *pp,
                                      double *maxdiff);
void smacofMaxDistanceDifference(const double *x, const double *y,
                                 const int *pm, double *maxdiff);
void smacofRMSDifference(const double *x, double *y, const int *pn,
                         const int *pp, double *diff);
void smacofMakeIIandJJ(const int *pn, int *ii, int *jj);
void smacofDoubleCenter(const double *delta, double *cross, const int *pn);

// smacofSort.c

int smacofComparison(const void *px, const void *py);
void smacofWeightedSort(double *delta, double *weights, int *knum, int *irow,
                        int *icol, const int *pn);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

// smacofWeightedInitial.c

void smacofWeightedInitTorgerson(const double *delta, const double *weights,
                                 const int *irow, const int *icol, double *xold,
                                 const int *pn, const int *pp, const int *pm);
void smacofWeightedInitMaximumSum(const double *delta, const double *weights,
                                  const int *irow, const int *icol,
                                  double *xini, const int *pn, const int *pp,
                                  const int *pm);
void smacofWeightedInitial(const double *delta, const double *weights,
                           const int *irow, const int *icol, double *xini,
                           const int *pinit, const int *pn, const int *pp,
                           const int *pm);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const double *delta, double *xold,
                                   const int *pn, const int *pp);
void smacofUnweightedInitMaximumSum(const double *delta, double *xini,
                                    const int *pn, const int *pp);
void smacofUnweightedInitial(const double *delta, double *xini,
                             const int *pinit, const int *pn, const int *pp);

// smacofLinearAlgebra.c

void smacofJacobi(double *a, double *evec, double *eval, const int *pn,
                  const int *pm, const int *itmax, const int *eps,
                  const bool *verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const int *pn);
void smacofGramSchmidt(double *x, double *r, int *pn, int *pp);
void smacofMultiplySymmetricMatrix(const double *a, const double *x, double *y,
                                   const int *pn, const int *pp);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv, const int *ndim);
void smacofDistance(const double *x, double *d, const int *pn, const int *pp);
void smacofCenter(double *x, const int *np, const int *pp);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const int *pn, const int *pp);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const int *pn,
                           const int *pp);

// smacofCommonPrintRead.c

void smacofPrintSDCLMatrix(const double *v, const int *ndim, const int *pw,
                           const int *pr);
void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp,
                          const int *pw, const int *pr);
void smacofPrintSHMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr);
void smacofPrintSHMatrixIJ(const double *d, const int *pn, int *pm,
                           const int *irow, const int *icol, const int *pw,
                           const int *pr);
void smacofPrintLTMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr);
void smacofPrintSLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr);
void smacofPrintSymmetricMatrix(const double *x, const int *pn, const int *pw,
                                const int *pr);
void smacofWeightedReadInputFile(char *fname, int *irow, int *icol,
                                 double *delta, double *weights);
void smacofUnweightedReadInputFile(char *fname, double *delta);

// smacofCommonAccelerate.c

void smacofRelax(const double *xold, double *xnew, const double *pechange,
                 const double *ppchange, const int *pnp, const int *pitel,
                 const bool *prelax, double *prate);

// smacofCommonIndices.c

int SINDEX(const int i, const int j, const int n);
int TINDEX(const int i, const int j, const int n);
int PINDEX(const int i, const int j, const int n);
int UINDEX(const int i, const int j, const int n);

// structures

struct data {
    int order;
    int irow;
    int icol;
    double delta;
    double weight;
};

#endif /* SMACOF_H */
