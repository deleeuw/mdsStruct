#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG true
#define NaN (nan(""))
#define PI (2.0 * asin(1.0))
#define EPS(i) ((double)pow(10.0, -(i)))

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofRCTranslation.c

void smacofFromAnyRtoC(const double *rmatrix, const size_t nrow,
                       const size_t ncol, double cmatrix[nrow][ncol]);
void smacofFromSymmetricHollowRtoC(const double *rmatrix, const size_t order,
                                   double cmatrix[order][order]);
void smacofFromSymmetricRtoC(const double *rmatrix, const size_t order,
                             double cmatrix[order][order]);

// smacofEngine.c

void smacofSSMWEngine(double *delta, double *weights, const size_t *irow,
                      const size_t *icol, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const size_t *pinit,
                      const size_t *pn, const size_t *pp, const size_t *pm, size_t *pitel,
                      const size_t *pitmax, const size_t *peps1, const size_t *peps2,
                      const bool *pverbose, const bool *prelax);

void smacofSSMUEngine(const size_t n, const size_t p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const size_t init, const size_t itmax,
                      const size_t eps1, const size_t eps2,
                      const bool verbose, const bool relax, size_t itel,
                      double *psnew);

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dold, double *bmat,
                               const size_t *irow, const size_t *icol, const size_t *pn,
                               const size_t *pm);
void smacofWeightedMakeVMatrix(const double *weights, double *vmat,
                               const size_t *irow, const size_t *icol, const size_t *pn,
                               const size_t *pm);
void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const size_t *pn,
                           const size_t *pp);
void smacofWeightedMakeStress(const double *delta, const double *weights,
                              const double *dist, const size_t *irow,
                              const size_t *icol, const size_t *pn, const size_t *pm,
                              double *stress);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const size_t n, const double delta[n][n],
                                 const double dmat[n][n], double bmat[n][n]);
void smacofUnweightedGuttman(const size_t n, const size_t p,
                             const double bmat[n][n], const double xold[n][p],
                             double xnew[n][p]);
void smacofUnweightedMakeStress(const size_t n, const double delta[n][n],
                                const double dist[n][n], double *stress);

// smacofWeightedNorm.c

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const size_t *pn,
                         const size_t *pp, const size_t *pm);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const size_t *pm);
void smacofWeightedNormWeights(double *weights, const size_t *pm);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const size_t n, const size_t p,
                           const double delta[n][n], double dold[n][n],
                           double xold[n][p]);
void smacofUnweightedNormDelta(const size_t n, double delta[n][n]);

// smacofCommonUtils.c

void smacofInitRandom(const size_t n, const size_t p, double xini[n][p]);
void smacofMaxConfigurationDifference(const size_t n, const size_t p,
                                      const double x[n][p],
                                      const double y[n][p], double *maxdiff);
void smacofMaxDistanceDifference(const size_t n, const double dold[n][n],
                                 const double dnew[n][n], double *pdchange);
void smacofRMSDifference(const size_t n, const size_t p,
                         const double x[n][p], double y[n][p], double *diff);
void smacofDoubleCenter(const size_t n, const double delta[n][n],
                        double cross[n][n]);

// smacofSort.c

int smacofComparison(const void *px, const void *py);
void smacofWeightedSort(double *delta, double *weights, size_t *knum, size_t *irow,
                        size_t *icol, const size_t *pn);
void smacofTieBlocks(const double *x, size_t *it, double *eps, const size_t *n);

// smacofWeightedInitial.c

void smacofWeightedInitTorgerson(const double *delta, const double *weights,
                                 const size_t *irow, const size_t *icol, double *xold,
                                 const size_t *pn, const size_t *pp, const size_t *pm);
void smacofWeightedInitMaximumSum(const double *delta, const double *weights,
                                  const size_t *irow, const size_t *icol,
                                  double *xini, const size_t *pn, const size_t *pp,
                                  const size_t *pm);
void smacofWeightedInitial(const double *delta, const double *weights,
                           const size_t *irow, const size_t *icol, double *xini,
                           const size_t *pinit, const size_t *pn, const size_t *pp,
                           const size_t *pm);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const size_t n, const size_t p,
                                   const double delta[n][n], double xold[n][p]);
void smacofUnweightedInitMaximumSum(const size_t n, const size_t p,
                                    const double delta[n][n],
                                    double xold[n][p]);
void smacofUnweightedInitial(const size_t n, const size_t p,
                             const size_t init, const double delta[n][n],
                             double xold[n][p]);

// smacofCommonLinearAlgebra.c

void smacofJacobi(const size_t n, const size_t p, double a[n][n],
                  double evec[n][n], double *eval, const size_t itmax,
                  const size_t eps, const bool verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const size_t *pn);
void smacofGramSchmidt(const size_t n, const size_t p, double x[n][p],
                       double r[n]);
void smacofMultiplySymmetricMatrix(const size_t n, const size_t p,
                                   const double a[n][n], const double x[n][p],
                                   double y[n][p]);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv, const size_t *ndim);
void smacofDistance(const size_t n, const size_t p, const double x[n][p],
                    double d[n][n]);
void smacofCenter(const size_t n, const size_t p, double x[n][p]);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const size_t *pn, const size_t *pp);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const size_t *pn,
                           const size_t *pp);

// smacofCommonPrintRead.c

void smacofPrintSDCLMatrix(const double *v, const size_t *ndim, const size_t *pw,
                           const size_t *pr);
void smacofPrintAnyMatrix(const size_t n, const size_t p,
                          const int width, const int precision,
                          double x[n][p]);
void smacofPrintSHMatrix(const double *d, const size_t *pn, const size_t *pw,
                         const size_t *pr);
void smacofPrintSHMatrixIJ(const double *d, const size_t *pn, size_t *pm,
                           const size_t *irow, const size_t *icol, const size_t *pw,
                           const size_t *pr);
void smacofPrintLTMatrix(const double *d, const size_t *pn, const size_t *pw,
                         const size_t *pr);
void smacofPrintSLTMatrix(const double *d, const size_t *pn, const size_t *pw,
                          const size_t *pr);
void smacofPrintSymmetricMatrix(const double *x, const size_t *pn, const size_t *pw,
                                const size_t *pr);
void smacofWeightedReadInputFile(char *fname, size_t *irow, size_t *icol,
                                 double *delta, double *weights);
void smacofUnweightedReadInputFile(char *fname, double *delta);

// smacofCommonAccelerate.c

void smacofRelax(const size_t n, const size_t p, const double xold[n][p],
                 double xnew[n][p], const double echange, const double pchange,
                 const size_t itel, const bool relax, double *prate);

// smacofCommonIndices.c

size_t SINDEX(const size_t i, const size_t j, const size_t n);
size_t TINDEX(const size_t i, const size_t j, const size_t n);
size_t PINDEX(const size_t i, const size_t j, const size_t n);
size_t UINDEX(const size_t i, const size_t j, const size_t n);

// structures

struct data {
    size_t order;
    size_t irow;
    size_t icol;
    double delta;
    double weight;
};

#endif /* SMACOF_H */
