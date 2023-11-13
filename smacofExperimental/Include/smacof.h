#ifndef SMACOF_H
#define SMACOF_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG true
#define NaN (nan(""))
#define PI (2.0 * asin(1.0))

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofCommonRCTranslation.c

void smacofFromAnyRtoC(const int nrow, const int ncol, const double *rmatrix,
                       double (*cmatrix)[nrow]);
void smacofFromSymmetricHollowRtoC(const int n, const double *rmatrix,
                                   double (*cmatrix)[n]);
void smacofFromSymmetricRtoC(const int n, const double *rmatrix,
                             double (*cmatrix)[n]);
void smacofFromCtoR(const int nrow, const int ncol,
                    const double (*cmatrix)[nrow], double *rmatrix);

// smacofEngine.c

void smacofSSMWEngine(double *delta, double *weights, const int *irow,
                      const int *icol, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax);

void smacofSSMUEngine(const int n, const int p, double *delta, double *xini,
                      double *xnew, double *dnew, double *bnew, const int init,
                      const int itmax, const int ieps1, const int ieps2,
                      const bool verbose, const bool relax, int *pitel,
                      double *psnew);

// smacofCommonUtils.c

void smacofInitRandom(const int n, const int p, double (*xini)[n]);
void smacofDoubleCenter(const int n, const double (*delta)[n],
                        double (*cross)[n]);
double smacofMaxConfigurationDifference(const int n, const int p,
                                        const double (*x)[n],
                                        const double (*y)[n]);
double smacofMaxDistanceDifference(const int n, const double (*dold)[n],
                                   const double (*dnew)[n]);
double smacofRMSDifference(const int n, const int p, const double x[n][p],
                           double (*y)[n]);

// smacofCommonSort.c

int smacofComparison(const void *px, const void *py);
void smacofWeightedSort(double *delta, double *weights, int *knum, int *irow,
                        int *icol, const int *pn);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

// smacofCommonLinearAlgebra.c

void smacofJacobi(const int n, const int p, double (*a)[n], double (*evec)[n],
                  double *eval, const int itmax, const int eps,
                  const bool verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const int *pn);
void smacofGramSchmidt(const int n, const int p, double (*x)[n], double r[n]);
void smacofMultiplySymmetricMatrix(const int n, const int p,
                                   const double (*a)[n], const double (*x)[n],
                                   double (*y)[n]);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv, const int *ndim);
void smacofDistance(const int n, const int p, const double (*x)[n],
                    double (*d)[n]);
void smacofCenter(const int n, const int p, double (*x)[n]);

// smacofCommonPrintRead.c

void smacofPrintSDCLMatrix(const double *v, const int *ndim, const int *pw,
                           const int *pr);
void smacofPrintAnyMatrix(const int n, const int p, const int width,
                          const int precision, double (*x)[n]);
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

void smacofRelax(const int n, const int p, const double (*xold)[n],
                 double (*xnew)[n], const double echange, const double pchange,
                 const int itel, const bool relax);
double smacofRate(const double echange, const double pchange, const int itel);

// smacofCommonIndices.c

int SINDEX(const int i, const int j, const int n);
int TINDEX(const int i, const int j, const int n);
int PINDEX(const int i, const int j, const int n);
int UINDEX(const int i, const int j, const int n);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const int n, const int p, const double (*delta)[n],
                           double (*dold)[n], double (*xold)[n]);
void smacofUnweightedNormDelta(const int n, double (*delta)[n]);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const int n, const double (*delta)[n],
                                 const double (*dmat)[n], double (*bmat)[n]);
void smacofUnweightedGuttman(const int n, const int p, const double (*bmat)[n],
                             const double (*xold)[n], double (*xnew)[n]);
void smacofUnweightedMakeStress(const int n, const double (*delta)[n],
                                const double (*dist)[n], double *stress);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const int n, const int p,
                                   const double delta[n][n], double xold[n][p]);
void smacofUnweightedInitMaximumSum(const int n, const int p,
                                    const double delta[n][n],
                                    double xold[n][p]);
void smacofUnweightedInitial(const int n, const int p, const int init,
                             const double delta[n][n], double xold[n][p]);

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

// smacofWeightedNorm.c

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const int *pn,
                         const int *pp, const int *pm);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const int *pm);
void smacofWeightedNormWeights(double *weights, const int *pm);

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

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const int *pn, const int *pp);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const int *pn,
                           const int *pp);

// structures

struct data {
  int order;
  int irow;
  int icol;
  double delta;
  double weight;
};

#endif /* SMACOF_H */
