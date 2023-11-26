#ifndef SMACOF_H
#define SMACOF_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG false
#define PI (2.0 * asin(1.0))

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofCommonRCTranslation.c  // done

void smacofFromAnyRtoC(const unsigned nrow, const unsigned ncol,
                       const double *rmatrix, double (*cmatrix)[nrow][ncol]);
void smacofFromSymmetricHollowRtoC(const unsigned n, const double *rmatrix,
                                   double (*cmatrix)[n][n]);
void smacofFromSymmetricRtoC(const unsigned n, const double *rmatrix,
                             double (*cmatrix)[n][n]);
void smacofFromCtoR(const unsigned nrow, const unsigned ncol,
                    const double (*cmatrix)[nrow][ncol], double *rmatrix);

// smacofEngine.c

void smacofSSMWEngine(double *delta, double *weights, const unsigned *irow,
                      const unsigned *icol, double *xini, double *xnew,
                      double *dnew, double *bnew, double *psnew,
                      const unsigned *pinit, const unsigned *pn,
                      const unsigned *pp, const unsigned *pm, unsigned *pitel,
                      const unsigned *pitmax, const unsigned *peps1,
                      const unsigned *peps2, const bool *pverbose,
                      const bool *prelax);
void smacofSSMUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      const unsigned init, const unsigned itmax, unsigned ieps1,
                      unsigned ieps2, const bool verbose, const bool relax,
                      unsigned *pitel, double *psnew);
void smacofSSIUEngine(const unsigned n, const unsigned p, double *delta,
                      double *xini, double *xnew, double *dnew, double *dhat,
                      double *bnew, const unsigned init, const unsigned itmax,
                      const unsigned ieps1, const unsigned ieps2,
                      const bool verbose, const bool relax, unsigned *pitel,
                      double *psnew, const bool interval);

// smacofCommonUtils.c

void smacofInitRandom(const unsigned n, const unsigned p, double (*xini)[n][p]);
void smacofDoubleCenter(const unsigned n, const double (*delta)[n][n],
                        double (*cross)[n][n]);
double smacofMaxConfigurationDifference(const unsigned n, const unsigned p,
                                        const double (*x)[n][p],
                                        const double (*y)[n][p]);
double smacofMaxDistanceDifference(const unsigned n, const double (*dold)[n][n],
                                   const double (*dnew)[n][n]);
double smacofRMSDifference(const unsigned n, const unsigned p,
                           const double (*x)[n][p], const double (*y)[n][p]);
void smacofZeroAnyMatrix(const unsigned n, const unsigned p, double (*x)[n][p]);

// smacofCommonLinearAlgebra.c

void smacofJacobi(const unsigned n, const unsigned p, double (*a)[n][n],
                  double (*evec)[n][p], double (*eval)[n], const unsigned itmax,
                  const unsigned eps, const bool verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const unsigned *pn);
void smacofGramSchmidt(const unsigned n, const unsigned p, double (*x)[n][p],
                       double (*r)[n]);
void smacofMultiplySymmetricMatrix(const unsigned n, const unsigned p,
                                   const double (*a)[n][n],
                                   const double (*x)[n][p], double (*y)[n][p]);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv,
                               const unsigned *ndim);
void smacofDistance(const unsigned n, const unsigned p, const double (*x)[n][p],
                    double (*d)[n][n]);
void smacofCenter(const unsigned n, const unsigned p, double (*x)[n][p]);

// smacofCommonPrintRead.c

void smacofPrintAnyMatrix(const unsigned n, const unsigned p,
                          const unsigned width, const unsigned precision,
                          const double (*x)[n][p]);
void smacofUnweightedReadInputFile(const char *fname, double *delta);

void smacofUnweightedInterval(const unsigned n, const double (*delta)[n][n],
                              const double (*dmat)[n][n], double (*dhat)[n][n]);

/*
void smacofPrintSDCLMatrix(const double *v, const unsigned *ndim, const unsigned
*pw, const unsigned *pr); void smacofPrintSHMatrix(const double *d, const
unsigned *pn, const unsigned *pw, const unsigned *pr); void
smacofPrintSHMatrixIJ(const double *d, const unsigned *pn, unsigned *pm, const
unsigned *irow, const unsigned *icol, const unsigned *pw, const unsigned *pr);
void smacofPrintLTMatrix(const double *d, const unsigned *pn, const unsigned
*pw, const unsigned *pr); void smacofPrintSLTMatrix(const double *d, const
unsigned *pn, const unsigned *pw, const unsigned *pr); void
smacofPrintSymmetricMatrix(const double *x, const unsigned *pn, const unsigned
*pw, const unsigned *pr); void smacofWeightedReadInputFile(char *fname, unsigned
*irow, unsigned *icol, double *delta, double *weights);
*/

// smacofCommonAccelerate.c

void smacofRelax(const unsigned n, const unsigned p,
                 const double (*cxold)[n][p], double (*cxnew)[n][p],
                 const double chold, double *pchnew, double *prate,
                 const unsigned itel, const bool relax);

// smacofCommonIndices.c

/*
 * int SINDEX(const unsigned i, const unsigned j, const unsigned n);
int TINDEX(const unsigned i, const unsigned j, const unsigned n);
int PINDEX(const unsigned i, const unsigned j, const unsigned n);
int UINDEX(const unsigned i, const unsigned j, const unsigned n);
 */

// smacofUnweightedNorm.c

void smacofUnweightedScale(const unsigned n, const unsigned p,
                           const double (*delta)[n][p], double (*dold)[n][n],
                           double (*xold)[n][p]);
void smacofUnweightedNormDelta(const unsigned n, double (*delta)[n][n]);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const unsigned n, const double (*delta)[n][n],
                                 const double (*dmat)[n][n],
                                 double (*bmat)[n][n]);
void smacofUnweightedGuttman(const unsigned n, const unsigned p,
                             const double (*bmat)[n][n],
                             const double (*xold)[n][p], double (*xnew)[n][p]);
double smacofUnweightedMakeStress(const unsigned n, const double (*delta)[n][n],
                                  const double (*dist)[n][n]);

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const unsigned n, const unsigned p,
                                   const double (*delta)[n][n],
                                   double (*xold)[n][p]);
void smacofUnweightedInitMaximumSum(const unsigned n, const unsigned p,
                                    const double (*delta)[n][n],
                                    double (*xold)[n][p]);
void smacofUnweightedInitial(const unsigned n, const unsigned p,
                             const unsigned init, const double (*delta)[n][n],
                             double (*xold)[n][p]);

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

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const unsigned *pn,
                         const unsigned *pp, const unsigned *pm);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const unsigned *pm);
void smacofWeightedNormWeights(double *weights, const unsigned *pm);

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

// structures

struct data {
    unsigned order;
    unsigned irow;
    unsigned icol;
    double delta;
    double weight;
};

#endif /* SMACOF_H */
