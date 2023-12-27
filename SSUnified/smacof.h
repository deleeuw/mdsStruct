#ifndef SMACOF_H
#define SMACOF_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <plot.h>

#define PI (2.0 * asin(1.0))
#define SSIZE 80
#define OSIZE 256
#define HAVE_INIT_CONFIGURATION 4
#define HAVE_INNER_KNOTS 1

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

void smacofSSEngine(const int n, const int p, double **delta, double **w,
                    double **vmat, double **vinv, double **xold, double **xnew,
                    double **dmat, double **dhat, double **basis,
                    const int haveinit, const int typeinit, const int itmax,
                    const int ieps1, const int ieps2, const bool verbose,
                    const bool relax, const int ditmax, const int ieps3,
                    const bool dverbose, const bool weights, const int degree,
                    const int ordinal, const int transform, char *iterstring);

// smacofIndices.c

int VINDEX(const int i);
int MINDEX(const int i, const int j, const int n);
int SINDEX(const int i, const int j, const int n);
int PINDEX(const int i, const int j, const int n);
int TINDEX(const int i, const int j, const int n);
int UINDEX(const int i, const int j, const int n);

// smacofCCD.c

void smacofCCD(const int n, const int degree, double *dmatvec, double *wvec,
               double *b, double *dhatvec, double **z, const int itmax,
               const int eps, const bool verbose, const bool weights,
               const bool ordinal);

// smacofDerivatives.c

void smacofGradient(const int n, const int p, const double *delta,
                    const double *dold, const double *xold, double *xnew,
                    double *gradient);
void smacofHessian(const int n, const int p, const double *delta,
                   const double *xconf, const double *dmat, const double *bmat);

// smacofRCTranslation.c

void smacofAnyRtoC(const int nrow, const int ncol, const double *rmatrix,
                   double **cmatrix);
void smacofFromSymmetricHollowRtoC(const int n, const double *rmatrix,
                                   double **cmatrix);
void smacofSymmetricRtoC(const int n, const double *rmatrix, double **cmatrix);
void smacofSymmetricCtoR(const int n, double **x, double *y);
void smacofFromLowerTriangularRtoC(const int n, const double *rvector,
                                   double **cmatrix);
void smacofAnyCtoR(const int nrow, const int ncol, const double **cmatrix,
                   double *rmatrix);

// smacofUtils.c

void smacofInitRandom(const int n, const int p, double **xini);
void smacofDoubleCenter(const int n, double **delta, double **cross);
double smacofEtaSquareDifference(const int n, const int p, const bool weights,
                                 double **vmat, double **x, double **y);

// smacofMatrixUtils.c

double **smacofMakeAnyMatrix(const int n, const int p);
double **smacofMakeSymmetricMatrix(const int n);
double *smacofMakeVector(const int n);
void smacofCopyAnyMatrix(const int n, const int p, double **x, double **y);
void smacofCopySymmetricMatrix(const int n, double **x, double **y);
void smacofCopyVector(const int n, double *x, double *y);
void smacofFreeMatrix(const int n, double **x);
void smacofFreeVector(double *x);

// smacofLinearAlgebra.c

void smacofJacobi(const int n, const int p, double **a, double **evec,
                  double *eval, const int itmax, const int eps,
                  const bool verbose);
void smacofDoubleJacobi(const int n, double **a, double **b, double **evec,
                        double *eval, const int itmax, const int ieps,
                        const bool verbose);
void smacofInvertPositiveDefiniteMatrix(const int n, double **x, double **xinv);
void smacofGramSchmidt(const int n, const int p, double **x, double **q);
void smacofMultiplyAnyAnyMatrix(const int n, const int p, const int m,
                                double **a, double **x, double **y);
void smacofCrossprodAnyAnyMatrix(const int n, const int p, const int m,
                                 double **x, double **y, double **z);
void smacofDistance(const int n, const int p, double **x, double **d);
void smacofCenter(const int n, const int p, double **x);
void smacofScaleMatrixColumns(const int n, const int m, const double p,
                              double **x, double *y, double **v);
void smacofMultiplyAnySymmetricAnyMatrix(const int n, const int m, double **x,
                                         double **a, double **u);

// smacofPrintRead.c

void smacofPrintAnyMatrix(FILE *stream, const int n, const int p,
                          const int width, const int precision, double **x);
void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                const int precision, double **x);
void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x);
void smacofReadInputFile(FILE *stream, double *delta);
void smacofReadParameterFile(FILE *stream, int *n, int *p, int *itmax,
                             int *haveinit, int *typeinit, int *feps, int *ceps,
                             int *width, int *precision, int *verbose,
                             int *relax, int *ditmax, int *deps, int *dverbose,
                             int *degree, int *ordinal, int *weights,
                             int *haveknots, double *lowend, double *highend,
                             int *anchor, int *ninner, int *percentiles,
                             int *writefile, int *makeplots, int *transform);
void smacofWriteOutputFile(FILE *stream, const int n, const int p,
                           const bool weights, const int width,
                           const int precision, double **delta, double **w,
                           double **dhat, double **xnew, double **dmat,
                           char *iterstring);
void smacofWriteEvalBmat(FILE *stream, const int n, const int width,
                         const int precision, double **bmat, double **vmat);
double smacofMaxWeights(const int n, double **w);

// smacofAccelerate.c

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew);

// smacofBases.c

void smacofBernsteinBasis(const int n, const int m, const double *x,
                          const bool ordinal, double **z);
void smacofCumsumMatrix(const int n, const int m, double **x);
bool smacofCheckIncreasing(double *innerknots, const double lowend,
                           const double highend, const int ninner);
void smacofExtendPartition(const double *innerknots, const int *multiplicities,
                           const int *order, const int *ninner,
                           const double *lowend, const double *highend,
                           double *extended);
void smacofBisect(const double *x, const double *knots, const int *lowindex,
                  const int *highindex, int *index);
void smacofBsplines(const double *x, const double *knots, const int *order,
                    const int *nknots, int *index, double *q);
void smacofBsplineBasis(const double *x, double *knots, const int *order,
                        const int *nknots, const int *nvalues, double *result);
void smacofMakeInnerKnots(const int ninner, const bool percentiles, const int n,
                          const double lowend, const double highend,
                          double *dsort, double *innerknots);
double smacofPercentile(const int n, const double q, double *dsort);
int smacofDoubleComparison(const void *a, const void *b);

// smacofCore.c

void smacofMakeVMatrix(const int n, double **w, double **vmat);
void smacofInverseVMatrix(const int n, double **vmat, double **vinv);
void smacofGuttmanTransform(const int n, const int p, const bool weights,
                            double **delta, double **dmat, double **w,
                            double **vinv, double **xold, double **xnew);
double smacofRho(const int n, const bool weights, double **w, double **dhat,
                 double **dmat);
double smacofEta2(const int n, const bool weights, double **w, double **dmat);
double smacofDeltaEta2(const int n, const bool weights, double **w,
                       double **dhat);
double smacofStress(const int n, const bool weights, double **w, double **dhat,
                    double **dmat);

// smacofNorm.c

void smacofScale(const int n, const int p, const bool weights, double **delta,
                 double **w, double **dmat, double **xold);
void smacofNormDelta(const int n, const bool weights, double **delta,
                     double **w);
void smacofNormWeights(const int n, double **w);

// smacofInitial.c

void smacofInitTorgerson(const int n, const int p, double **delta,
                         double **xold);
void smacofInitMaximumSum(const int n, const int p, const bool weights,
                          double **delta, double **w, double **xold);
void smacofInitial(const int n, const int p, const int haveinit,
                   const bool weights, double **delta, double **w,
                   double **xini);

// smacofQP.c

void smacofDykstra(double *x, double *z, const double *amat, const int *pn,
                   const int *pm, const int *pitmax, const int *peps,
                   const bool *pverbose);

void smacofHildreth(double *x, const double *y, double *lbd, const double *amat,
                    const int *pn, const int *pm, const int *pitmax,
                    const int *peps, const bool *pverbose);

// smacofTransforms.c

void smacofInterval(const int n, double **delta, double **dmat, double **dhat);

// smacofPlot.c

void smacofShepardPlot(const int n, double **delta, double **dhat,
                       double **dmat);
int smacofCoComp(const void *px, const void *py);
void smacofCoSort(double *x, int *k, const int n);
void smacofSqueeze(const int n, const double lowend, const double highend,
                   double *x);

#endif /* SMACOF_H */
