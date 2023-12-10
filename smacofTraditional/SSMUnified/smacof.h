#ifndef SMACOF_H
#define SMACOF_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#define PI (2.0 * asin(1.0))
#define SSIZE 80
#define OSIZE 256

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

void smacofSSMEngine(const int n, const int p, double **delta,
                      double **w, double **xold, double **xnew,
                      double **dmat, double **bmat, double **vmat,
                      double **vinv, const int init, const int itmax,
                      const int ieps1, const int ieps2, const bool verbose,
                      const bool relax, char *iterstring);

// smacofRCTranslation.c

void smacofFromAnyRtoC(const int nrow, const int ncol, const double *rmatrix,
                       double **cmatrix);
void smacofFromSymmetricHollowRtoC(const int n, const double *rmatrix,
                                   double **cmatrix);
void smacofFromSymmetricRtoC(const int n, const double *rmatrix,
                             double **cmatrix);
void smacofFromLowerTriangularRtoC(const int n, const double *rvector,
                                   double **cmatrix);
void smacofFromAnyCtoR(const int nrow, const int ncol, const double **cmatrix,
                       double *rmatrix);
void smacofFromSymmetricCtoR(const int n, double **x, double *y);

// smacofUtils.c

void smacofInitRandom(const int n, const int p, double **xini);
void smacofDoubleCenter(const int n, double **delta, double **cross);
double smacofMaxConfigurationDifference(const int n, const int p, double **x,
                                        double **y);
double smacofMaxDistanceDifference(const int n, double **dold, double **dnew);
double smacofRMSDifference(const int n, const int p, double **x, double **y);
void smacofZeroAnyMatrix(const int n, const int p, double **x);

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
void smacofMultipleAnySymmetricAnyMatrix(const int n, const int m, double **x,
                                         double **a, double **u);

// smacofPrintRead.c

void smacofPrintAnyMatrix(FILE *stream, const int n, const int p,
                           const int width, const int precision, double **x);
void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                 const int precision, double **x);
void smacofPrintAnyVector(FILE *stream, const int n, const int width,
                           const int precision, double *x);
void smacofReadInputFile(const char *fname, double *delta);
void smacofReadParameterFile(const char *fname, int *n, int *p, int *itmax,
                             int *init, int *feps, int *ceps, int *width,
                             int *precision, int *verbose, int *relax,
                             int *interval, int *degree, int *ordinal, int *weights);

// smacofAccelerate.c

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew);

// smacofBases.c

void smacofBernsteinBase(const int n, const int m, const double *x,
                         const bool ordinal, double **z);
void smacofCumsumMatrix(const int n, const int m, double **x);
void smacofCheckIncreasing(const double *, const double *, const double *,
                           const int *, bool *);
void smacofExtendPartition(const double *, const int *, const int *,
                           const int *, const double *, const double *,
                           double *);
void smacofBisect(const double *, const double *, const int *, const int *,
                  int *);
void smacofBsplines(const double *, const double *, const int *, const int *,
                    int *, double *);
void smacofBsplineBasis(const double *, const double *, const int *,
                        const int *, const int *, double *);


void smacofWriteOutputFile(FILE *stream, const int n, const int p,
                           const bool weights, const int width,
                           const int precision, double **delta, double **w,
                           double **xnew, double **dmat, double **bmat,
                           char *iterstring);

void smacofWriteEvalBmat(FILE *stream, const int n, const int width,
                                 const int precision, double **bmat,
                                 double **vmat);

double smacofMaxWeights(const int n, double **w);

// smacofCore.c


void smacofMakeVMatrix(const int n, double **w, double **vmat);
void smacofInverseVMatrix(const int n, double **vmat, double **vinv);
void smacofGuttmanTransform(const int n, const int p, const bool weights,
                            double **delta, double **dmat, double **w,
                            double **vinv, double **xold, double **xnew);
double smacofStress(const int n, const bool weights, double **delta, double **w,
                    double **dmat);

// smacofNorm.c

void smacofScale(const int n, const int p, const bool weights, double **delta,
                         double **w, double **dmat, double **xold);
void smacofNormDelta(const int n, const bool weights, double **delta, double **w);
void smacofNormWeights(const int n, double **w);
void smacofUnweighting(const int n, double wmax, double **delta, double **w, 
    double **dmat, double **dhat);

// smacofInitial.c

void smacofInitTorgerson(const int n, const int p, const bool weights, double **delta, double **w,
                                 double **xold);
void smacofInitMaximumSum(const int n, const int p, const bool weights, double **delta, double **w,
                                 double **xold);
void smacofInitial(const int n, const int p, const int init, const bool weights,
                           double **delta, double **w, double **xini);

// smacofDerivatives.c

void smacofGradient(const double **delta, const double **w,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const unsigned *pn, const unsigned *pp);
void smacofHessian(const int n, const int p, double **delta,
                           double **w, double **dmat, double **bmat,
                           double **vmat, double **x, double ****hessian);

#endif /* SMACOF_H */
