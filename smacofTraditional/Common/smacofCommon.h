#ifndef SMACOF_COMMON_H
#define SMACOF_COMMON_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PI (2.0 * asin(1.0))

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofCommonRCTranslation.c

void smacofFromAnyRtoC(const int nrow, const int ncol, const double *rmatrix,
                       double **cmatrix);
void smacofFromSymmetricHollowRtoC(const int n, const double *rmatrix,
                                   double **cmatrix);
void smacofFromSymmetricRtoC(const int n, const double *rmatrix,
                             double **cmatrix);
void smacofFromAnyCtoR(const int nrow, const int ncol, const double **cmatrix,
                       double *rmatrix);

// smacofCommonUtils.c

void smacofInitRandom(const int n, const int p, double **xini);
void smacofDoubleCenter(const int n, double **delta, double **cross);
double smacofMaxConfigurationDifference(const int n, const int p, double **x,
                                        double **y);
double smacofMaxDistanceDifference(const int n, double **dold, double **dnew);
double smacofRMSDifference(const int n, const int p, double **x, double **y);
void smacofZeroAnyMatrix(const int n, const int p, double **x);

// smacofCommonMatrixUtils.c

double **smacofMakeAnyMatrix(const int n, const int p);
double *smacofMakeAnyVector(const int n);
double **smacofMakeLowerDiagonalMatrix(const int n);
void smacofCopyAnyMatrix(const int n, const int p, double **x, double **y);
void smacofCopyAnyVector(const int n, double *x, double *y);
void smacofFreeAnyMatrix(const int n, double **x);
void smacofFreeAnyVector(double *x);

// smacofCommonLinearAlgebra.c

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

// smacofCommonPrintRead.c

void smacofPrintAnyMatrix(const int n, const int p, const int width,
                          const int precision, double **x);
void smacofFPrintAnyMatrix(FILE *stream, const int n, const int p,
                           const int width, const int precision, double **x);
void smacofPrintSymmetricMatrix(const int n, const int width,
                                const int precision, double **x);
void smacofPrintFSymmetricMatrix(FILE *stream, const int n, const int width,
                                 const int precision, double **x);
void smacofPrintAnyVector(const int n, const int width, const int precision,
                          double *x);
void smacofFPrintAnyVector(FILE *stream, const int n, const int width,
                           const int precision, double *x);
void smacofReadInputFile(const char *fname, double *delta);
void smacofReadParameterFile(const char *fname, int *n, int *p, int *itmax,
                             int *init, int *feps, int *ceps, int *width,
                             int *precision, int *verbose, int *relax);
void smacofUnweightedInterval(const int n, const double **delta,
                              const double **dmat, double **dhat);

// smacofCommonAccelerate.c

void smacofRelax(const int n, const int p, const double rate, double **xold,
                 double **xnew);

// smacofCommonBases.c

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

#endif /* SMACOF_COMMON_H */
