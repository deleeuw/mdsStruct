#ifndef SMACOF_COMMON_H
#define SMACOF_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define DEBUG false
#define PI (2.0 * asin(1.0))
#define NaN (0.0/0.0)

#define SQUARE(x) ((x) * (x))
#define THIRD(x) ((x) * (x) * (x))
#define FOURTH(x) ((x) * (x) * (x) * (x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KDELTA(i, j) (((i) == (j)) ? 1 : 0)

// smacofCommonUtils.c

void smacofInitRandom(const unsigned n, const unsigned p, double *x);
double smacofMaxConfigurationDifference(const unsigned n, const unsigned p,
                                      const double *x, const double *y);
double smacofMaxDistanceDifference(const unsigned m, const double *x, const double *y);
double smacofRMSDifference(const unsigned n, const unsigned p, const double *x, double *y);
void smacofDoubleCenter(const unsigned n, const double *delta, double *cross);

// smacofCommonPrint.c

void smacofPrintSDCLMatrix(const unsigned n, const double *x, const unsigned width,
                           const unsigned precision);
void smacofPrintAnyMatrix(const unsigned n, const unsigned p, const double *x,
                          const unsigned width, const unsigned precision);
void smacofPrintSHMatrix(const unsigned n, const double *x, const unsigned width,
                         const unsigned precision);
void smacofPrintLTMatrix(const unsigned n, const double *x, const unsigned width,
                         const unsigned precision);
void smacofPrintSLTMatrix(const unsigned n, const double *x, const unsigned width,
                          const unsigned precision);
void smacofPrintSymmetricMatrix(const unsigned n, const double *x, const unsigned width,
                                const unsigned precision);


// smacofCommonAccelerate.c

double smacofRelax(const double *xold, double *xnew, const double echange,
                 const double pchange, const unsigned np, const unsigned itel,
                 const bool relax);

// smacofCommonLinearAlgebra.c

void smacofJacobi(const unsigned n, const unsigned p, double *a, double *evec, double *eval, 
                  const unsigned itmax, const unsigned eps, const bool verbose);
void smacofInvertPDMatrix(const unsigned n, const double *x);
void smacofGramSchmidt(unsigned n, unsigned p, double *x, double *r);
void smacofMultiplySymmetricMatrix(const unsigned n, const unsigned p, const double *a, 
    const double *x, double *y);
void smacofMPInverseSDCLMatrix(const unsigned n, const double *w, double *vinv);
void smacofDistance(const unsigned n, const unsigned p, const double *x, double *d);
void smacofCenter(const unsigned n, const unsigned p, double *x);


// smacofCommonIndices.c 

unsigned VINDEX(const unsigned i);
unsigned MINDEX(const unsigned i, const unsigned j, const unsigned n);
unsigned SINDEX(const unsigned i, const unsigned j, const unsigned n);
unsigned PINDEX(const unsigned i, const unsigned j, const unsigned n);
unsigned TINDEX(const unsigned i, const unsigned j, const unsigned n);
unsigned UINDEX(const unsigned i, const unsigned j, const unsigned n);

// smacofCommonSort.c

unsigned smacofComparison(const void *px, const void *py);
void smacofMakeIIandJJ(const unsigned n, unsigned *ii, unsigned *jj);
void smacofSort(double *delta, double *weights, unsigned *knum, unsigned *irow,
                        unsigned *icol, const unsigned n);
void smacofTieBlocks(const double *x, unsigned *it, double *eps, const unsigned *n);


#endif /* SMACOF_COMMON_H */
