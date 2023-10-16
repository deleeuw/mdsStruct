#ifndef SMACOF_H
#define SMACOF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define DEBUG false

// smacofEngine.c

void smacofEngine(double *delta, double *weights, double * xold, const int *ii,
                  const int *jj, const int *pm, const int *pn, const int *pp,
                  const int *itmax, const int *peps1, const int *peps2, 
                  const bool *verbose);

// smacofCore.c

void smacofGuttman(const double *delta, const double *weights, const double *vinv,
                   const double *dist, double *xold, double *xnew, const int *pn,
                   const int *pp);
void smacofLoss(const double *delta, const double *w, const double *d,
                const int *m, double *loss);
void smacofDist(const double *x, double *d, const int *ii, const int *jj,
                const int *pm, const int *pn, const int *pp);

// smacofNorm.c

void smacofScale(const double *delta, const double *weights, double *dold,
                 double *xold, const int *pm, const int *pn, const int *pp);
void smacofCenter(double *x, const int *np, const int *pp);
void smacofNormDelta(double *delta, const double *weights, const int *pm);
void smacofNormWeights(double *weights, const int *pm); 

// smacofUtils.c

void smacofPrintSDCMatrix(const double *v, const int *ndim, const int *width, const int *precision);
void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp, const int *pw, const int *pr);
void smacofPrintSymmetricHollowMatrix(const double *d, const int *pn, const int *pw,
      const int *pr);
void smacofMultiplySDCMatrix(const double *a, const double *x, double *y, const int *pn,
      const int *pp);
void smacofMPInverseSDCMatrix(const double *w, double *vinv, const int *ndim);
void smacofRootMeanSquare(const double *x, const double *y, const double *w, const int *pn,
      const int *pp, double *change);
void smacofMaxDifference(const double *x, const double *y, double *maxdiff, const int *pn, const int *pp);

// smacofSort.c

int  smacofComparison(const void *px, const void *py);
void smacofSort(double *delta, double *weight, int *row, int *col, int *index,
                const int *n);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

// smacofInitial.c

void smacofDoubleCenter(const double *delta, double *cross, const int *pn);
void smacofInitial(const double *delta, double *xold, const int *pn, const int *pp);
void smacofGramSchmidt(double *x, double *r, int *pn, int *pp);
void smacofSimultaneousIteration(double *cross, double *xold, const int *pn, const int *pp,
                                 const int *itmax, const double *eps);
// smacofDerivatives.c


struct fiveTuple {
    int index;
    int row;
    int col;
    double delta;
    double weight;
};

static inline int VINDEX(const int i);
static inline int SINDEX(const int i, const int j, const int n);
static inline int MINDEX(const int i, const int j, const int n);
static inline int PINDEX(const int i, const int j, const int n);

static inline double SQUARE(const double);
static inline double MAX(const double, const double);
static inline double MIN(const double, const double);
static inline int IMIN(const int, const int);
static inline int IMAX(const int, const int);


// VINDEX takes 1,...,n to 0,...,n-1

static inline int VINDEX(const int i) {
 return (i - 1); 
}

// MINDEX retrieves element (i,j) from an n x m matrix in
// column-major-order storage

static inline int MINDEX(const int i, const int j, const int n) {
    return ((i - 1) + (j - 1) * n);
}

// SINDEX retrieves element (i, j) from a strict lower triangular matrix
// of order n. Thus always i > j.

static inline int SINDEX(const int i, const int j, const int n) {
    if (i <= j) {
        printf("SINDEX error\n");
        return EXIT_FAILURE;
    }
    return (((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1);
}

// PINDEX retrieves element (max(i, j), min(i,j)) from a strict lower
// triangularmatrix of order n. Thus always i ≠ j.

static inline int PINDEX(const int i, const int j, const int n) {
    if (i == j) {
        printf("PINDEX error\n");
        return EXIT_FAILURE;
    }
    if (i > j) {
        return (SINDEX(i, j, n));
    } else {
        return (SINDEX(j, i, n));
    }
}

static inline double SQUARE(const double x) { return (x * x); }

static inline double MAX(const double x, const double y) {
    return ((x > y) ? x : y);
}

static inline double MIN(const double x, const double y) {
    return ((x < y) ? x : y);
}

static inline int IMAX(const int x, const int y) {
 return ((x > y) ? x : y); 
}

static inline int IMIN(const int x, const int y) {
 return ((x < y) ? x : y); 
}

#endif /* SMACOF_H */

