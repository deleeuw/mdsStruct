#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG false

// smacofEngine.c

void smacofEngine(double *delta, double *weights, double *xini, double *xnew,
                  double *dini, double *dnew, double *bmat, double *psnew,
                  const int *pinit, const int *pn, const int *pp, int *pitel,
                  const int *pitmax, const int *peps1, const int *peps2,
                  const bool *pverbose);

// smacofCore.c

void smacofMakeBMatrix(const double *delta, const double *weights,
                       const double *dold, double *bmat, const int *pm);
void smacofMakeVMatrix(const double *weights, double *v, const int *pn);
void smacofEtaSquare(const double *weights, const double *dist, const int *pm,
                     double *etasquare);
void smacofRho(const double *delta, const double *weights, const double *dist,
               const int *pm, double *rho);

void smacofGuttman(const double *vinv, const double *bmat, const double *xold,
                   double *xnew, const int *pn, const int *pp);
void smacofStress(const double *delta, const double *w, const double *d,
                  const int *m, double *stress);
void smacofDistance(const double *x, double *d, const int *pn, const int *pp);

// smacofNorm.c

void smacofScale(const double *delta, const double *weights, double *dold,
                 double *xold, const int *pn, const int *pp);
void smacofCenter(double *x, const int *np, const int *pp);
void smacofNormDelta(double *delta, const double *weights, const int *pm);
void smacofNormWeights(double *weights, const int *pm);

// smacofUtils.c

void smacofMultiplySDCMatrix(const double *a, const double *x, double *y,
                             const int *pn, const int *pp);
void smacofMPInverseSDCMatrix(const double *w, double *vinv, const int *ndim);
void smacofVChange(const double *x, const double *y, const double *v,
                   const int *pn, const int *pp, double *change);
void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      double *maxdiff, const int *pn,
                                      const int *pp);
void smacofMaxDistanceDifference(const double *x, const double *y,
                                 double *maxdiff, const int *pm);
void smacofAddSDCDiagonal(const double *a, double *b, const int *pn);
void smacofMakeIIandJJ(const int *pn, int *ii, int *jj);

// smacofSort.c

int smacofComparison(const void *px, const void *py);
void smacofSort(double *delta, double *weight, int *row, int *col, int *index,
                const int *n);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

// smacofInitial.c

void smacofInitRandom(double *xini, const int *np, const int *pp);
void smacofInitSDCmatrix(const double *delta, const double *weights,
                         double *xini, const int *pn, const int *pp,
                         const int *pitmax_j, const int *peps_j,
                         const bool *pverbose_j);
void smacofInitial(const double *delta, const double *weights, double *xini,
                   double *dini, const int *pinit, const int *pn,
                   const int *pp);

// smacofTorgerson.c

void smacofDoubleCenter(const double *delta, double *cross, const int *pn);
void smacofTorgerson(const double *delta, double *xold, const int *pn,
                     const int *pp, const int *pitmax, const int *peps,
                     const bool *pverbose);

// smacofElegant.c

void smacofElegant(const double *delta, const double *weights, const int *pn,
                   const int *pm, const double *plbd, double *proot,
                   const int *pitmax, const double *peps, const bool *pverbose);
void smacofPerronRoot(double *a, const int *pn, const double *plbd,
                      double *proot, const int *pitmax, const double *peps,
                      const bool *verbose);

// smacofLinearAlgebra.c

void smacofGramSchmidt(double *x, double *r, int *pn, int *pp);
void smacofSimultaneousIteration(double *cross, double *xold, const int *pn,
                                 const int *pp, const int *itmax,
                                 const double *eps, const bool *verbose);
void smacofJacobi(double *a, double *evec, double *eval, const int *pn,
                  const int *itmax, const int *eps, const bool *verbose);

// smacofDerivatives.c

void smacofGradient(const double *delta, const double *weights,
                    const double *vinv, const double *dold, const double *xold,
                    double *xnew, double *gradient, const int *pn,
                    const int *pp);
void smacofHessian(const double *delta, const double *weights,
                   const double *xconf, const double *dmat, const double *bmat,
                   const double *v, const int *pn, const int *pp);

// smacofPrint.c

void smacofPrintSDCMatrix(const double *v, const int *ndim, const int *width,
                          const int *precision);
void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp,
                          const int *pw, const int *pr);
void smacofPrintSHMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr);
void smacofPrintLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr);
void smacofPrintSLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr);

// structures

struct fiveTuple {
    int index;
    int row;
    int col;
    double delta;
    double weight;
};

struct allEpsilon {
    int stressEps;
    int confEps;
    int distEps;
    int jacobiEps;
    int blockEps;
};

struct allItmax {
    int mainItmax;
    int jacobiItmax;
};

// VINDEX takes 1,...,n to 0,...,n-1

static inline int VINDEX(const int i) { return (i - 1); }

// MINDEX retrieves element (i,j) from an n x m matrix in
// column-major-order storage

static inline int MINDEX(const int i, const int j, const int n) {
    return ((i - 1) + (j - 1) * n);
}

// SINDEX retrieves element (i, j) from a strictly lower triangular matrix
// of order n. Thus always i > j.

static inline int SINDEX(const int i, const int j, const int n) {
    if (i <= j) {
        printf("SINDEX error\n");
        return EXIT_FAILURE;
    }
    return (((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1);
}

// PINDEX retrieves element (max(i, j), min(i,j)) from a strictly lower
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

// TINDEX retrieves element (i, j) from a lower triangular matrix
// of order n. Thus always i >= j.

static inline int TINDEX(const int i, const int j, const int n) {
    if (i < j) {
        printf("TINDEX error\n");
        return EXIT_FAILURE;
    }
    return ((j - 1) * n) - ((j - 1) * (j - 2) / 2) + (i - (j - 1)) - 1;
}

static inline double SQUARE(const double x) { return (x * x); }

static inline double THIRD(const double x) { return (x * x * x); }

static inline double MAX(const double x, const double y) {
    return ((x > y) ? x : y);
}

static inline double MIN(const double x, const double y) {
    return ((x < y) ? x : y);
}

static inline int IMAX(const int x, const int y) { return ((x > y) ? x : y); }

static inline int IMIN(const int x, const int y) { return ((x < y) ? x : y); }

static inline int KDELTA(const int i, const int j) {
    return ((i == j) ? 1 : 0);
}

#endif /* SMACOF_H */
