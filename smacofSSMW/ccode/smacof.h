#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG false
#define NaN (0.0 / 0.0)

// smacofEngine.c

void smacofEngine(double *delta, double *weights, const int *irow,
                  const int *icol, double *xini, double *xnew, double *dini,
                  double *dnew, double *bnew, double *psnew, const int *pinit,
                  const int *pn, const int *pp, const int *pm, int *pitel,
                  const int *pitmax, const int *peps1, const int *peps2,
                  const bool *pverbose, const bool *prelax,
                  const bool *padjust);

// smacofCore.c

void smacofMakeBMatrix(const double *delta, const double *weights,
                       const double *dold, double *bmat, const int *pm);
void smacofMakeVMatrix(const double *weights, double *v, const int *pm);
void smacofEtaSquare(const double *weights, const double *dist, const int *pm,
                     double *etasquare);
void smacofRho(const double *delta, const double *weights, const double *dist,
               const int *pm, double *rho);
void smacofGuttman(const double *vinv, const double *bmat, const double *xold,
                   double *xnew, const int *pn, const int *pp);
void smacofStress(const double *delta, const double *w, const double *d,
                  const int *m, double *stress);
void smacofDistance(const double *x, double *d, const int *irow,
                    const int *icol, const int *pn, const int *pp,
                    const int *pm);

// smacofNorm.c

void smacofScale(const double *delta, const double *weights, double *dold,
                 double *xold, const int *pn, const int *pp, const int *pm);
void smacofCenter(double *x, const int *np, const int *pp);
void smacofNormDelta(double *delta, const double *weights, const int *pm);
void smacofNormWeights(double *weights, const int *pm);

// smacofUtils.c

void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      const int *pn, const int *pp,
                                      double *maxdiff);
void smacofMaxDistanceDifference(const double *x, const double *y,
                                 const int *pm, double *maxdiff);
void smacofRMSDifference(const double *x, double *y, const int *pn,
                         const int *pp, double *diff);
void smacofMakeIIandJJ(const int *pn, int *ii, int *jj);
void smacofDoubleCenter(const double *delta, double *cross, const int *pn);

// smacofSDCLMatrix.c

void smacofMultiplySDCLMatrix(const double *a, const double *x, double *y,
                              const int *pn, const int *pp);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv, const int *ndim);
void smacofAddSDCLDiagonal(const double *a, double *b, const int *pn);

// smacofSort.c

int smacofComparison(const void *px, const void *py);
void smacofSort(double *delta, double *weights, int *knum, int *irow, int *icol,
                const int *pn);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

// smacofInitial.c

void smacofInitRandom(double *xini, const int *np, const int *pp);
void smacofInitTorgerson(const double *delta, const double *weights,
                         const int *irow, const int *icol, double *xold,
                         const int *pn, const int *pp, const int *pm);
void smacofInitMaximumSum(const double *delta, const double *weights,
                          double *xini, const int *pn, const int *pp,
                          const int *pitmax_j, const int *peps_j,
                          const bool *pverbose_j);
void smacofInitial(const double *delta, const double *weights, const int *irow,
                   const int *icol, double *xini, const int *pinit,
                   const int *pn, const int *pp, const int *pm,
                   const bool *padjust);
void smacofDiagonalAdjust(const double *delta, const double *weights, double *x,
                          const int *pn, const int *pp, const int *pitmax,
                          const int *peps, const bool *pverbose);

// smacofElegant.c

void smacofInitElegant(const double *delta, const double *weights,
                       const int *pn, const int *pm, const double *plbd,
                       double *proot, const int *pitmax, const double *peps,
                       const bool *pverbose);
void smacofPerronRoot(double *a, const int *pn, const double *plbd,
                      double *proot, const int *pitmax, const double *peps,
                      const bool *verbose);

// smacofLinearAlgebra.c

void smacofJacobi(double *a, double *evec, double *eval, const int *pn,
                  const int *pm, const int *itmax, const int *eps,
                  const bool *verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const int *pn);
void smacofSimultaneousIteration(double *cross, double *xold, const int *pn,
                                 const int *pp, const int *itmax,
                                 const double *eps, const bool *verbose);
void smacofGramSchmidt(double *x, double *r, int *pn, int *pp);

// smacofDerivatives.c

void smacofGradient(const double *delta, const double *weights,
                    const double *vinv, const double *dold, const double *xold,
                    double *xnew, double *gradient, const int *pn,
                    const int *pp);
void smacofHessian(const double *delta, const double *weights,
                   const double *xconf, const double *dmat, const double *bmat,
                   const double *v, const int *pn, const int *pp);

// smacofPrint.c

void smacofPrintSDCLMatrix(const double *v, const int *ndim, const int *pw,
                           const int *pr);
void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp,
                          const int *pw, const int *pr);
void smacofPrintSHMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr);
void smacofPrintLTMatrix(const double *d, const int *pn, const int *pw,
                         const int *pr);
void smacofPrintSLTMatrix(const double *d, const int *pn, const int *pw,
                          const int *pr);
void smacofPrintSymmetricMatrix(const double *x, const int *pn, const int *pw,
                                const int *pr);

// smacofAccelerate.c

void smacofRelax(const double *xold, double *xnew, const double *pechange,
                 const double *ppchange, const int *pnp, const int *pitel,
                 const bool *prelax, double *prate);

// structures

struct data {
    int order;
    int irow;
    int icol;
    double delta;
    double weight;
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

static inline int UINDEX(const int i, const int j, const int n) {
    if (i >= j) {
        return (TINDEX(i, j, n));
    } else {
        return (TINDEX(j, i, n));
    }
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
