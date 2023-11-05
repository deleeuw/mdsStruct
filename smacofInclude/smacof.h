#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG true
#define NaN (nan(""))

// smacofEngine.c

void smacofSSMWEngine(double *delta, double *weights, const int *irow,
                      const int *icol, double *xini, double *xnew, double *dnew,
                      double *bnew, double *psnew, const int *pinit,
                      const int *pn, const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax);
void smacofSSMUEngine(double *delta, const int *irow, const int *icol,
                      double *xini, double *xnew, double *dnew, double *bnew,
                      double *psnew, const int *pinit, const int *pn,
                      const int *pp, const int *pm, int *pitel,
                      const int *pitmax, const int *peps1, const int *peps2,
                      const bool *pverbose, const bool *prelax);

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const double *delta, const double *weights,
                               const double *dold, double *bmat,
                               const int *irow, const int *icol, const int *pn,
                               const int *pm);
void smacofWeightedMakeVMatrix(const double *weights, double *vmat, const int *irow,
                       const int *icol, const int *pn, const int *pm);
void smacofWeightedEtaSquare(const double *weights, const double *dist,
                             const int *pm, double *etasquare);
void smacofWeightedRho(const double *delta, const double *weights,
                       const double *dist, const int *pm, double *rho);
void smacofWeightedGuttman(const double *vinv, const double *bmat,
                           const double *xold, double *xnew, const int *pn,
                           const int *pp);
void smacofWeightedMakeStress(const double *delta, const double *w, const double *d,
                          const int *m, double *stress);

// smacofUnweightedCore.c

void smacofUnweightedMakeBMatrix(const double *delta, const double *dold,
                                 double *bmat, const int *irow, const int *icol,
                                 const int *pn, const int *pm);
void smacofUnweightedGuttman(const double *bmat, const double *xold,
                             double *xnew, const int *pn, const int *pp);
void smacofUnweightedMakeStress(const double *delta, const double *d, const int *m,
                            double *stress);


// smacofWeightedNorm.c

void smacofWeightedScale(const double *delta, const double *weights,
                         double *dold, double *xold, const int *pn,
                         const int *pp, const int *pm);
void smacofWeightedNormDelta(double *delta, const double *weights,
                             const int *pm);
void smacofWeightedNormWeights(double *weights, const int *pm);

// smacofUnweightedNorm.c

void smacofUnweightedScale(const double *delta, 
                         double *dold, double *xold, const int *pn,
                         const int *pp, const int *pm);
void smacofUnweightedNormDelta(double *delta,  const int *pm);

// smacofCommonUtils.c

void smacofInitRandom(double *xini, const int *pn, const int *pp);
void smacofMaxConfigurationDifference(const double *x, const double *y,
                                      const int *pn, const int *pp,
                                      double *maxdiff);
void smacofMaxDistanceDifference(const double *x, const double *y,
                                 const int *pm, double *maxdiff);
void smacofRMSDifference(const double *x, double *y, const int *pn,
                         const int *pp, double *diff);
void smacofMakeIIandJJ(const int *pn, int *ii, int *jj);
void smacofDoubleCenter(const double *delta, double *cross, const int *pn);


// smacofSort.c

int smacofComparison(const void *px, const void *py);
void smacofWeightedSort(double *delta, double *weights, int *knum, int *irow,
                        int *icol, const int *pn);
void smacofTieBlocks(const double *x, int *it, double *eps, const int *n);

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

// smacofUnweightedInitial.c

void smacofUnweightedInitTorgerson(const double *delta, 
                                 const int *irow, const int *icol, double *xold,
                                 const int *pn, const int *pp, const int *pm);
void smacofUnweightedInitMaximumSum(const double *delta, 
                                  const int *irow, const int *icol,
                                  double *xini, const int *pn, const int *pp,
                                  const int *pm);
void smacofUnweightedInitial(const double *delta,
                           const int *irow, const int *icol, double *xini,
                           const int *pinit, const int *pn, const int *pp,
                           const int *pm);

// smacofLinearAlgebra.c

void smacofJacobi(double *a, double *evec, double *eval, const int *pn,
                  const int *pm, const int *itmax, const int *eps,
                  const bool *verbose);
void smacofInvertPDMatrix(const double *x, double *xinv, const int *pn);
void smacofSimultaneousIteration(double *cross, double *xold, const int *pn,
                                 const int *pp, const int *itmax,
                                 const double *eps, const bool *verbose);
void smacofGramSchmidt(double *x, double *r, int *pn, int *pp);
void smacofMultiplySymmetricMatrix(const double *a, const double *x, double *y,
                                   const int *pn, const int *pp);
void smacofMPInverseSDCLMatrix(const double *w, double *vinv, const int *ndim);
void smacofDistance(const double *x, double *d, const int *pn, const int *pp);
void smacofCenter(double *x, const int *np, const int *pp);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const int *pn, const int *pp);
void smacofWeightedHessian(const double *delta, const double *weights,
                           const double *xconf, const double *dmat,
                           const double *bmat, const double *v, const int *pn,
                           const int *pp);

// smacofCommonPrintRead.c

void smacofPrintSDCLMatrix(const double *v, const int *ndim, const int *pw,
                           const int *pr);
void smacofPrintAnyMatrix(const double *x, const int *pn, const int *pp,
                          const int *pw, const int *pr);
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
void smacofWeightedReadInputFile(char *fname, int *irow, int *icol, double *delta,
                                 double *weights);
void smacofUnweightedReadInputFile(char *fname, int *irow, int *icol, double *delta);

// smacofCommonAccelerate.c

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
