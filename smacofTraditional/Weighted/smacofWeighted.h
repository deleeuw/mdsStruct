#ifndef SMACOF_WEIGHTED_H
#define SMACOF_WEIGHTED_H

#include "../Common/smacofCommon.h"

void smacofSSMWEngine(const int n, const int p, double **delta,
                      double **weights, double **xold, double **xnew,
                      double **dmat, double **bmat, double **vmat,
                      double **vinv, const int init, const int itmax,
                      const int ieps1, const int ieps2, const bool verbose,
                      const bool relax, char *iterstring);

void smacofWeightedWriteOutputFile(FILE *stream, const int n, const int p,
                                   const int width, const int precision,
                                   double **delta, double **weights,
                                   double **xnew, double **dmat, double **bmat,
                                   char *iterstring);

void smacofWeightedWriteEvalBmat(FILE *stream, const int n, const int width,
                                 const int precision, double **bmat,
                                 double **vmat);

double smacofWeightedMaxWeights(const int n, double **weights);

// smacofWeightedCore.c

void smacofWeightedMakeBMatrix(const int n, double **delta, double **weights,
                               double **dmat, double **bmat);
void smacofWeightedMakeVMatrix(const int n, double **weights, double **vmat);
void smacofWeightedInverseVMatrix(const int n, double **vmat, double **vinv);
void smacofWeightedGuttman(const int n, const int p, double **vinv,
                           double **bmat, double **xold, double **xnew);
double smacofWeightedMakeStress(const int n, double **delta, double **weights,
                                double **dist);

// smacofWeightedNorm.c

void smacofWeightedScale(const int n, const int p, double **delta,
                         double **weights, double **dmat, double **xold);
void smacofWeightedNormDelta(const int n, double **delta, double **weights);
void smacofWeightedNormWeights(const int n, double **weights);
void smacofUnweighting(const int n, double wmax, double **delta, double **weights, 
    double **dmat, double **dhat);

// smacofWeightedInitial.c

void smacofWeightedInitTorgerson(const int n, const int p, double **delta,
                                 double **xold);
void smacofWeightedInitMaximumSum(const int n, const int p, double **delta,
                                  double **weights, double **xold);
void smacofWeightedInitial(const int n, const int p, const int init,
                           double **delta, double **weights, double **xini);

// smacofWeightedDerivatives.c

void smacofWeightedGradient(const double *delta, const double *weights,
                            const double *vinv, const double *dold,
                            const double *xold, double *xnew, double *gradient,
                            const unsigned *pn, const unsigned *pp);
void smacofWeightedHessian(const int n, const int p, double **delta,
                           double **weights, double **dmat, double **bmat,
                           double **vmat, double **x, double ****hessian);

#endif /* SMACOF_WEIGHTED_H */
