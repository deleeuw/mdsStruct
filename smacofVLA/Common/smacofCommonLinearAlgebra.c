#include "../Include/smacof.h"

void smacofJacobi(const unsigned n, const unsigned ndim, double (*a)[n][n],
                  double (*evec)[n][n], double (*eval)[n], const unsigned itmax,
                  const unsigned ieps, const bool verbose) {
    unsigned itel = 1;
    double d = 0.0, s = 0.0, t = 0.0, u = 0.0, v = 0.0, p = 0.0, q = 0.0,
           r = 0.0;
    double fold = 0.0, fnew = 0.0, eps = pow(10.0, -(double)ieps);
    double(*oldi)[n] = malloc(sizeof *oldi);
    double(*oldj)[n] = malloc(sizeof *oldj);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            (*evec)[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (unsigned i = 0; i < ndim; i++) {
        fold += SQUARE((*a)[i][i]);
    }
    while (true) {
        for (unsigned j = 0; j < ndim; j++) {
            for (unsigned i = j + 1; i < n; i++) {
                p = (*a)[i][j];
                q = (*a)[i][i];
                r = (*a)[j][j];
                if (fabs(p) < 1e-10) {
                    continue;
                }
                d = (q - r) / 2.0;
                s = (p < 0) ? -1.0 : 1.0;
                t = -d / sqrt(SQUARE(d) + SQUARE(p));
                u = sqrt((1 + t) / 2);
                v = s * sqrt((1 - t) / 2);
                for (unsigned k = 0; k < n; k++) {
                    unsigned ik = MIN(i, k);
                    unsigned ki = MAX(i, k);
                    unsigned jk = MIN(j, k);
                    unsigned kj = MAX(j, k);
                    (*oldi)[k] = (*a)[ki][ik];
                    (*oldj)[k] = (*a)[kj][jk];
                }
                for (unsigned k = 0; k < n; k++) {
                    unsigned ik = MIN(i, k);
                    unsigned ki = MAX(i, k);
                    unsigned jk = MIN(j, k);
                    unsigned kj = MAX(j, k);
                    (*a)[ki][ik] = u * (*oldi)[k] - v * (*oldj)[k];
                    (*a)[kj][jk] = v * (*oldi)[k] + u * (*oldj)[k];
                }
                for (unsigned k = 0; k < n; k++) {
                    (*oldi)[k] = (*evec)[k][i];
                    (*oldj)[k] = (*evec)[k][j];
                    (*evec)[k][i] = u * (*oldi)[k] - v * (*oldj)[k];
                    (*evec)[k][j] = v * (*oldi)[k] + u * (*oldj)[k];
                }
                (*a)[i][i] = SQUARE(u) * q + SQUARE(v) * r - 2 * u * v * p;
                (*a)[j][j] = SQUARE(v) * q + SQUARE(u) * r + 2 * u * v * p;
                (*a)[i][j] = u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
                (*a)[j][i] = u * v * (q - r) + (SQUARE(u) - SQUARE(v)) * p;
            }
        }
        fnew = 0.0;
        for (unsigned i = 0; i < ndim; i++) {
            fnew += SQUARE((*a)[i][i]);
        }
        if (verbose) {
            printf("itel %3d fold %15.10f fnew %15.10f\n", itel, fold, fnew);
        }
        if (((fnew - fold) < eps) || (itel == itmax)) {
            break;
        }
        fold = fnew;
        itel++;
    }
    for (unsigned i = 0; i < n; i++) {
        (*eval)[i] = (*a)[i][i];
    }
    free(oldi);
    free(oldj);
    return;
}

void smacofGramSchmidt(const unsigned n, const unsigned p, double (*x)[n][p],
                       double (*r)[n]) {
    unsigned s = 0;
    while (s < p) {
        for (unsigned t = 0; t < s; t++) {
            double sum = 0.0;
            for (unsigned i = 1; i < n; i++) {
                sum += (*x)[i][t] * (*x)[i][s];
            }
            for (unsigned i = 1; i < n; i++) {
                (*x)[i][s] -= sum * (*x)[i][t];
            }
        }
        double sum = 0.0;
        for (unsigned i = 0; i < n; i++) {
            sum += SQUARE((*x)[i][s]);
        }
        sum = sqrt(sum);
        (*r)[s] = sum;  // update all of r
        for (unsigned i = 1; i <= n; i++) {
            (*x)[i][s] /= sum;
        }
        s++;
    }
    return;
}

void smacofCenter(const unsigned n, const unsigned p, double (*x)[n][p]) {
    for (unsigned s = 0; s < p; s++) {
        double sum = 0.0;
        for (unsigned i = 0; i < n; i++) {
            sum += (*x)[i][s];
        }
        sum /= (double)n;
        for (unsigned i = 0; i < n; i++) {
            (*x)[i][s] -= sum;
        }
    }
    return;
}

/*
void smacofInvertPDMatrix(const double *x, double *xinv, const unsigned *pn) {
    unsigned n = *pn, m = n * (n + 1) / 2, ik = 0, jk = 0, ij = 0;
    for (unsigned k = 1; k <= m; k++) {
        xinv[VINDEX(k)] = x[VINDEX(k)];
    }
    for (unsigned k = 1; k <= n; k++) {
        double piv = xinv[TINDEX(k, k, n)];
        for (unsigned j = 1; j <= n; j++) {
            if (j == k) {
                continue;
            }
            jk = UINDEX(j, k, n);
            for (unsigned i = j; i <= n; i++) {
                if (i == k) {
                    continue;
                }
                ik = UINDEX(i, k, n);
                ij = TINDEX(i, j, n);
                xinv[ij] = xinv[ij] - xinv[ik] * xinv[jk] / piv;
            }
        }
        for (unsigned i = 1; i <= n; i++) {
            if (i == k) {
                continue;
            }
            ik = UINDEX(i, k, n);
            xinv[ik] = xinv[ik] / piv;
        }
        xinv[TINDEX(k, k, n)] = -1 / piv;
    }
    for (unsigned k = 1; k <= m; k++) {
        xinv[VINDEX(k)] = -xinv[VINDEX(k)];
    }
    return;
}

void smacofMPInverseSDCLMatrix(const double *vmat, double *vinv,
                               const unsigned *pn) {
    unsigned n = *pn, nn = n * (n + 1) / 2;
    double add = 1.0 / ((double)n);
    double *vadd = malloc((sizeof *vadd) * nn);
    for (unsigned j = 1; j <= n; j++) {
        for (unsigned i = j; i <= n; i++) {
            vadd[TINDEX(i, j, n)] = vmat[TINDEX(i, j, n)] + add;
        }
    }
    (void)smacofInvertPDMatrix(vadd, vinv, pn);
    for (unsigned k = 1; k <= nn; k++) {
        vinv[VINDEX(k)] -= add;
    }
    free(vadd);
    return;
}
*/

void smacofMultiplySymmetricMatrix(const unsigned n, const unsigned p,
                                   const double (*a)[n][n],
                                   const double (*x)[n][p], double (*y)[n][p]) {
    for (unsigned s = 0; s < p; s++) {
        for (unsigned i = 0; i < n; i++) {
            double sum = 0.0;
            for (unsigned j = 0; j < n; j++) {
                sum += (*a)[i][j] * (*x)[j][s];
            }
            (*y)[i][s] = sum;
        }
    }
    return;
}

void smacofDistance(const unsigned n, const unsigned p, const double (*x)[n][p],
                    double (*d)[n][n]) {
    for (unsigned j = 0; j < n; j++) {
        for (unsigned i = 0; i < n; i++) {
            double sum = 0.0;
            for (unsigned s = 0; s < p; s++) {
                sum += SQUARE((*x)[i][s] - (*x)[j][s]);
            }
            (*d)[i][j] = sqrt(fabs(sum));
        }
    }
    return;
}

/*int main(void) {
  unsigned n = 3;
  unsigned p = 2;
  unsigned width = 6;
  unsigned precision = 2;
  double rr[6] = {6.0, -4.0, -2.0, 6.0, -2.0, 4.0};
  double rx[6] = {1.0, 2.0, 3.0, 3.0, 2.0, 1.0};
  double(*cr)[n][n] = malloc(sizeof(*cr));
  assert(!(cr == NULL));
  double(*cx)[n][p] = malloc(sizeof(*cx));
  assert(!(cx == NULL));
  double(*cy)[n][p] = malloc(sizeof(*cy));
  assert(!(cy == NULL));
  (void)smacofFromSymmetricRtoC(n, rr, cr);
  (void)smacofPrintAnyMatrix(n, n, width, precision, cr);
  (void)smacofFromAnyRtoC(n, p, rx, cx);
  (void)smacofPrintAnyMatrix(n, p, width, precision, cx);
  (void)smacofMultiplySymmetricMatrix(n, p, cr, cx, cy);
  (void)smacofPrintAnyMatrix(n, p, width, precision, cy);
  return EXIT_SUCCESS;
}*/
