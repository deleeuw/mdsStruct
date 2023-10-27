#include "smacof.h"

// Computes the Moore-Penrose inverse of an SDCL
// matrix of order n and rank n-1 from only the elements
// below the diagonal

void smacofMPInverseSDCMatrix(const double *w, double *vinv, const int *pn) {
  int n = *pn, m = n * (n - 1) / 2, ik = 0, jk = 0, ij = 0;
  double add = 1.0 / ((double)n), piv = 0.0;
  double *d = (double *)calloc((size_t)n, (size_t)sizeof(double));
  for (int i = 1; i <= m; i++) {
    vinv[VINDEX(i)] = -w[VINDEX(i)];
  }
  for (int i = 1; i <= n; i++) {
    d[VINDEX(i)] = 0.0;
    for (int j = 1; j <= n; j++) {
      if (i == j) {
        continue;
      }
      d[VINDEX(i)] += vinv[PINDEX(i, j, n)];
    }
    d[VINDEX(i)] = -d[VINDEX(i)] + add;
  }
  for (int i = 1; i <= m; i++) {
    vinv[VINDEX(i)] += add;
  }
  for (int k = 1; k <= n; k++) {
    piv = d[VINDEX(k)];
    for (int j = 1; j <= (n - 1); j++) {
      if (j == k) {
        continue;
      }
      jk = PINDEX(j, k, n);
      for (int i = (j + 1); i <= n; i++) {
        if (i == k) {
          continue;
        }
        ik = PINDEX(i, k, n);
        ij = SINDEX(i, j, n);
        vinv[ij] = vinv[ij] - vinv[ik] * vinv[jk] / piv;
      }
    }
    for (int i = 1; i <= n; i++) {
      if (i == k) {
        continue;
      }
      ik = PINDEX(i, k, n);
      d[VINDEX(i)] = d[VINDEX(i)] - vinv[ik] * vinv[ik] / piv;
    }
    for (int i = 1; i <= n; i++) {
      if (i == k) {
        continue;
      }
      ik = PINDEX(i, k, n);
      vinv[ik] = vinv[ik] / piv;
    }
    d[VINDEX(k)] = -1 / piv;
  }
  for (int j = 1; j <= (n - 1); j++) {
    for (int i = (j + 1); i <= n; i++) {
      ij = SINDEX(i, j, n);
      vinv[ij] = vinv[ij] + add;
    }
  }
  for (int i = 1; i <= n; i++) {
    d[VINDEX(i)] = d[VINDEX(i)] + add;
  }
  free(d);
  return;
}

void smacofMultiplySDCMatrix(const double *a, const double *x, double *y,
                             const int *pn, const int *pp) {
  int n = *pn, p = *pp;
  for (int s = 1; s <= p; s++) {
    for (int i = 1; i <= n; i++) {
      double sum = 0.0;
      for (int j = 1; j <= n; j++) {
        if (j == i) {
          continue;
        }
        int ij = PINDEX(i, j, n);
        sum += -a[ij] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
      }
      y[MINDEX(i, s, n)] = sum;
    }
  }
  return;
}
