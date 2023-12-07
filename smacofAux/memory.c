#include <stdio.h>
#include <stdlib.h>

int main(void) {
  int n = 10, p = 4;
  double **x = (double **)calloc((size_t)n, sizeof(double *));
  double **z = (double **)calloc((size_t)p, sizeof(double *));
  double y[10][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
                     {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
                     {0, 0, 0, 0}, {0, 0, 0, 0}};
  for (int i = 0; i < n; i++) {
    x[i] = (double *)calloc((size_t)p, sizeof(double));
  }
  for (int s = 0; s < p; s++) {
    z[s] = (double *)calloc((size_t)n, sizeof(double));
  }
  printf("x[i][j] (dynamic allocation on heap, row-major)\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      printf("%ld ", (uintptr_t)&x[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  printf("x[i] and &x[i] and sizeof x[i]\n");
  for (int i = 0; i < n; i++) {
    printf("%ld %ld %ld\n", (uintptr_t)x[i], (uintptr_t)&x[i], sizeof x[i]);
  }
  printf("\n\n");
  printf("x\n");
  printf("%ld\n", (uintptr_t)x);
  printf("\n\n");
  printf("*******************************************\n");
  printf("\n\n");
  printf("y[i][j] (static allocation on stack)\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      printf("%ld ", (uintptr_t)&y[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  printf("y[i] and &y[i] and sizeof y[i]\n");
  for (int i = 0; i < n; i++) {
    printf("%ld %ld %ld\n", (uintptr_t)y[i], (uintptr_t)&y[i], sizeof y[i]);
  }
  printf("\n\n");
  printf("y\n");
  printf("%ld\n", (uintptr_t)y);
  printf("\n\n");
  for (int i = 0; i < n; i++) {
    free(x[i]);
  }
  free(x);
  return EXIT_SUCCESS;
}