#include "smacof.h"

int smacofComparison(const void *px, const void *py) {
    double x = ((struct data *)px)->delta;
    double y = ((struct data *)py)->delta;
    return (int)copysign(1.0, x - y);
}

void smacofSort(double *delta, double *weights, int *order, int *irow, int *icol,
                const int *pn) {
    int n = *pn;
    struct data *xi =
        (struct data *)calloc((size_t)n, (size_t)sizeof(struct data));
    for (int i = 1; i <= n; i++) {
        int ii = VINDEX(i);
        xi[ii].order = ii + 1;
        xi[ii].irow = irow[ii];
        xi[ii].icol = icol[ii];
        xi[ii].delta = delta[ii];
        xi[ii].weight = weights[ii];
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(struct data), smacofComparison);
    for (int i = 1; i <= n; i++) {
        int ii = VINDEX(i);
        order[ii] = xi[ii].order;
        irow[ii] = xi[ii].irow;
        icol[ii] = xi[ii].icol;
        delta[ii] = xi[ii].delta;
        weights[ii] = xi[ii].weight;
    }
    free(xi);
}

void smacofTieBlocks(const double *delta, int *block, double *eps,
                     const int *pn) {
    int n = *pn;
    block[0] = 1;
    for (int i = 1; i < n; i++) {
        if (fabs(delta[i] - delta[i - 1]) < *eps) {
            block[i] = block[i - 1];
        } else {
            block[i] = block[i - 1] + 1;
        }
    }
    return;
}

int main() {
    int n = 6;
    double eps = 1e-6;
    int irow[6] = {2, 3, 4, 3, 4, 4};
    int icol[6] = {1, 1, 1, 2, 2, 3};
    int order[6] = {1, 2, 3, 4, 5, 6};
    int block[6] = {0, 0, 0, 0, 0, 0};
    double weights[6] = {1, 2, 3, 4, 5, 6};
    double delta[6] = {3, 1, 5, 5, 6, 1};
    for (int i = 1; i <= n; i++) {
        int ii = VINDEX(i);
        printf("%3d %3d %3d %6.2f %6.2f\n",
            order[ii], irow[ii], icol[ii], delta[ii], weights[ii]);
    }
    printf("\n\n");
    (void)smacofSort(delta, weights, order, irow, icol, &n);
    (void)smacofTieBlocks(delta, block, &eps, &n);
    for (int i = 1; i <= n; i++) {
        int ii = VINDEX(i);
        printf("%3d %3d %3d %6.2f %6.2f %3d\n",
            order[ii], irow[ii], icol[ii], delta[ii], weights[ii], block[ii]);
    }
}
