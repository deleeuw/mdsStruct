#include <stdio.h>
#include <stdlib.h>

void fill(size_t n, size_t p, double (*mat)[n]);
void display(size_t n, size_t p, double (*mat)[n]);
void center(size_t n, size_t p, double (*mat)[n]);

int main(void) {
    size_t n = 10, p = 4;
    double(*mat)[n] = malloc((sizeof *mat) * p);
    if (mat == NULL) {
        fprintf(stderr, "MALLOC FAILURE\n");
        exit(EXIT_FAILURE);
    }
    (void)fill(n, p, mat);
    (void)display(n, p, mat);
    (void)center(n, p, mat);
    (void)display(n, p, mat);
    free(mat);
    return EXIT_SUCCESS;
}

void fill(size_t n, size_t p, double (*mat)[n]) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            mat[i][j] = (double)(i + j);
        }
    }
    return;
}

void display(size_t n, size_t p, double (*mat)[n]) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            printf("%6.2f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void center(size_t n, size_t p, double (*mat)[n]) {
    for (size_t j = 0; j < p; j++) {
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) {
            sum += mat[i][j];
        }
        sum /= (double)n;
        for (size_t i = 0; i < n; ++i) {
            mat[i][j] -= sum;
        }
    }
    return;
}
