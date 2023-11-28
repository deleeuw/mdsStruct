#include "../Include/smacof.h"

void smacofFromAnyCtoR(const unsigned nrow, const unsigned ncol,
                       const double (*cmatrix)[nrow][ncol], double *rmatrix) {
    unsigned k = 0;
    for (unsigned j = 0; j < ncol; j++) {
        for (unsigned i = 0; i < nrow; i++) {
            *(rmatrix + k) = (*cmatrix)[i][j];
            k++;
        }
    }
    return;
}

void smacofFromAnyRtoC(const unsigned nrow, const unsigned ncol,
                       const double *rmatrix, double (*cmatrix)[nrow][ncol]) {
    for (unsigned i = 0; i < nrow; i++) {
        for (unsigned j = 0; j < ncol; j++) {
            (*cmatrix)[i][j] = *(rmatrix + (j * nrow + i));
        }
    }
    return;
}

void smacofFromSymmetricHollowRtoC(const unsigned n, const double *rmatrix,
                                   double (*cmatrix)[n][n]) {
    for (unsigned j = 0; j < (n - 1); j++) {
        for (unsigned i = (j + 1); i < n; i++) {
            double cell = rmatrix[(j * n) - (j * (j + 1) / 2) + (i - j) - 1];
            (*cmatrix)[i][j] = cell;
            (*cmatrix)[j][i] = cell;
        }
    }
    for (unsigned i = 0; i < n; i++) {
        (*cmatrix)[i][i] = 0.0;
    }
    return;
}

void smacofFromSymmetricRtoC(const unsigned n, const double *rmatrix,
                             double (*cmatrix)[n][n]) {
    for (unsigned j = 0; j < n; j++) {
        for (unsigned i = j; i < n; i++) {
            double cell = rmatrix[(j * n) - (j * (j - 1) / 2) + (i - j)];
            (*cmatrix)[i][j] = cell;
            (*cmatrix)[j][i] = cell;
        }
    }
    return;
}

/*
double rshmatrix[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
double rsymatrix[10] = {10.0, 1.0, 2.0, 3.0, 10.0, 4.0, 5.0, 10.0, 6.0, 10.0};
double raymatrix[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
unsigned n = 4, nrow = 4, ncol = 2;
unsigned width = 6, precision = 2;

int main(void) {
    double(*cshmatrix)[nrow] = malloc(sizeof(*cshmatrix) * ncol);
    assert(!(cshmatrix == NULL));
    double(*csymatrix)[nrow] = malloc(sizeof(*csymatrix) * ncol);
    assert(!(cshmatrix == NULL));
    double(*caymatrix)[nrow] = malloc(sizeof(*csymatrix) * ncol);
    assert(!(cshmatrix == NULL));
    (void)smacofFromSymmetricHollowRtoC(n, rshmatrix, cshmatrix);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            printf(" %+*.*f ", width, precision, cshmatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free(cshmatrix);
    (void)smacofFromSymmetricRtoC(n, rsymatrix, csymatrix);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            printf(" %+*.*f ", width, precision, csymatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free(csymatrix);
    (void)smacofFromAnyRtoC(nrow, ncol, raymatrix, caymatrix);
    for (unsigned i = 0; i < nrow; i++) {
        for (unsigned j = 0; j < ncol; j++) {
            printf(" %+*.*f ", width, precision, caymatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    (void)smacofFromCtoR(nrow, ncol, caymatrix, raymatrix);
    for (unsigned i = 0; i < nrow; i++) {
        for (unsigned j = 0; j < ncol; j++) {
            printf(" %+*.*f ", width, precision, raymatrix[j * nrow + i]);
        }
        printf("\n");
    }
    printf("\n\n");
    free(caymatrix);
    return EXIT_SUCCESS;
}
*/
