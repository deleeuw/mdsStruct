#include "../Include/smacof.h"

/*
void smacofFromCtoR(const double *cmatrix, const int *pnrow, const int *pncol,
                    double *rmatrix) {
    int nrow = *pnrow, ncol = *pncol, k = 0;
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
            *(rmatrix + k) = *(cmatrix + (i * ncol + j));
            k++;
        }
    }
    return;
}
*/

void smacofFromAnyRtoC(const double *rmatrix, const size_t nrow,
                       const size_t ncol, double cmatrix[nrow][ncol]) {
    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            cmatrix[i][j] = *(rmatrix + (j * nrow + i));
        }
    }
    return;
}

void smacofFromSymmetricHollowRtoC(const double *rmatrix, const size_t order,
                                   double cmatrix[order][order]) {
    for (size_t j = 0; j < (order - 1); j++) {
        for (size_t i = (j + 1); i < order; i++) {
            double cell =
                rmatrix[(j * order) - (j * (j + 1) / 2) + (i - j) - 1];
            cmatrix[i][j] = cell;
            cmatrix[j][i] = cell;
        }
    }
    for (size_t i = 0; i < order; i++) {
        cmatrix[i][i] = 0.0;
    }
    return;
}

void smacofFromSymmetricRtoC(const double *rmatrix, const size_t order,
                             double cmatrix[order][order]) {
    for (size_t j = 0; j < (order - 1); j++) {
        for (size_t i = (j + 1); i < order; i++) {
            double cell =
                rmatrix[(j * order) - (j * (j + 1) / 2) + (i - j) - 1];
            cmatrix[i][j] = cell;
            cmatrix[j][i] = cell;
        }
    }
    for (size_t i = 0; i < order; i++) {
        cmatrix[i][i] = 0.0;
    }
    return;
}

/*
double rmatrix[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
size_t order = 4;
size_t width = 4, precision = 2;

int main(void) {
    size_t nrow = 4, ncol = 4, size = nrow * ncol;
    double(*cmatrix)[nrow][ncol] =
        calloc((size_t)size, (size_t)sizeof(*cmatrix));
    (void)smacofFromSHRtoC(rmatrix, nrow, *cmatrix);
    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            printf(" %+*.*f ", width, precision, (*cmatrix)[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free(cmatrix);
    return EXIT_SUCCESS;
}
*/
