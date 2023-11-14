#include <stdio.h>
#include <stdlib.h>

int n = 4, p = 2;

int main() {
    int fill[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    int(*mot)[p] = malloc((sizeof *mot) * n);
    int(*mat)[n][p] = malloc(sizeof *mat);
    printf("sizeof *mot %3ld sizeof mot %3ld\n", sizeof *mot, sizeof mot);
    printf("sizeof *mat %3ld sizeof mat %3ld\n", sizeof *mat, sizeof mat);
    int k = 0;
    printf("\n\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            mot[i][j] = fill[k];
            printf("%3d", mot[i][j]);
            k++;
        }
        printf("\n");
    }
    printf("\n\n");
    k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            (*mat)[i][j] = fill[k];
            printf("%3d", (*mat)[i][j]);
            k++;
        }
        printf("\n");
    }
    printf("\n\n");
    k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            *(mot + j)[i] = fill[k];
            printf("%3d", *(mot + j)[i]);
            k++;
        }
        printf("\n");
    }
    printf("\n\n");
    k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            *(*(mot + j) + i) = fill[k];
            printf("%3d", *(*(mot + j) + i));
            k++;
        }
        printf("\n");
    }
    return EXIT_SUCCESS;
}
