#include <stdio.h>
#include <stdlib.h>

int main(void) {
    int n = 5, p = 3, width = 4;
    unsigned mit[][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {5.0, 6.0, 7.0}, {7.0, 8.0, 9.0}, {9.0, 10.0, 11.0}}; // matrix
    unsigned mut[n][p]; // n x p VAM, cannot be initialized
    unsigned (*mot)[p] = malloc((sizeof *mot) * n); // mot = pointer to an array of n VAM's of length p
    unsigned (*mat)[n][p] = malloc(sizeof *mat); // mat is pointer to n x p VAM
    printf("sizeof  mit %3ld\n", sizeof mit);
    printf("sizeof  mut %3ld\n", sizeof mut);
    printf("sizeof *mot %3ld sizeof mot %3ld\n", sizeof *mot, sizeof mot);
    printf("sizeof *mat %3ld sizeof mat %3ld\n", sizeof *mat, sizeof mat);
    printf("\n\n mit\n\n");
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            printf("%*d", width, mit[i][j]);
        }
        printf("\n");
    }
    printf("\n\n mot\n\n");
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            mot[i][j] = mit[i][j];
            printf("%*d", width, mot[i][j]);
        }
        printf("\n");
    }
    printf("\n\n mat\n\n");
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            (*mat)[i][j] = mit[i][j];
            printf("%*d", width, (*mat)[i][j]);
        }
        printf("\n");
    }
    printf("\n\n mut\n\n");
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < p; j++) {
            mut[i][j] = mit[i][j];
            printf("%*d", width, mut[i][j]);
        }
        printf("\n");
    }
    free(mot);
    free(mat);
    return EXIT_SUCCESS;
}
