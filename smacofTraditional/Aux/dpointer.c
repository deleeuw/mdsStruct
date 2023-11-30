#include <stdlib.h>
#include <stdio.h>

void pr(const int n, const int m, double **x);
double **matrix(const int n, const int m);
void fl(const int n, const int m, double **x);
void fr(const int n, const int m, double **x);

int main(void){
    int n = 5, m = 10;    
    double **x = matrix(n, m);
    (void)fl(n, m, x);    
    (void)pr(n, m, x);
    (void)fr(n, m, x);
    return EXIT_SUCCESS;
}

// make a matrix

double **matrix(const int n, const int m) {
    double **x = (double **)calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        x[i] = (double *)calloc(m, sizeof(double));
    }
    return x;
}

// fill a matrix

void fl(const int n, const int m, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            x[i][j] = (double)(i + j);
        }
    }
}

// print a matrix

void fr(const int n, const int m, double **x) {
    for (int i = 0; i < n; i++) {
        free(x[i]);
    }
    free(x);
}



void pr(const int n, const int m, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%5.2f ", x[i][j]);
        }
        printf("\n");
    }
}
