#include "smacof.h"

void smacofPrintAnyMatrix(FILE *stream, const int n, const int p,
                          const int width, const int precision, double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            fprintf(stream, "%+*.*f", width, precision, x[i][s]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    return;
}

void smacofPrintSymmetricMatrix(FILE *stream, const int n, const int width,
                                const int precision, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            fprintf(stream, "%+*.*f", width, precision, x[i][j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n\n");
    return;
}

void smacofPrintVector(FILE *stream, const int n, const int width,
                       const int precision, double *x) {
    for (int i = 0; i < n; i++) {
        fprintf(stream, "%+*.*f", width, precision, x[i]);
    }
    fprintf(stream, "\n\n");
    return;
}

// read data from a file into a one-dimensional array

void smacofReadInputFile(FILE *stream, double *delta) {
    int k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    return;
}

void smacofReadParameterFile(FILE *stream, int *n, int *p, int *itmax,
                             int *init, int *feps, int *ceps, int *width,
                             int *precision, int *verbose, int *relax,
                             int *interval, int *degree, int *ordinal,
                             int *weights) {
    char *line = (char *)malloc(80 * sizeof(char));
    char aux1[20], aux2[20];
    int par = 0;
    size_t len = 256;
    ssize_t read;
    while ((read = getline(&line, &len, stream)) != -1) {
        sscanf(line, "%s %s %d", aux1, aux2, &par);
        if (strcmp("n", aux1) == 0) {
            *n = par;
        }
        if (strcmp("p", aux1) == 0) {
            *p = par;
        }
        if (strcmp("itmax", aux1) == 0) {
            *itmax = par;
        }
        if (strcmp("init", aux1) == 0) {
            *init = par;
        }
        if (strcmp("feps", aux1) == 0) {
            *feps = par;
        }
        if (strcmp("ceps", aux1) == 0) {
            *ceps = par;
        }
        if (strcmp("width", aux1) == 0) {
            *width = par;
        }
        if (strcmp("precision", aux1) == 0) {
            *precision = par;
        }
        if (strcmp("verbose", aux1) == 0) {
            *verbose = par;
        }
        if (strcmp("relax", aux1) == 0) {
            *relax = par;
        }
        if (strcmp("interval", aux1) == 0) {
            *interval = par;
        }
        if (strcmp("degree", aux1) == 0) {
            *degree = par;
        }
        if (strcmp("ordinal", aux1) == 0) {
            *ordinal = par;
        }
        if (strcmp("weights", aux1) == 0) {
            *weights = par;
        }
    }
    free(line);
    return;
}

void smacofWriteEvalBmat(FILE *stream, const int n, const int width,
                         const int precision, double **bmat, double **vmat) {
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeVector(n);
    double **bvec = smacofMakeAnyMatrix(n, n);
    (void)smacofDoubleJacobi(n, bmat, vmat, evec, eval, 100, 15, false);
    fprintf(stream, "Eigenvalues V^+B(X)\n\n");
    (void)smacofPrintVector(stream, n, width, precision, eval);
    (void)smacofFreeVector(eval);
    (void)smacofFreeMatrix(n, evec);
    (void)smacofFreeMatrix(n, bvec);
}

void smacofWriteOutputFile(FILE *stream, const int n, const int p,
                           const bool weights, const int width,
                           const int precision, double **delta, double **w,
                           double **dhat, double **xnew, double **dmat,
                           double **bmat, char *iterstring) {
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    fprintf(stream, "%s\n\n", asctime(tm));
    fprintf(stream, "Final Iteration\n\n");
    fprintf(stream, "%s\n\n", iterstring);
    fprintf(stream, "Data Delta\n\n");
    smacofPrintSymmetricMatrix(stream, n, width, precision, delta);
    if (weights) {
        fprintf(stream, "Weights\n\n");
        smacofPrintSymmetricMatrix(stream, n, width, precision, w);
    }
    fprintf(stream, "Transformed Data Dhat\n\n");
    smacofPrintSymmetricMatrix(stream, n, width, precision, dhat);
    fprintf(stream, "Configuration X\n\n");
    smacofPrintAnyMatrix(stream, n, p, width, precision, xnew);
    fprintf(stream, "Distance Matrix D(X)\n\n");
    smacofPrintSymmetricMatrix(stream, n, width, precision, dmat);
    fprintf(stream, "Fit Matrix B(X)\n\n");
    smacofPrintSymmetricMatrix(stream, n, width, precision, bmat);
    return;
}
