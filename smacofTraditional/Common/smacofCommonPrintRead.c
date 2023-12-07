#include "smacofCommon.h"

void smacofPrintAnyMatrix(const int n, const int p, const int width,
                          const int precision, double **x) {
    for (int i = 0; i < n; i++) {
        for (int s = 0; s < p; s++) {
            printf("%+*.*f", width, precision, x[i][s]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofFPrintAnyMatrix(FILE *stream, const int n, const int p,
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

void smacofPrintSymmetricMatrix(const int n, const int width,
                                const int precision, double **x) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            printf("%+*.*f", width, precision, x[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofFPrintSymmetricMatrix(FILE *stream, const int n, const int width,
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

void smacofPrintAnyVector(const int n, const int width, const int precision,
                          double *x) {
    for (int i = 0; i < n; i++) {
        printf("%+*.*f", width, precision, x[i]);
    }
    printf("\n\n");
    return;
}

void smacofFPrintAnyVector(FILE *stream, const int n, const int width,
                           const int precision, double *x) {
    for (int i = 0; i < n; i++) {
        fprintf(stream, "%+*.*f", width, precision, x[i]);
    }
    fprintf(stream, "\n\n");
    return;
}

// read data from a file into a one-dimensional array

void smacofReadInputFile(const char *fname, double *delta) {
    FILE *stream = fopen(fname, "r");
    if (stream == NULL) {
        printf("Error: cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    int k = 0;
    fscanf(stream, "%lf", &delta[k]);
    while (!feof(stream)) {
        k++;
        fscanf(stream, "%lf", &delta[k]);
    }
    fclose(stream);
    return;
}

void smacofReadParameterFile(const char *fname, int *n, int *p, int *itmax,
                             int *init, int *feps, int *ceps, int *width,
                             int *precision, int *verbose, int *relax,
                             int *evalbmat, int *comphessian) {
    FILE *stream = fopen(fname, "return");
    assert(stream != NULL);
    char *line = (char *)malloc(80 * sizeof(char));
    char aux1[20], aux2[20];
    int par = 0;
    size_t len = 0;
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
        if (strcmp("evalbmat", aux1) == 0) {
            *evalbmat = par;
        }
        if (strcmp("comphessian", aux1) == 0) {
            *comphessian = par;
        }
    }
    fclose(stream);
    free(line);
    return;
}

void smacofUnweightedWriteOutputFile(FILE *stream, const int n, const int p,
                                     const int width, const int precision,
                                     double **delta, double **xnew,
                                     double **dmat, double **bmat,
                                     char *iterstring) {
    fprintf(stream, "Final Iteration\n\n");
    fprintf(stream, "%s\n\n", iterstring);
    fprintf(stream, "Data Delta\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, delta);
    fprintf(stream, "Configuration X\n\n");
    smacofFPrintAnyMatrix(stream, n, p, width, precision, xnew);
    fprintf(stream, "Distance Matrix D(X)\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, dmat);
    fprintf(stream, "Fit Matrix B(X)\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, bmat);
    return;
}

void smacofUnweightedWriteEvalBmat(FILE *stream, const int n, const int width,
                                   const int precision, double **bmat) {
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    double **btmp = smacofMakeAnyMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            btmp[i][j] = bmat[i][j] * (double)(n - 1);
        }
    }
    (void)smacofJacobi(n, n, btmp, evec, eval, 100, 15, false);
    fprintf(stream, "Eigenvalues V^+B(X)\n\n");
    (void)smacofFPrintAnyVector(stream, n, width, precision, eval);
    (void)smacofFreeAnyVector(eval);
    (void)smacofFreeAnyMatrix(n, evec);
    (void)smacofFreeAnyMatrix(n, btmp);
}

void smacofWeightedWriteEvalBmat(FILE *stream, const int n, const int width,
                                 const int precision, double **bmat,
                                 double **vmat) {
    double **evec = smacofMakeAnyMatrix(n, n);
    double *eval = smacofMakeAnyVector(n);
    double **bvec = smacofMakeAnyMatrix(n, n);
    (void)smacofDoubleJacobi(n, bmat, vmat, evec, eval, 100, 15, false);
    fprintf(stream, "Eigenvalues V^+B(X)\n\n");
    (void)smacofFPrintAnyVector(stream, n, width, precision, eval);
    (void)smacofFreeAnyVector(eval);
    (void)smacofFreeAnyMatrix(n, evec);
    (void)smacofFreeAnyMatrix(n, bvec);
}

void smacofWeightedWriteOutputFile(FILE *stream, const int n, const int p,
                                   const int width, const int precision,
                                   double **delta, double **weights,
                                   double **xnew, double **dmat, double **bmat,
                                   char *iterstring) {
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    fprintf(stream, "%s\n\n", asctime(tm));
    fprintf(stream, "Final Iteration\n\n");
    fprintf(stream, "%s\n\n", iterstring);
    fprintf(stream, "Data Delta\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, delta);
    fprintf(stream, "Weights W\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, weights);
    fprintf(stream, "Configuration X\n\n");
    smacofFPrintAnyMatrix(stream, n, p, width, precision, xnew);
    fprintf(stream, "Distance Matrix D(X)\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, dmat);
    fprintf(stream, "Fit Matrix B(X)\n\n");
    smacofFPrintSymmetricMatrix(stream, n, width, precision, bmat);
    return;
}
