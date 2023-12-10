#include "smacof.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1], *progname = &(argv[0])[6];
    char parname[SSIZE], deltaname[SSIZE], xoldname[SSIZE], outname[SSIZE];
    int n = 0, p = 0, init = 0, itmax = 0, feps = 0, ceps = 0, verbose = 0,
        width = 0, precision = 0, relax = 0, interval = 0, degree = 0,
        ordinal = 0, hasweights = 0;
    char *iterstring = (char *)malloc((size_t)OSIZE * sizeof(char));
    strcat(strcpy(parname, name), "Parameters.txt");
    (void)smacofReadParameterFile(parname, &n, &p, &itmax, &init, &feps, &ceps,
                                  &width, &precision, &verbose, &relax,
                                  &interval, &degree, &ordinal, &hasweights);
    int nn = n * n, np = n * p, m = n * (n - 1) / 2;
    strcat(strcpy(deltaname, name), "Delta.txt");
    double *deltavec = smacofMakeAnyVector(m);
    (void)smacofReadInputFile(deltaname, deltavec);
    double **delta = smacofMakeSymmetricMatrix(n);
    (void)smacofFromSymmetricRtoC(n, deltavec, delta);
    if (hasweights) {
        strcat(strcpy(deltaname, name), "Weights.txt");
        double *weightsvec = smacofMakeAnyVector(m);
        (void)smacofReadInputFile(weightsname, weightsvec);
        double **weights = smacofMakeSymmetricMatrix(n);
        (void)smacofFromSymmetricRtoC(n, weightsvec, weights);
    } else {
        double **weights = smacofMakeSymmetricMatrix(0);
    }
    strcat(strcpy(xoldname, name), "Xold.txt");
    double *xoldvec = smacofMakeAnyVector(n * p);
    (void)smacofReadInputFile(xoldname, xoldvec);
    double **xold = smacofMakeAnyMatrix(n, p);
    (void)smacofFromAnyRtoC(n, p, xoldvec, xold);
    (void)smacofFreeAnyVector(deltavec);
    (void)smacofFreeAnyVector(xoldvec);
    double **dmat = smacofMakeAnyMatrix(n, n);
    double **dhat = smacofMakeAnyMatrix(n, n);
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    // now we are getting serious
    (void)smacofSSMEngine(n, p, delta, weights, xold, xnew, dmat, dhat, bmat,
                          init, itmax, feps, ceps, verbose, relax, hasweights,
                          iterstring);
    // phew
    strcat(strcat(strcpy(outname, name), progname), "Output.txt");
    FILE *stream = fopen(outname, "w");
    assert(stream != NULL);
    (void)smacofWriteOutputFile(stream, n, p, width, precision, delta, dhat,
                                xnew, dmat, bmat, iterstring);
    (void)smacofFreeAnyMatrix(n, delta);
    (void)smacofFreeAnyMatrix(n, dmat);
    (void)smacofFreeAnyMatrix(n, dhat);
    (void)smacofFreeAnyMatrix(n, bmat);
    (void)smacofFreeAnyMatrix(n, xold);
    (void)smacofFreeAnyMatrix(n, xnew);
    free(iterstring);
    fclose(stream);
    return EXIT_SUCCESS;
}
