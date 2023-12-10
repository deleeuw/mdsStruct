#include "../Unweighted/smacofUnweighted.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    char *progname = &(argv[0])[6];
    char parname[80], deltaname[80], xoldname[80], outname[80];
    int n = 0, p = 0, init = 0, itmax = 0, feps = 0, ceps = 0, verbose = 0,
        width = 0, precision = 0, relax = 0, interval = 0, degree = 0,
        ordinal = 0;
    char *iterstring = (char *)malloc((size_t)256 * sizeof(char));
    strcat(strcpy(parname, name), "Parameters.txt");
    (void)smacofReadParameterFile(parname, &n, &p, &itmax, &init, &feps, &ceps,
                                  &width, &precision, &verbose, &relax,
                                  &interval, &degree, &ordinal);
    int nn = n * n, np = n * p, m = n * (n - 1) / 2;
    strcat(strcpy(deltaname, name), "Delta.txt");
    double *deltavec = smacofMakeAnyVector(nn);
    (void)smacofReadInputFile(deltaname, deltavec);
    double **delta = smacofMakeAnyMatrix(n, n);
    (void)smacofFromAnyRtoC(n, n, deltavec, delta);
    strcat(strcpy(xoldname, name), "Xold.txt");
    double *xoldvec = smacofMakeAnyVector(np);
    (void)smacofReadInputFile(xoldname, xoldvec);
    double **xold = smacofMakeAnyMatrix(n, p);
    (void)smacofFromAnyRtoC(n, p, xoldvec, xold);
    printf("%3d %3d %3d %3d\n", n, degree, np, nn);
    double **basis = smacofMakeAnyMatrix(nn, degree);
    (void)smacofBernsteinBase(nn, degree, deltavec, ordinal, basis);
    (void)smacofPrintAnyMatrix(nn, degree, 15, 10, basis);
    (void)smacofFreeAnyVector(deltavec);
    (void)smacofFreeAnyVector(xoldvec);
    exit(EXIT_SUCCESS);
    double **dmat = smacofMakeAnyMatrix(n, n);
    double **dhat = smacofMakeAnyMatrix(n, n);
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    double *coef = smacofMakeAnyVector(degree);
    exit(EXIT_SUCCESS);
    // now we are getting serious
    (void)smacofSSPUEngine(n, p, delta, basis, coef, xold, xnew, dmat, dhat, bmat, init,
                           itmax, feps, ceps, verbose, relax, degree,
                           ordinal, iterstring);
    // phew
    strcat(strcat(strcpy(outname, name), progname), "Output.txt");
    FILE *stream = fopen(outname, "w");
    assert(stream != NULL);
    (void)smacofUnweightedWriteOutputFile(stream, n, p, width, precision, delta,
                                          dhat, xnew, dmat, bmat, iterstring);
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
