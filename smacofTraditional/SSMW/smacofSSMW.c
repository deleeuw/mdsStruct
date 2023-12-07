#include "../Weighted/smacofWeighted.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    printf("%s\n", *argv);
    char parname[80], deltaname[80], xoldname[80], weightname[80], outname[80];
    int n = 0, p = 0, init = 0, itmax = 0, feps = 0, ceps = 0, verbose = 0,
        width = 0, precision = 0, relax = 0, evalbmat = 0, comphessian = 0;
    strcat(strcpy(parname, name), "Parameters.txt");
    char *iterstring = (char *)malloc((size_t)256 * sizeof(char));
    (void)smacofReadParameterFile(parname, &n, &p, &itmax, &init, &feps, &ceps,
                                  &width, &precision, &verbose, &relax,
                                  &evalbmat, &comphessian);
    strcat(strcpy(deltaname, name), "Delta.txt");
    double *deltavec = smacofMakeAnyVector(n * n);
    (void)smacofReadInputFile(deltaname, deltavec);
    double **delta = smacofMakeAnyMatrix(n, n);
    (void)smacofFromAnyRtoC(n, n, deltavec, delta);
    strcat(strcpy(weightname, name), "Weights.txt");
    double *weightvec = smacofMakeAnyVector(n * n);
    (void)smacofReadInputFile(weightname, weightvec);
    double **weights = smacofMakeAnyMatrix(n, n);
    (void)smacofFromAnyRtoC(n, n, weightvec, weights);
    strcat(strcpy(xoldname, name), "Xold.txt");
    double *xoldvec = smacofMakeAnyVector(n * p);
    (void)smacofReadInputFile(xoldname, xoldvec);
    double **xold = smacofMakeAnyMatrix(n, p);
    (void)smacofFromAnyRtoC(n, p, xoldvec, xold);
    (void)smacofFreeAnyVector(deltavec);
    (void)smacofFreeAnyVector(xoldvec);
    (void)smacofFreeAnyVector(weightvec);
    double **dmat = smacofMakeAnyMatrix(n, n);
    double **bmat = smacofMakeAnyMatrix(n, n);
    double **vmat = smacofMakeAnyMatrix(n, n);
    double **vinv = smacofMakeAnyMatrix(n, n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    // now we are getting serious
    (void)smacofSSMWEngine(n, p, delta, weights, xold, xnew, dmat, bmat, vmat,
                           vinv, init, itmax, feps, ceps, verbose, relax,
                           iterstring);
    // phew
    strcat(strcpy(outname, name), "Output.txt");
    FILE *stream = fopen(outname, "w");
    assert(stream != NULL);
    (void)smacofWeightedWriteOutputFile(stream, n, p, width, precision, delta,
                                        weights, xnew, dmat, bmat, iterstring);
    if (evalbmat) {
        (void)smacofWeightedWriteEvalBmat(stream, n, width, precision, bmat,
                                          vmat);
    }
    if (comphessian) {
        double ****hessian = NULL;
        (void)smacofWeightedHessian(n, p, delta, weights, dmat, bmat, vmat,
                                    xnew, hessian);
        printf("I got here\n");
        (void)smacofPrintAnyMatrix(n, n, width, precision, hessian[0][0]);
    }
    (void)smacofFreeAnyMatrix(n, delta);
    (void)smacofFreeAnyMatrix(n, weights);
    (void)smacofFreeAnyMatrix(n, dmat);
    (void)smacofFreeAnyMatrix(n, bmat);
    (void)smacofFreeAnyMatrix(n, xold);
    (void)smacofFreeAnyMatrix(n, xnew);
    (void)smacofFreeAnyMatrix(n, vmat);
    (void)smacofFreeAnyMatrix(n, vinv);
    free(iterstring);
    fclose(stream);
    return EXIT_SUCCESS;
}
