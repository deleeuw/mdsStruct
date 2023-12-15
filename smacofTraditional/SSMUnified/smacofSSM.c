#include "smacof.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1], *progname = &(argv[0])[6];
    char parname[SSIZE], deltaname[SSIZE], xoldname[SSIZE], outname[SSIZE],
        weightsname[SSIZE];
    int n = 0, p = 0, init = 0, itmax = 0, feps = 0, ceps = 0, verbose = 0,
        width = 0, precision = 0, relax = 0, transform = 0, degree = 0,
        ordinal = 0, weights = 0;
    char *iterstring = (char *)malloc((size_t)OSIZE * sizeof(char));
    strcat(strcpy(parname, name), "Parameters.txt");
    FILE *parameterfile = fopen(parname, "r");
    assert(smacofReadParameterFile != NULL);
    (void)smacofReadParameterFile(parameterfile, &n, &p, &itmax, &init, &feps,
                                  &ceps, &width, &precision, &verbose, &relax,
                                  &transform, &degree, &ordinal, &weights);
    fclose(parameterfile);
    int m = n * (n - 1) / 2;
    strcat(strcpy(deltaname, name), "Delta.txt");
    FILE *deltafile = fopen(deltaname, "r");
    assert(deltafile != NULL);
    double *deltavec = smacofMakeVector(m);
    (void)smacofReadInputFile(deltafile, deltavec);
    double **delta = smacofMakeSymmetricMatrix(n);
    (void)smacofSymmetricRtoC(n, deltavec, delta);
    fclose(deltafile);
    double **w = NULL;
    if (weights) {
        w = smacofMakeSymmetricMatrix(n);
        strcat(strcpy(weightsname, name), "Weights.txt");
        FILE *weightsfile = fopen(weightsname, "r");
        assert(weightsfile != NULL);
        double *weightsvec = smacofMakeVector(m);
        (void)smacofReadInputFile(weightsfile, weightsvec);
        (void)smacofSymmetricRtoC(n, weightsvec, w);
        (void)smacofFreeVector(weightsvec);
        fclose(weightsfile);
    } else {
        w = smacofMakeSymmetricMatrix(0);
    }
    double **xold = smacofMakeAnyMatrix(n, p);
    if (init == HAVE_INIT) {
        strcat(strcpy(xoldname, name), "Xold.txt");
        FILE *xoldfile = fopen(xoldname, "r");
        assert(xoldfile != NULL);
        double *xoldvec = smacofMakeVector(n * p);
        (void)smacofReadInputFile(xoldfile, xoldvec);
        (void)smacofFromAnyRtoC(n, p, xoldvec, xold);
        (void)smacofFreeVector(xoldvec);
        fclose(xoldfile);
    }
    double **basis = NULL;
    if (transform == POLYNOMIAL) {
        double *y = smacofMakeVector(m);
        double dmax = -INFINITY, dmin = INFINITY;
        for (int i = 0; i < m; i++) {
            dmax = MAX(dmax, deltavec[i]);
            dmin = MIN(dmin, deltavec[i]);
        }
        for (int i = 0; i < m; i++) {
            y[i] = (deltavec[i] - dmin) / (dmax - dmin);
        }
        basis = smacofMakeAnyMatrix(m, degree);
        (void)smacofBernsteinBase(m, degree, y, ordinal, basis);
        (void)smacofFreeVector(y);
    }
    double **dmat = smacofMakeAnyMatrix(n, n);
    double **dhat = smacofMakeAnyMatrix(n, n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    // now we are getting serious
    (void)smacofSSMEngine(n, p, delta, w, xold, xnew, dmat, dhat, basis, init, itmax,
                          feps, ceps, verbose, relax, weights, transform, iterstring);
    // phew
    strcat(strcat(strcpy(outname, name), progname), "Output.txt");
    FILE *stream = fopen(outname, "w");
    assert(stream != NULL);
    (void)smacofWriteOutputFile(stream, n, p, weights, width, precision, delta,
                                w, dhat, xnew, dmat, iterstring);
    (void)smacofFreeMatrix(n, delta);
    (void)smacofFreeMatrix(n, dmat);
    (void)smacofFreeMatrix(n, dhat);
    (void)smacofFreeMatrix(n, xold);
    (void)smacofFreeMatrix(n, xnew);
    if (weights) {
        (void)smacofFreeMatrix(n, w);
    }
    free(iterstring);
    fclose(stream);
    return EXIT_SUCCESS;
}
