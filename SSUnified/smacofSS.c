#include "smacof.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    char parname[SSIZE], deltaname[SSIZE], xoldname[SSIZE], outname[SSIZE],
        weightsname[SSIZE], iknotname[SSIZE], timename[SSIZE];
    int n = 0, p = 0, haveinit = 0, typeinit = 0, itmax = 0, feps = 0, ceps = 0,
        verbose = 0, width = 0, precision = 0, relax = 0, degree = 0,
        ordinal = 0, weights = 0, haveknots = 0, anchor = 0, ninner = 0,
        percentiles = 0, ditmax = 0, deps = 0, dverbose = 0, writefile = 0,
        makeplots = 0, transform = 0;
    double lowend = 0.0, highend = 0.0;
    time_t t = time(NULL);
    char *iterstring = (char *)malloc((size_t)OSIZE * sizeof(char));
    strcat(strcpy(parname, name), "Parameters.txt");
    FILE *parameterfile = fopen(parname, "r");
    assert(smacofReadParameterFile != NULL);
    (void)smacofReadParameterFile(
        parameterfile, &n, &p, &itmax, &haveinit, &typeinit, &feps, &ceps,
        &width, &precision, &verbose, &relax, &ditmax, &deps, &dverbose,
        &degree, &ordinal, &weights, &haveknots, &lowend, &highend, &anchor,
        &ninner, &percentiles, &writefile, &makeplots, &transform);
    fclose(parameterfile);
    int m = n * (n - 1) / 2, order = degree + 1, nknots = ninner + (2 * order),
        ncol = nknots - order;
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
    double **vmat = NULL, **vinv = NULL;
    if (weights) {
        vmat = smacofMakeSymmetricMatrix(n);
        vinv = smacofMakeSymmetricMatrix(n);
        (void)smacofMakeVMatrix(n, w, vmat);
        (void)smacofInverseVMatrix(n, vmat, vinv);
    } else {
        vmat = smacofMakeSymmetricMatrix(0);
        vinv = smacofMakeSymmetricMatrix(0);
    }
    double **xold = smacofMakeAnyMatrix(n, p);
    if (haveinit) {
        strcat(strcpy(xoldname, name), "Xold.txt");
        FILE *xoldfile = fopen(xoldname, "r");
        assert(xoldfile != NULL);
        double *xoldvec = smacofMakeVector(n * p);
        (void)smacofReadInputFile(xoldfile, xoldvec);
        (void)smacofAnyRtoC(n, p, xoldvec, xold);
        (void)smacofFreeVector(xoldvec);
        fclose(xoldfile);
    }
    double *innerknots = smacofMakeVector(ninner);
    int *multiplicities = (int *)calloc((size_t)ninner, (size_t)sizeof(int));
    for (int i = 0; i < ninner; i++) {
        multiplicities[i] = 1;
    }
    if (haveknots == HAVE_INNER_KNOTS) {
        strcat(strcpy(iknotname, name), "Knots.txt");
        FILE *iknotfile = fopen(iknotname, "r");
        assert(iknotfile != NULL);
        (void)smacofReadInputFile(iknotfile, innerknots);
        fclose(iknotfile);
    } else {
        double *dsort = smacofMakeVector(m);
        (void)smacofMakeInnerKnots(ninner, percentiles, n, lowend, highend,
                                   dsort, innerknots);
        (void)smacofFreeVector(dsort);
    }
    if (smacofCheckIncreasing(innerknots, lowend, highend, ninner)) {
        printf("The knots file is not strictly ordered\n");
        exit(EXIT_FAILURE);
    }
    double *knots = smacofMakeVector(nknots);
    (void)smacofExtendPartition(innerknots, multiplicities, &order, &ninner,
                                &lowend, &highend, knots);
    double *basisvector = smacofMakeVector(m * ncol);
    (void)smacofBsplineBasis(deltavec, knots, &order, &nknots, &m, basisvector);
    double **basis = smacofMakeAnyMatrix(m, ncol);
    (void)smacofAnyRtoC(m, ncol, basisvector, basis);
    (void)smacofFreeVector(basisvector);
    (void)smacofFreeVector(deltavec);
    (void)smacofPrintAnyMatrix(stdout, m, ncol, 8, 4, basis);
    if (ordinal) {
        (void)smacofCumsumMatrix(m, ncol, basis);
    }
    (void)smacofPrintAnyMatrix(stdout, m, ncol, 8, 4, basis);
    double **dmat = smacofMakeSymmetricMatrix(n);
    double **dhat = smacofMakeSymmetricMatrix(n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    // now we are getting serious
    (void)smacofSSEngine(n, p, delta, w, vmat, vinv, xold, xnew, dmat, dhat,
                         basis, haveinit, typeinit, itmax, feps, ceps, verbose,
                         relax, ditmax, deps, dverbose, weights, ncol, ordinal,
                         transform, iterstring);
    // phew
    if (writefile) {
        struct tm *tm = localtime(&t);
        strftime(timename, sizeof(timename), "_%F_%H:%M", tm);
        strcat(strcat(strcpy(outname, name), timename), ".txt");
        FILE *output = fopen(outname, "w");
        assert(output != NULL);
        //(void)smacofWriteParameterFile();
        (void)smacofWriteOutputFile(output, n, p, weights, width, precision,
                                    delta, w, dhat, xnew, dmat, iterstring);
        fclose(output);
    }
    if (makeplots) {
        (void)smacofShepardPlot(n, delta, dhat, dmat);
        (void)smacofConfigurationPlot(n, 0, 1, xnew);
    }
    (void)smacofFreeMatrix(n, delta);
    (void)smacofFreeMatrix(n, dmat);
    (void)smacofFreeMatrix(n, dhat);
    (void)smacofFreeMatrix(n, xold);
    (void)smacofFreeMatrix(n, xnew);
    if (weights) {
        (void)smacofFreeMatrix(n, w);
    }
    free(iterstring);
    return EXIT_SUCCESS;
}
