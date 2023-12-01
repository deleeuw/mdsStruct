#include "../Unweighted/smacofUnweighted.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("I need a job name on the command line as argument\n");
        exit(EXIT_FAILURE);
    }
    char *name = argv[1];
    char parname[80], deltaname[80], xininame[80];
    int n = 0, p = 0, init = 0, itmax = 0, feps = 0, ceps = 0, verbose = 0,
        width = 0, precision = 0, relax = 0, itel = 1;
    double snew = 0.0;
    strcat(strcpy(parname, name), "Parameters.txt");
    (void)smacofReadParameterFile(parname, &n, &p, &itmax, &init, &feps, &ceps,
                                  &width, &precision, &verbose, &relax);
    strcat(strcpy(deltaname, name), "Delta.txt");
    double *deltavec = smacofMakeAnyVector(n * n);
    (void)smacofReadInputFile(deltaname, deltavec);
    double **delta = smacofMakeAnyMatrix(n, n);
    (void)smacofFromAnyRtoC(n, n, deltavec, delta);
    strcat(strcpy(xininame, name), "Xini.txt"); 
    double *xinivec = smacofMakeAnyVector(n * p);
    (void)smacofReadInputFile(xininame, xinivec);
    double **xini = smacofMakeAnyMatrix(n, p);
    (void)smacofFromAnyRtoC(n, p, xinivec, xini);
    (void)smacofFreeAnyVector(deltavec);
    (void)smacofFreeAnyVector(xinivec);
    double **dnew = smacofMakeAnyMatrix(n, n);
    double **bnew = smacofMakeAnyMatrix(n, n);
    double **xnew = smacofMakeAnyMatrix(n, p);
    (void)smacofSSMUEngine(n, p, delta, xini, xnew, dnew, bnew, init, itmax,
                         feps, ceps, verbose, relax, &itel, &snew);
    (void)smacofFreeAnyMatrix(n, delta);
    (void)smacofFreeAnyMatrix(n, dnew);
    (void)smacofFreeAnyMatrix(n, bnew);
    (void)smacofFreeAnyMatrix(n, xini);
    (void)smacofFreeAnyMatrix(n, xnew);
    return EXIT_SUCCESS;
}
