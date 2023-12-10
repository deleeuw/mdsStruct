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
      width = 0, precision = 0, relax = 0, interval = 0;
  char *iterstring = (char *)malloc((size_t)256 * sizeof(char));
  strcat(strcpy(parname, name), "Parameters.txt");
  (void)smacofReadParameterFile(parname, &n, &p, &itmax, &init, &feps, &ceps,
                                &width, &precision, &verbose, &relax, &interval);
  strcat(strcpy(deltaname, name), "Delta.txt");
  double *deltavec = smacofMakeAnyVector(n * n);
  (void)smacofReadInputFile(deltaname, deltavec);
  double **delta = smacofMakeAnyMatrix(n, n);
  (void)smacofFromAnyRtoC(n, n, deltavec, delta);
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
  (void)smacofSSIUEngine(n, p, delta, xold, xnew, dmat, dhat, bmat, init, itmax, feps,
                         ceps, verbose, relax, interval, iterstring);
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
