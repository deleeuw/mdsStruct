void smacofSimultaneousIteration(double *cross, double *xold, const int *pn,
                                 const int *pp, const int *itmax,
                                 const double *eps, const bool *verbose) {
    int n = *pn, p = *pp, np = n * p, itel = 1;
    double oldsum = 0.0, newsum = 0.0, maxdiff = 0.0;
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            if (i == s) {
                xold[MINDEX(i, s, n)] = 1.0;
            } else {
                xold[MINDEX(i, s, n)] = 0.0;
            }
        }
    }
    double *xnew = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *r = (double *)calloc((size_t)p, (size_t)sizeof(double));
    (void)smacofGramSchmidt(xold, r, &n, &p);
    oldsum = 0.0;
    while (true) {
        (void)smacofMultiplySDCLMatrix(cross, xold, xnew, &n, &p);
        (void)smacofGramSchmidt(xnew, r, &n, &p);
        (void)smacofMaxConfigurationDifference(xold, xnew, &n, &p, &maxdiff);
        newsum = 0.0;
        for (int i = 1; i <= p; i++) {
            newsum += r[VINDEX(i)];
        }
        if (verbose) {
            printf("itel %3d oldsum %15.10f newsum %15.10f\n", itel, oldsum,
                   newsum);
        }
        if ((itel == *itmax) || ((newsum - oldsum) < *eps)) {
            break;
        }
        itel++;
        oldsum = newsum;
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
    }
    for (int i = 1; i <= n; i++) {
        for (int s = 1; s <= p; s++) {
            xold[MINDEX(i, s, n)] = xnew[MINDEX(i, s, n)] * sqrt(r[VINDEX(s)]);
        }
    }
    free(xnew);
    free(r);
    return;
}
