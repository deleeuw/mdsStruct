#include "smacof.h"

struct couple {
    double value;
    int index;
};

int smacofCoComp(const void *px, const void *py) {
    double x = ((struct couple *)px)->value;
    double y = ((struct couple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void smacofCoSort(double *x, int *k, const int n) {
    struct couple *xi =
        (struct couple *)calloc((size_t)n, (size_t)sizeof(struct couple));
    for (int i = 0; i < n; i++) {
        xi[i].value = x[i];
        xi[i].index = i;
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(struct couple), smacofCoComp);
    for (int i = 0; i < n; i++) {
        x[i] = xi[i].value;
        k[i] = xi[i].index;
    }
    free(xi);
}

void smacofSqueeze(const int n, const double lowend, const double highend,
                   double *x) {
    double min = INFINITY, max = -INFINITY;
    for (int i = 0; i < n; i++) {
        max = MAX(max, x[i]);
        min = MIN(min, x[i]);
    }
    for (int i = 0; i < n; i++) {
        double a = (highend - lowend) / (max - min);
        x[i] = a * (x[i] - min) + lowend;
    }
    return;
}

void smacofShepardPlot(const int n, double **delta, double **dhat,
                       double **dmat) {
    int m = n * (n - 1) / 2;
    double *deltavec = smacofMakeVector(m + 1);
    double *dmatvec = smacofMakeVector(m + 1);
    double *dhatvec = smacofMakeVector(m + 1);
    int *indexvec = (int *)calloc((size_t)m, sizeof(int));
    (void)smacofSymmetricCtoR(n, dmat, dmatvec);
    (void)smacofSymmetricCtoR(n, dhat, dhatvec);
    (void)smacofSymmetricCtoR(n, delta, deltavec);
    dmatvec[m] = 0.0;
    dhatvec[m] = 0.0;
    deltavec[m] = 0.0;
    (void)smacofSqueeze(m, 0.0, 10.0, deltavec);
    (void)smacofSqueeze(m, 0.0, 10.0, dhatvec);
    (void)smacofSqueeze(m, 0.0, 10.0, dmatvec);
    (void)smacofCoSort(deltavec, indexvec, m + 1);

    plPlotter *plotter;
    plPlotterParams *plotter_params;
    FILE *outfile = fopen("shepardplot.ps", "w");

    plotter_params = pl_newplparams();
    pl_setplparam(
        plotter_params, "PAGESIZE",
        "letter, xsize = 8.5in, ysize = 11in, xorigin = 0in, yorigin = 0in");
    plotter = pl_newpl_r("ps", stdin, outfile, stderr, plotter_params);

    pl_openpl_r(plotter);
    pl_fspace_r(plotter, 0.0, 0.0, 10.0, 10.0);

    pl_fontname_r(plotter, "Times-Bold");
    pl_ffontsize_r(plotter, 0.2);
    pl_fline_r(plotter, 0, 0, 0, 10);
    pl_fline_r(plotter, 0, 10, 10, 10);
    pl_fline_r(plotter, 10, 10, 10, 0);
    pl_fline_r(plotter, 10, 0, 0, 0);
    for (int i = 0; i < m; i++) {
        pl_fmove_r(plotter, deltavec[i], dmatvec[indexvec[i]]);
        pl_alabel_r(plotter, 'c', 'c', "*");
    }
    for (int i = 1; i < m; i++) {
        pl_pencolorname_r(plotter, "red");
        pl_flinewidth_r(plotter, 0.05);
        pl_fline_r(plotter, deltavec[i - 1], dhatvec[indexvec[i - 1]],
                   deltavec[i], dhatvec[indexvec[i]]);
    }

    pl_closepl_r(plotter);
    pl_deletepl_r(plotter);
    (void)smacofFreeVector(dhatvec);
    (void)smacofFreeVector(dmatvec);
    (void)smacofFreeVector(deltavec);
    free(indexvec);
    return;
}


void smacofConfigurationPlot(const int n, const int dim1, const int dim2, double **xnew) {  
    double max1 = -INFINITY, max2 = -INFINITY;
    double min1 = INFINITY, min2 = INFINITY;  
    for (int i = 0; i < n; i++) {
        max1 = MAX(max1, xnew[i][dim1]);
        min1 = MIN(min1, xnew[i][dim1]);
        max2 = MAX(max2, xnew[i][dim2]);
        min2 = MIN(min2, xnew[i][dim2]);
    }
    plPlotter *plotter;
    plPlotterParams *plotter_params;
    FILE *outfile = fopen("configurationplot.ps", "w");

    plotter_params = pl_newplparams();
    pl_setplparam(
        plotter_params, "PAGESIZE",
        "letter, xsize = 8.5in, ysize = 11in, xorigin = 0in, yorigin = 0in");
    plotter = pl_newpl_r("ps", stdin, outfile, stderr, plotter_params);

    pl_openpl_r(plotter);
    pl_fspace_r(plotter, min1, min2, max1, max2);

    pl_fontname_r(plotter, "Times-Bold");
    pl_ffontsize_r(plotter, 0.01);
    pl_fline_r(plotter, min1, min2, max1, min2);
    pl_fline_r(plotter, max1, min2, max1, max2);
    pl_fline_r(plotter, max1, max2, min1, max2);
    pl_fline_r(plotter, min1, max2, min1, min2);
    char k[5];
    for (int i = 0; i < n; i++) {
        sprintf(k, "%3d", i + 1);
        pl_fmove_r(plotter, .9 * xnew[i][dim1], .9 * xnew[i][dim2]);
        pl_alabel_r(plotter, 'c', 'c', k);
    }
    pl_closepl_r(plotter);
    pl_deletepl_r(plotter);
    return;
}
