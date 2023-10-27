#include "smacof.h"

// to be called from R

void smacofEngine(double *delta, double *weights, double *xini, double *xnew,
                  double *dini, double *dnew, double *bnew, double *psnew,
                  const int *pinit, const int *pn, const int *pp, int *pitel,
                  const int *pitmax, const int *peps1, const int *peps2,
                  const bool *pverbose, const bool *prelax,
                  const bool *padjust) {
    int n = *pn, p = *pp, np = p * n, m = n * (n - 1) / 2, itel = *pitel,
        itmax = *pitmax, itmax_j = 100, init = *pinit, eps_j = 15;
    int width = 15, precision = 10;
    bool verbose = *pverbose, verbose_j = false, verbose_e = false,
         relax = *prelax, adjust = *padjust;
    double sold = 0.0, snew = *psnew, cchange = 0.0, dchange = 0.0,
           pchange = 1.0, echange = 1.0, rate = 1.0, eopt = 1.0;
    double rho = 0.0, etaold = 0.0, etanew = 0.0, chch = 0.0;
    double eps1 = pow(10.0, -(double)*peps1), eps2 = pow(10.0, -(double)*peps2),
           eps_e = pow(10.0, -15.0);
    double *xold = (double *)calloc((size_t)np, (size_t)sizeof(double));
    double *dold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *bold = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vmat = (double *)calloc((size_t)m, (size_t)sizeof(double));
    double *vinv = (double *)calloc((size_t)m, (size_t)sizeof(double));
    (void)smacofNormWeights(weights, &m);
    (void)smacofNormDelta(delta, weights, &m);
    (void)smacofMakeVMatrix(weights, vmat, pn);
    (void)smacofMPInverseSDCLMatrix(weights, vinv, pn);
    (void)smacofInitial(delta, weights, xini, pinit, pn, pp, padjust);
    (void)smacofDistance(xini, dini, pn, pp);
    (void)memcpy(xold, xini, (size_t)np * sizeof(double));
    (void)memcpy(dold, dini, (size_t)m * sizeof(double));
    (void)smacofStress(delta, weights, dold, &m, &sold);
    (void)smacofMakeBMatrix(delta, weights, dold, bold, &m);
    while (true) {
        (void)smacofGuttman(vinv, bold, xold, xnew, pn, pp);
        (void)smacofRMSDifference(xold, xnew, pn, pp, &echange);
        (void)smacofRelax(xold, xnew, &echange, &pchange, &np, &itel, &relax,
                          &rate);
        (void)smacofDistance(xnew, dnew, pn, pp);
        (void)smacofMakeBMatrix(delta, weights, dnew, bnew, &m);
        (void)smacofStress(delta, weights, dnew, &m, &snew);
        (void)smacofMaxConfigurationDifference(xold, xnew, pn, pp, &cchange);
        (void)smacofMaxDistanceDifference(dold, dnew, &m, &dchange);
        if (verbose) {
            printf(
                "itel %3d sold %12.10f sdif %+12.10f cvdf %+12.10f cchg "
                "%12.10f dchg %12.10f rate %12.10f\n",
                itel, sold, sold - snew, echange, cchange, dchange, rate);
        }
        if ((itel == itmax) || (((sold - snew) < eps1) && (cchange < eps2))) {
            break;
        }
        itel++;
        sold = snew;
        pchange = echange;
        (void)memcpy(xold, xnew, (size_t)np * sizeof(double));
        (void)memcpy(dold, dnew, (size_t)m * sizeof(double));
        (void)memcpy(bold, bnew, (size_t)m * sizeof(double));
    }
    *psnew = snew;
    *pitel = itel;
    free(vinv);
    free(vmat);
    free(dold);
    free(xold);
    free(bold);
    return;
}

int main() {
    double delta[91] = {
        0.002744, 0.195112, 0.195112, 0.551368, 0.830584, 0.804357, 0.884736,
        0.941192, 0.804357, 0.753571, 0.681472, 0.658503, 0.592704, 0.125,
        0.175616, 0.474552, 0.753571, 0.804357, 0.804357, 0.941192, 0.884736,
        0.804357, 0.704969, 0.658503, 0.636056, 0.006859, 0.148877, 0.571787,
        0.729,    0.778688, 0.941192, 0.970299, 0.941192, 0.970299, 0.857375,
        0.912673, 0.097336, 0.421875, 0.729,    0.753571, 0.941192, 0.970299,
        1,        0.970299, 0.941192, 0.884736, 0.059319, 0.328509, 0.405224,
        0.804357, 0.941192, 0.941192, 0.970299, 0.941192, 1,        0.054872,
        0.166375, 0.636056, 0.778688, 0.941192, 0.941192, 0.941192, 0.970299,
        0.019683, 0.474552, 0.636056, 0.857375, 0.941192, 0.941192, 1,
        0.300763, 0.531441, 0.884736, 0.912673, 0.941192, 0.941192, 0.074088,
        0.250047, 0.389017, 0.512,    0.456533, 0.017576, 0.125,    0.205379,
        0.373248, 0.013824, 0.054872, 0.091125, 0.003375, 0.032768, 0.013824};
    double weights[91] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double bnew[91] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double dnew[91] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double xini[28] = {0.0714285714285714,
                       0.142857142857143,
                       0.214285714285714,
                       0.285714285714286,
                       0.357142857142857,
                       0.428571428571429,
                       0.5,
                       0.571428571428571,
                       0.642857142857143,
                       0.714285714285714,
                       0.785714285714286,
                       0.857142857142857,
                       0.928571428571429,
                       1,
                       0.00510204081632653,
                       0.0204081632653061,
                       0.0459183673469388,
                       0.0816326530612245,
                       0.127551020408163,
                       0.183673469387755,
                       0.25,
                       0.326530612244898,
                       0.413265306122449,
                       0.510204081632653,
                       0.61734693877551,
                       0.73469387755102,
                       0.862244897959184,
                       1};
    double dini[91] = {
        0.0730501074656957, 0.148573671209806,  0.227541868581406,
        0.310847881871996,  0.399297853124962,  0.493607617242768,
        0.59440417775333,   0.702230634044296,  0.817553287412672,
        0.940769842580907,  1.07221788224475,   1.21218305346265,
        1.36090662690246,   0.0758472895271352, 0.155423940935998,
        0.239578711874977,  0.329071744828512,  0.42457441268095,
        0.526672975533222,  0.635874779098745,  0.752615874064726,
        0.877269176382067,  1.01015254455221,   1.15153637660977,
        1.30165050768179,   0.0798595706249925, 0.164535872414256,
        0.25474464760857,   0.351115317022148,  0.454196270784818,
        0.564461905548544,  0.682320470519386,  0.808122035641769,
        0.942166126317085,  1.0847087564015,    1.23596872959693,
        0.08491488253619,   0.175557658511074,  0.272517762470891,
        0.376307937032363,  0.487371764656704,  0.606091526731326,
        0.7327958760244,    0.867767005121196,  1.01124714239749,
        1.16344431132565,   0.0908392541569635, 0.188153968516181,
        0.292423058794022,  0.404061017820884,  0.523425625731714,
        0.650825253840897,  0.786525555198044,  0.930755449060521,
        1.08371238255973,   0.0974743529313408, 0.202030508910442,
        0.314055375439028,  0.433883502560598,  0.561803967998603,
        0.698066586795391,  0.842887408657568,  0.996453519651731,
        0.104685125146343,  0.216941751280299,  0.337082380799162,
        0.46537772453026,   0.602062434755405,  0.747340139738798,
        0.901387818865997,  0.112360793599721,  0.23268886226513,
        0.361237460853243,  0.498226759825865,  0.643848442047141,
        0.798269682269712,  0.120412486951081,  0.249113379912933,
        0.386309065228285,  0.532179788179808,  0.686883266216136,
        0.128769688409428,  0.266089894089904,  0.412129959729681,
        0.567038530181625,  0.137376653243227,  0.283519265090812,
        0.43856781985391,   0.146189273284637,  0.3013229195966,
        0.155172513528016};
    double xnew[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double snew = 0.0;
    int n = 14, p = 2, itel = 1, itmax = 100, init = 4, peps1 = 15, peps2 = 10;
    bool verbose = true, relax = true, adjust = true;
    (void)smacofEngine(delta, weights, xini, xnew, dini, dnew, bnew, &snew,
                       &init, &n, &p, &itel, &itmax, &peps1, &peps2, &verbose,
                       &relax, &adjust);
}
