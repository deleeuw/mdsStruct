

int smacofComparison(const void *px, const void *py) {
    double x = ((struct smacofData *)px)->delta;
    double y = ((struct smacofData *)py)->delta;
    return (int)copysign(1.0, x - y);
}

void smacofSort(double *delta, double *weight, int *row, int *col, int *index,
                const int *ndata) {
    int n = *ndata;
    struct smacofData *xi = (struct smacofData *)calloc(
        (size_t)nn, (size_t)sizeof(struct smacofData));
    for (int i = 0; i < n; i++) {
        xi[i].index = i;
        xi[i].row = row[i];
        xi[i].col = col[i];
        xi[i].delta = delta[i];
        xi[i].weight = weight[i];
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(struct smacofData),
                smacofComparison);
    for (int i = 0; i < n; i++) {
        index[i] = xi[i].index;
        row[i] = xi[i].row;
        col[i] = xi[i].col;
        delta[i] = xi[i].delta;
        weight[i] = xi[i].weight;
    }
    free(xi);
    return;
}

void smacofTieBlocks(const double *delta, int *block, double *eps,
                     const int *ndata) {
    int n = *ndata;
    block[0] = 1;
    for (int i = 1; i < n; i++) {
        if (fabs(delta[i] - delta[i - 1]) < *eps) {
            block[i] = block[i - 1];
        } else {
            block[i] = block[i - 1] + 1;
        }
    }
    return;
}
