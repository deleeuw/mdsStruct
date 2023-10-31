---
title: "Notes on the C Version of Smacof"
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started October 10 2023, Version of October 30, 2023'
output:
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
    toc: true
    toc_depth: 4
    number_sections: yes
  bookdown::pdf_document2:
    latex_engine: lualatex
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 4
    number_sections: yes
graphics: yes
mainfont: Times New Roman
fontsize: 12pt
bibliography: ["mypubs.bib","total.bib"]
abstract: TBD
---





**Note:** This is a working paper which will be expanded/updated frequently. All suggestions for improvement are welcome. All Rmd, tex, html, pdf, R, and C files are in the public domain. Attribution
will be appreciated, but is not required. The files can be found at
https://github.com/deleeuw/mdsStruct. 

# Introduction

The loss function in (metric, least squares, Euclidean, symmetric) Multidimensional Scaling (MDS)
is
$$
\sigma(X):=\frac12\jis w_{ij}(\delta_{ij}-d_{ij}(X))^2.
$$

This assumes symmetry and it uses all elements below the diagonal of both $W$, $\Delta$, and $D(X)$.
For missing data we set $w_{ij}=0$. 


# Example

Here is a small input example.



# smacofStructure

We use this example as input for the R version of *smacofSort()*, which results in the
*smacofStructure*


The *dist* element in the data frame is zero, because we have not computed distances yet.
Given a configuration $X$ the *row* and *col* elements in the data frame allow from
straightforward computation of distances. In the case of nonmetric MDS, or more generally
in MDS with transformed dissimilarities, the *dhat* column will be filled as well.
The fact that preprocessing with *smacofSort()* gives the ordered dissimilarities, as
well as the tie blocks, is especially useful in the ordinal case, both with monotone
polynomials and monotone splines. In the metric (ratio) case there is no need for the
*dhat* column.

It is obvious how this *smacofStructure* can be adapted if there are multiway data.If we have row-conditional or matrix-conditional data, then we use one of these smacofStructures for each
row or each matrix.

# Appendix: Code

We use *qsort* to sort the rows of the input data frame by increasing delta. The sorting
is done in C, the R version is a wrapper around the compiled C code. The C code also
contains a main which analyzes the same small example as we have used in the text.
Compile with "clang -o runner -O2 smacofSort.c" and then start "runner" in the shell.
The C code uses the .C() interface in R, which can be improved using .Call(), 
but probably in this case with little gain.

# smacofHildreth

Consider the QP problem of minimizing $f(x)=\frac12(x-y)'W(x-y)$
over all $x\in\mathbb{R}^n$ satisfying $Ax\geq 0$, where $A$ is $m\times n$. The Lagrangian is
$$
\mathcal{L}(x,\lambda)=\frac12(x-y)'W(x-y)-\lambda'Ax
$$
amd we want
$$
\min_x\max_{\lambda\geq 0}\mathcal{L}(x,\lambda)=\max_{\lambda\geq 0}\min_x\mathcal{L}(x,\lambda)
$$
The inner minimum over $x$ is attained for 
$$
x=y+W^{-1}A'\lambda,
$$ 
and is equal to
$$
\min_x\mathcal{L}(x,\lambda)=-\frac12\lambda'AW^{-1}A'\lambda-\lambda'Ay
$$
Thus
$$
\max_{\lambda\geq 0}\min_x\mathcal{L}(x,\lambda)=
-\frac12\min_{\lambda\geq 0}\left\{(Wy+A'\lambda)'W^{-1}(Wy+A'\lambda)-y'Wy\right\}
$$
Same treatment for the more general problem of minimizing $g(x)=\frac12(Bx-y)'W(Bx-y)$ over
$Ax\geq 0$. Now $x=(B'WB)^{-1}(B'Wy+A'\lambda)$ and
$$
\min_x\mathcal{L}(x,\lambda)=\frac12((B(B'WB)^{-1}B'W-I)y+B(B'WB)^{-1}A'\lambda)'W(B(B'WB)^{-1}(B'Wy+A'\lambda)-y)-\lambda'A(B'WB)^{-1}(B'Wy+A'\lambda)=
$$
$$
=\frac12y'WB(B'WB)^{-1}B'Wy-2
$$
# smacofDykstra

# smacofJacobi

# smacofSort


# Code

## smacofSort.R


## smacofSort.c


```c
#include "smacof.h"

int smacofComparison(const void *px, const void *py) {
    double x = ((struct fiveTuple *)px)->delta;
    double y = ((struct fiveTuple *)py)->delta;
    return (int)copysign(1.0, x - y);
}

void smacofSort(double *delta, double *weight, int *row, int *col, int *index,
                const int *ndata) {
    int n = *ndata;
    struct fiveTuple *xi =
        (struct fiveTuple *)calloc((size_t)n, (size_t)sizeof(struct fiveTuple));
    for (int i = 0; i < n; i++) {
        xi[i].index = i;
        xi[i].row = row[i];
        xi[i].col = col[i];
        xi[i].delta = delta[i];
        xi[i].weight = weight[i];
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(struct fiveTuple),
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
```

# References
