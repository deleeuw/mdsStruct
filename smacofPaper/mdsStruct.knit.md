---
title: "Notes on the C Version of Smacof"
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started October 10 2023, Version of December 01, 2023'
output:
  bookdown::pdf_document2:
    latex_engine: xelatex
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 4
    number_sections: yes
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
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

# Unweighted, full matrix

$W=E-I$ $V=(I-E)+(n-1)I=nI-E=nJ$. Thus $V^+=n^{-1}J$ and $V^+B=n^{-1}B$.

# Unweighting

\define{align*}
\sigma(d)&=\sum_{k=1}^Kw_k(\delta_k-d_k)^2\\
&=\sum_{k=1}^Kw_k(\delta_k-\overline{d}_k)^2-2\sum_{k=1}^Kw_k(\delta_k-\overline{d}_k)(d_k-\overline{d}_k)+\sum_{k=1}^Kw_k(d_k-\overline{d}_k)^2.
\end{align*}
Now suppose $w_k\leq w_\star$ and define
$$
r_k=\frac{w_k}{w_\star}(\delta_k-\overline{d}_k)
$$
\begin{align*}
\sigma(d)&\leq\sigma(\overline{d})+w_\star\left\{\sum_{k=1}^K(d_k-\overline{d}_k)^2-\sum_{k=1}^Kr_k(d_k-\overline{d}_k)\right\}\\&=\sigma(\overline{d})+w_\star\left\{\sum_{k=1}^K(d_k-(\overline{d}_k+r_k))^2-\sum_{k=1}^Kr_k^2\right\}.
\end{align*}
Note 
$$
\delta_k(X):=\frac{w_k}{w_\star}\delta_k+(1-\frac{w_k}{w_\star})d_k(X)
$$
So given $X$ computed the adjusted dissimilarities $\Delta(X)$ and perform one or more unweighted smacof
steps to compute the update $X^+$. Of course if all $w_k$ are the same then $\Delta(X)=\Delta$ and we
compute a regular unweighted smacof. It is of some interest to study how many smacof steps to make
before computing a new $\Delta(X)$.

# Thoughts on ALS

In nonmetric MDS we minimize
$$
\sigma(\Delta,X)=\sum_{k=1}^K w_k(\delta_k-d_k(X))^2
$$
over the configurations $X$ and the transformed dissimilarities $\Delta$,
where we assume $\Delta\in\mathfrak{D}$, with $\mathfrak{D}$ the intersection of
a convex cone and a sphere.

## The Single-Step approach

@kruskal_64a, @kruskal_64b defines
$$
\sigma_\star(X):=\min_{\Delta\in\mathfrak{D}}\sigma(\Delta,X)
$$
and
$$
\Delta(X)=\mathop{\text{argmin}}_{\Delta\in\mathfrak{D}}\sigma(\Delta,X)
$$
Thus
$$
\sigma_\star(X)=\sigma(\Delta(X),X)
$$
which is now a function of $X$ only. Under some conditions, which are usually true in MDS,
$$
\mathcal{D}\sigma_\star(X)=\mathcal{D}_2\sigma(\Delta(X),X)
$$
where $\mathcal{D}_2\sigma(\Delta(X),X)$ are the partials of $\sigma$ with respect to $X$.
Thus the partials of $\sigma_\star$ can be computed by evaluating the partials of $\sigma$ 
with repesct to $X$ at $(X,\Delta(X))$. This has created much confusion in the past. We can now solve the problem of minimizing $\sigma_\star$, which is a function of $X$ alone. 

I think Guttman calls this the *single step approach*. A variation of Kruskal's single-step approach defines
$$
\sigma_G(X)=\sum_{k=1}^Kw_k(\delta_k^\#(X)-d_k(X))^2
$$
where the $\delta_k^\#(X)$ are *Guttman's rank images*, i.e. the permutation of the
$d_k(X)$ that makes it monotone with the $\delta_k$ (@guttman_68). Or, alternatively, we define
$$
\sigma_S(X):=\sum_{k=1}^Kw_k(\delta_k^\%(X)-d_k(X))^2
$$
where the $\delta_k^\%(X)$ are *Shepard's rank images*, i.e. the permutation of
the $\delta_k$ that makes it monotone with the $d_k(X)$ (@shepard_62a, @shepard_62b).

The Shepard and Guttman alternatives are computationally more intricate and more complicated
than the Kruskal *monotone regression* approach, mostly because of problems with uniqueness and differentiation, but they are obviously both single step approaches.

## The Two-step Approach

The *two-step approach* or *alternating least squares* approach alternates minimization
of $\sigma(\Delta,X)$ over $X$ for our current best estimate of $\Delta$ with
minimization of $\sigma(\Delta,X)$ over $\Delta\in\mathfrak{D}$ for our current best
value of $X$. Thus an update looks like
\begin{align}
\Delta^{(k)}&=\mathop{\text{argmin}}_\Delta\sigma(\Delta,X^{(k)})(\#eq:step1),\\
X^{(k+1)}&=\mathop{\text{argmin}}_X\sigma(\Delta^{(k)},X)(\#eq:step2).
\end{align}
This approach to MDS was in the air since the early (unsuccessful) attempts around 1968 of Young and De Leeuw to combine Torgerson's classic metric MDS method with Kruskal's monotone regression transformation.

As formulated, however, there are some problems with the ALS algorithm.
Step \@ref(eq:step1) is easy to carry out, using monotone regression. Step \@ref(eq:step2) means solving a metric scaling problem,
which is an iterative proces that requires an infinite number of iterations. Thus, what is usually
implemented, is to combine step \@ref(eq:step1) with one of more iterations of a convergent iterative procedure
for metric MDS, such as smacof. If we take only one of these *inner iterations* the algorithm
becomes indistinguishable from Kruskal's single step method. This has also created much confusion in the
past.

It is somewhat worrisome that in the ALS approach we solve the first subproblem \@ref(eq:step1) exactly, while we take only a single step towards the solution for given $\Delta$ in the second subproblem \@ref(eq:step2). If we have an
infinite iterative procedure to compute the optial $\Delta\in\mathfrak{D}$ for given $X$, then
a more balanced approach is to take several inner iterations in the first step and several
innewr iterations in the second step. How many of each, nobody knows.

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

# smacofMaximumSum

Maximize
$$
\jis w_{ij}\delta_{ij}^2d_{ij}^2(X)
$$
over $X$ with $\text{tr}\ (X'X)^2=1$. This gives $BX=X\Lambda$, with
$$
B=\jis w_{ij}\delta_{ij}^2A_{ij}
$$
There is also a non-metric version. Maximize over $\text{tr}\ (X'X)^2$
$$
\jis\mathop{\sum\sum}_{1\leq k<l\leq n} w_{ij,kl}\text{sign}(\delta_{ij}-\delta_{kl})(d_{ij}^2(X)-d_{kl}^2(X))
$$
which simplifies to
$$
2\ \jis d_{ij}^2(X)\mathop{\sum\sum}_{1\leq k<l\leq n}w_{ij,kl}\text{sign}(\delta_{ij}-\delta_{kl})
$$
Simplifies more $w_{ij,kl}=w_{ij}$ or $_{ij,kl}=w_{ij}w_{kl}$

Also note
$$
\rho(X)=\jis w_{ij}\delta_{ij}d_{ij}(X)=
\jis \frac{w_{ij}}{\delta_{ij}d_{ij}(X)}\delta_{ij}^2d_{ij}^2(X)\approx\jis
\frac{w_{ij}}{\delta_{ij}^2}\delta_{ij}^2d_{ij}^2(X)=\eta^2(X)
$$

# smacofElegant

# smacofAdjustDiagonal

# smacofImpute

# smacofHildreth



@hildreth_57

Consider the QP problem of minimizing $f(x)=\frac12(x-y)'W(x-y)$
over all $x\in\mathbb{R}^n$ satisfying $Ax\geq 0$, where $A$ is $m\times n$. Wlg we can assume $a_j'Wa_j=1$. The Lagrangian is
$$
\mathcal{L}(x,\lambda)=\frac12(x-y)'W(x-y)-\lambda'Ax.
$$
$$
\max_{\lambda\geq 0}\mathcal{L}(x,\lambda)=\begin{cases}\frac12(x-y)'W(x-y)&\text{ if }Ax\geq 0,\\
+\infty&\text{ otherwise}.\end{cases}
$$
and thus
$$
\min_{x\in\mathbb{R}^n}\max_{\lambda\geq 0}\mathcal{L}(x,\lambda)=\min_{Ax\geq 0}\frac12(x-y)'W(x-y).
$$
By duality
$$
\min_{x\in\mathbb{R}^n}\max_{\lambda\geq 0}\mathcal{L}(x,\lambda)=\max_{\lambda\geq 0}\min_{x\in\mathbb{R}^n}\mathcal{L}(x,\lambda).
$$
The inner minimum over $x$ is attained for 
$$
x=y+W^{-1}A'\lambda,
$$ 
and is equal to
$$
\min_{x\in\mathbb{R}^n}\mathcal{L}(x,\lambda)=-\frac12\lambda'AW^{-1}A'\lambda-\lambda'Ay.
$$
Thus
$$
\max_{\lambda\geq 0}\min_{x\in\mathbb{R}^n}\mathcal{L}(x,\lambda)=
-\frac12\min_{\lambda\geq 0}\left\{(Wy+A'\lambda)'W^{-1}(Wy+A'\lambda)-y'Wy\right\}
$$
We minimize $h(\lambda)=(Wy+A'\lambda)'W^{-1}(Wy+A'\lambda)$ with coordinate descent. Let $\lambda_j(\eps)=\lambda+\eps e_j$. Then
$$
h(\lambda_j(\eps))=(Wy+A'\lambda+\eps a_j)'W^{-1}()=\eps^2a_j'W^{-1}a_j+2\eps a_j'W^{-1}(y+W^{-1}A'\lambda)+
$$
which must be minimized over $\eps\geq-\lambda_j$. So the minimum is attained at
$$
\eps=-\frac{a_j'W^{-1}x}{a_j'W^{-1}a_j}
$$
with $x=y+W^{-1}A'\lambda$ (cf ...), provided ... satisfies .. Otherwise $\eps=-\lambda_j$. Now update both $\lambda$ and $x$,
and go to the next $j$.

# smacofDykstra

# smacofJacobi

taken from @deleeuw_E_17o

partial jacobi

# smacofIndividualDifferenceModels.c

\begin{align}
X_k&=X,\\
X_k&=X\Lambda_k \text{ with }\Lambda_k\text{ diagonal},\\
X_k&=XC_k,\\
X_k&=X\Lambda_k Y' \text{ with }\Lambda_k\text{ diagonal},\\
X_k&=XC_kY',\\
X_k&=X\Lambda_k Y_k'\text{ with }\Lambda_k\text{ diagonal}\\
X_k&=XC_k Y'\text{ with } C_k=\sum_{s=1}^r z_{ks}H_s.
\end{align}

# smacofSort


# Code

## smacofSort.R


## smacofSort.c


```c

```

# References
