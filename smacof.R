delta <- rep(1.0, 6)
weights <- rep(1.0, 6)
xold <- matrix(c(1:4,4:1), 4, 2)
ii <- c(2, 3, 4, 3, 4, 4)
jj <- c(1, 1, 1, 2, 2, 3)
itmax = 100
eps1 = 10
eps2 = 6
verbose = TRUE

dyn.load("smacofCore.so")

smacof <- function(delta, weights, xold, ii, jj, itmax, eps1, eps2, verbose) {
  m <- length(delta)
  n <- nrow(xold)
  p <- ncol(xold)
  h <- .C("smacofEngine",
          delta = as.double(delta),
          weights = as.double(weights),
          xold = as.double(xold),
          ii = as.integer(ii),
          jj = as.integer(jj),
          m = as.integer(m),
          n = as.integer(n),
          p = as.integer(p),
          itmax = as.integer(itmax),
          eps1 = as.integer(eps1),
          eps2 = as.integer(eps2),
          verbose = as.integer(verbose))
  }
  
