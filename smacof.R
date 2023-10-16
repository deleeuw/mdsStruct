data(ekman, package = "smacof")
n <- attr(ekman, "Size")
delta <- as.vector(ekman)
weights <- rep(1.0, length(delta))
xold <- matrix(c(1:n, n:1), n, 2)
k <- 1
for (j in 1:(n - 1)) {
  for (i in (j  + 1):n) {
    ii[k] <- i
    jj[k] <- j
    k <- k + 1
  }
}
itmax = 1000
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
  
smacof(delta, weights, xold, ii, jj, itmax, eps1, eps2, verbose)

dyn.unload("smacofCore.so")