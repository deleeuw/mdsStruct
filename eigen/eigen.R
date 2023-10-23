dyn.load("smacofTorgerson.so")

hilbertMatrix <- function(n) {
  return(1 / (outer(1:n, 1:n, "+") - 1))
}

adhocJacobi <- function(n) {
  a <- hilbertMatrix(n)[outer(1:n, 1:n, ">=")]
  h <- .C("smacofJacobi",
          n = as.integer(n),
          a = as.double(a),
          evec = as.double(rep(0, n * n)),
          itmax = as.integer(100),
          eps = as.double(1e-6),
          verbose = as.integer(TRUE)
          )
  lbd <- a[c(1, 1 + cumsum(n:2))]
  return(list(eval = h$a[c(1, 1 + cumsum(n:2))], evec = matrix(h$evec, n, n)))
}

adhocBauer <- function(n, p = 2) {
  a <- hilbertMatrix(n)[outer(1:n, 1:n, ">=")]
  x <- rnorm(n * p)
  h <- .C("smacofSimultaneousIteration",
          a = as.double(a),
          x = as.double(x),
          n = as.integer(n),
          p = as.integer(p),
          itmax = as.integer(100),
          eps = double(1e-6),
          verbose = as.integer(FALSE)
          )
  return(h)
}