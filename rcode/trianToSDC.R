trianToSDC <- function(x, n) {
  k <- 1
  s <- matrix(0, n, n)
  for (j in 1:(n - 1)) {
    for(i in (j + 1):n) {
      s[i, j] <- s[j, i] <- -x[k]
      k <- k + 1
    }
  }
  diag(s) <- -rowSums(s)
  return(s)
}