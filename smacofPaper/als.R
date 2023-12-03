do <- function(x, y, n) {
  f <- rep(0, n + 1)
  for (k in 0:n) {
    f[k + 1] <- ((x ^ k) * (y ^ (n - k))) ^ (1 / n)
  }
  plot(0:n, f)
}