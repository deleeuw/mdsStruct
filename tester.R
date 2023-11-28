r <-
  matrix(c(+6.00,-4.00,-2.00,-4.00,+6.00,-2.00,-2.00,-2.00,+4.00),
         3,
         3)
x <- matrix(c(1, 2, 3, 3, 2, 1), 3, 2)


gs <- function(x) {
  p <- ncol(x)
  q <- matrix(0, p, p)
  for (s in 1:p) {
    if (s > 1) {
      for (t in 1:(s - 1)) {
        ss <- sum(x[, s] * x[, t])
        x[, s] <- x[, s] - ss * x[, t]
        q[t, s] <- ss
      }
    }
    ss <- sqrt(sum(x[, s] ^ 2))
    q[s, s] <- ss
    x[, s] <- x[, s] / ss
  }
  return(list (x = x, q = q))
}
