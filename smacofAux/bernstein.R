bernstein <- function(d, y) {
  y <- (y - min(y)) / (max(y) - min(y))
  return(outer(y, 0:5, function(x, y) dbinom(y, 5, x)))
}