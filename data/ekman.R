data(ekman, package = "smacof")
n <- attr(ekman, "Size")
ekmanRC <- as.vector((1 - ekman) ^ 3)
ekmanR <- as.matrix((1 - ekman) ^ 3)
fullIndex <- function(n) {
  ii <- c()
  jj <- c()
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      ii <- c(ii, i)
      jj <- c(jj, j)
    }
  }
  return(list(ii = ii, jj = jj))
}
h <- fullIndex(14)
g <- cbind(h$ii,h$jj,ekmanRC,1)

