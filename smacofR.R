library(MASS)

# double centers a symmetric matrix

doubleCenter <- function(x) {
  rs <- apply(x, 1, mean)
  ss <- mean(x)
  return(x - outer(rs, rs, "+") + ss)
}

# mPrint() formats a matrix (or vector, or scalar) of numbers
# for printing 

mPrint <- function(x,
                   digits = 6,
                   width = 8,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

# classical MDS

torgerson <- function(delta, p = 2) {
  e <- eigen(-.5 * doubleCenter(as.matrix(delta) ^ 2))
  l <- sqrt(pmax(0, e$values[1:p]))
  if (p == 1) {
    return(as.matrix(e$vectors[, 1] * l))
  } else {
    return(e$vectors[, 1:p] %*% diag(l))
  }
}


smacofR <-
  function(delta,
           weights = 1 - diag(nrow(delta)),
           p = 2,
           xold = NULL,
           itmax = 1000,
           eps1 = 15,
           eps2 = 10,
           verbose = FALSE) {
    if (is.null(xold)) {
      xold <- torgerson(delta, p)
    }
    n <- nrow(delta)
    deps1 <- 10 ^ -eps1
    deps2 <- 10 ^ -eps2
    weights <- 2 * weights / sum(weights)
    v <- -weights
    diag(v) <- -rowSums(v)
    vinv <- ginv(v)
    delta <- sqrt(2) * delta / sqrt(sum(weights * (delta ^ 2)))
    dold <- as.matrix(dist(xold))
    lbd <- sum(weights * delta * dold) / sum(weights * (dold ^ 2))
    xold <- lbd * xold
    dold <- lbd * dold
    sold <- sum(weights * (delta - dold) ^ 2) / 4.0
    itel <- 1
    repeat {
      b <- -weights * delta / (dold + diag(n))
      diag(b) <- -rowSums(b)
      xnew <- vinv %*% b %*% xold
      dnew <- as.matrix(dist(xnew))
      snew <- sum(weights * (delta - dnew) ^ 2) / 4.0
      cchange <- max(abs(xold - xnew))
      dchange <- max(abs(dold - dnew))
      diff <- sold - snew
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, digits = 3, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "sdif ",
          formatC(diff, digits = 10, format = "f"),
          "cchange ",
          formatC(cchange, digits = 10, format = "f"),
          "dchange ",
          formatC(dchange, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) ||
          ((diff < deps1) && (cchange < deps2))) {
        break
      }
      sold <- snew
      xold <- xnew
      dold <- dnew
      itel <- itel + 1
    }
    return(list(
      conf = xnew,
      dist = dnew,
      weights = weights,
      delta = delta,
      itel = itel
    ))
  }
