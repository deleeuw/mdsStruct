library(MASS)
jmat <- diag(4) - (1 / 4)
delta <- matrix(c(0, 1, sqrt(2), 1, 1, 0, 1, sqrt(2), sqrt(2), 1, 0, 1, 1, sqrt(2), 1, 0), 4, 4)
weights <- 1 - diag(4)
weights <- 2 * weights / sum(weights)
v <- -weights
diag(v) <- -rowSums(v)
vinv <- ginv(v)
delta <- sqrt(2) * delta / sqrt(sum(weights * delta ^ 2))
xini <- matrix(c(0.0, -2.0, -2.0, 4.0, 1.0, -1.0, -2.0, 2.0), 4, 2)
dini <- as.matrix(dist(xini))
lbd <- sum(weights * delta * dini) / sum(weights * dini ^ 2)
dold <- dini <- lbd * dini
xold <- xini <- lbd * xini
sold <- sum(weights * (delta - dold) ^ 2) / 2
bold <- -weights * delta / (dold + diag(4))
diag(bold) <- -rowSums(bold)
xnew <- vinv %*% bold %*% xold
dnew <- as.matrix(dist(xnew))
bnew <- -weights * delta / (dnew + diag(4))
diag(bnew) <- -rowSums(bnew)
snew <- sum(weights * (delta - dnew) ^ 2) / 2

