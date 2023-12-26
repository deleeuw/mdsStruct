library(MASS)

xold <- matrix(c(
   +1.9404,   -0.4662,
   -0.5406,   -0.0697,
   -1.1625,  -1.0741,
   -0.2372,   +1.6100), 4, 2, byrow = TRUE)

xnew <- matrix(c(
  +1.7741,   -0.5397,
  -0.6979,   -0.1611,
  -0.8667,   -1.0177,
  -0.2095,   +1.7186), 4, 2, byrow = TRUE);

dold <- as.matrix(dist(xold))

delta <-matrix(c(
   0,3,2,3,
   3,0,1,2,
   2,1,0,3,
   3,2,3,0), 4, 4);

w <- 1 - diag(4);

v <- -w
diag(v) <- -rowSums(v)
vinv <- ginv(v)

b <- -w * delta / (dold + diag(4))
diag(b) <- -rowSums(b)

xcor <- vinv %*% b %*% xold