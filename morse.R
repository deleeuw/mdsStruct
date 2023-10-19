data(morse, package = "smacof")
morseR <- as.matrix(1 - morse)
morseRC <- as.vector(morse)
n <- 36
xold <- matrix(rnorm(72), 36, 2)
