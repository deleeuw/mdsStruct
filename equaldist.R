n <- 15
eqdistR <- 1 - diag(n)
eqdistRC <- as.vector(as.dist(eqdistR))
eqdistold <- matrix(rnorm(30), n, 2)
eqdistold <- apply(eqdistold, 2, function(x) x - mean(x))