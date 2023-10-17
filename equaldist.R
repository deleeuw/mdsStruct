n <- 15
eqdistR <- 1 - diag(n)
eqdistRC <- as.vector(as.dist(deltaR))
eqdistold <- matrix(c(1:n, n:1), n, 2)
eqdistold <- apply(eqdistold, 2, function(x) x - mean(x))