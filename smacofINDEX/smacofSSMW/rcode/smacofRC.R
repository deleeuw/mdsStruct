dyn.load("/Users/deleeuw/Desktop/projects/mdsStruct/smacofBuild/libsmacof.so")

smacofRCW <-
  function(delta,
           p = 2,
           n = as.integer((1 + sqrt(1 + 8 * length(delta)) / 2)),
           xold = rbind(diag(p),matrix(0, n - p, p)),
           weights = rep(1, length(delta)),
           ii = fullIndex(n)$ii,
           jj = fullIndex(n)$ii,
           itmax = 10000,
           init = 1,
           eps1 = 15,
           eps2 = 10,
           verbose = FALSE) {
    h <- .C(
      "smacofSSMWEngine",
      delta = as.double(delta),
      weights = as.double(weights),
      xold = as.double(xold),
      xnew = as.double(array(0.0, dim(xold))),
      dnew = as.double(rep(0.0, length(delta))),
      snew = as.double(0.0),
      ii = as.integer(ii),
      jj = as.integer(jj),
      m = as.integer(length(delta)),
      n = as.integer(n),
      p = as.integer(p),
      itel = as.integer(1),
      itmax = as.integer(itmax),
      init = as.integer(init),
      eps1 = as.integer(eps1),
      eps2 = as.integer(eps2),
      verbose = as.integer(verbose)
    )
    return(h)
  }

