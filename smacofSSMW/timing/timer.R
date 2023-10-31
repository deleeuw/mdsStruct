


timer <- function(dataR,
                  dataRC,
                  itmax = 10000,
                  eps = 15) {
  hdataR <- smacofR(dataR, itmax = itmax, eps1 = eps)
  cat("************************************************************************\n")
  cat(
    "R  itel ",
    formatC(hdataR$itel, digits = 3, format = "d"),
    formatC(hdataR$loss, digits = 15, format = "f"),
    "\n"
  )
  hdataRC <- smacofRC(dataRC, itmax = itmax, eps1 = eps)
  cat(
    "RC itel ",
    formatC(hdataRC$itel, digits = 3, format = "d"),
    formatC(hdataRC$snew, digits = 15, format = "f"),
    "\n\n"
  )
  mdata <-
    microbenchmark(smacofR(dataR, itmax = itmax, eps1 = eps),
                   smacofRC(dataRC, itmax = itmax, eps1 = eps))
  mdataR <- median(split(mdata, mdata$expr)[[1]]$time)
  mdataRC <- median(split(mdata, mdata$expr)[[2]]$time)
  cat(
    "R median time ",
    formatC(mdataR, digits = 0, format = "f"),
    "RC median time ",
    formatC(mdataRC, digits = 0, format = "f"),
    "R/RC ratio",
    formatC(
      mdataR / mdataRC,
      width = 15,
      digits = 10,
      format = "f"
    ),
    "\n"
  )
  cat("************************************************************************\n\n")
}