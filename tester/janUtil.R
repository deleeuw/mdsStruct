mPrint <- function(x,
                   digits = 6,
                   width = 10,
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

butLast <- function(x, m = 1) {
  return(rev(rev(x)[-(1:m)]))
}

butFirst <- function(x, m = 1) {
  return(x[-(1:m)])
}