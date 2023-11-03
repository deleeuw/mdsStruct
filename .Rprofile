make <- function() {
  system("make")
}

ssmw <- function() {
  setwd("/Users/deleeuw/Desktop/projects/mdsStruct/smacofSSMW/ccode")
}

common <- function() {
  setwd("/Users/deleeuw/Desktop/projects/mdsStruct/smacofCommon/ccode")
}

project <- function() {
  setwd("/Users/deleeuw/Desktop/projects/mdsStruct")
}

format <- function() {
  system("make format")
}

clib <- function() {
  common()
  system("make clib")
}

rlib <- function() {
  common()
  system("make rlib")
}

pristine <- function() {
  system("make pristine")
}

run <- function() {
  ssmw()
  system("./smacofSSMW")
}
