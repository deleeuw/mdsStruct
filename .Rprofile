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
  ssmw()
  system("make rlib")
}

clean <- function() {
  system("make clean")
}

run <- function() {
  ssmw()
  system("./smacofSSMW")
}
