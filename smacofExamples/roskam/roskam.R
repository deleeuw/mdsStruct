roskam <-
structure(list(SOC = c("1", "1", "1", "1", "7", "6", "2", "4", 
"4", "3", "4", "3", "2", "2", "7", "5", "5", "9", "9", "8", "7", 
"8", "8", "8", "7", "4", "5", "1", "2", "2", "5", "4", "5", "6", 
"8", "2", "5", "8", "5"), EDU = c("5", "3", "6", "5", "1", "1", 
"1", "1", "1", "1", "1", "2", "9", "7", "2", "7", "9", "6", "6", 
"3", "2", "7", "6", "7", "3", "7", "6", "8", "5", "5", "3", "5", 
"7", "3", "5", "6", "8", "7", "6"), CLI = c("7", "2", "5", "4", 
"4", "2", "4", "2", "3", "2", "8", "1", "1", "1", "1", "8", "8", 
"5", "7", "7", "8", "6", "5", "5", "6", "9", "8", "9", "6", "4", 
"2", "6", "9", "7", "7", "5", "9", "3", "7"), MAT = c("3", "7", 
"3", "7", "3", "5", "5", "8", "5", "4", "3", "5", "6", "4", "3", 
"1", "1", "1", "2", "2", "5", "3", "2", "2", "2", "5", "2", "2", 
"4", "3", "9", "2", "3", "2", "4", "4", "2", "4", "2"), EXP = c("2", 
"6", "8", "6", "6", "3", "3", "3", "7", "6", "7", "6", "8", "3", 
"5", "3", "2", "3", "1", "1", "1", "1", "1", "1", "1", "1", "1", 
"3", "8", "6", "4", "8", "2", "8", "2", "3", "3", "2", "4"), 
    CUL = c("4", "4", "2", "2", "8", "7", "8", "5", "6", "8", 
    "6", "8", "3", "9", "8", "9", "7", "7", "8", "9", "9", "9", 
    "9", "9", "9", "8", "9", "7", "1", "1", "1", "7", "8", "5", 
    "9", "7", "7", "9", "9"), IND = c("6", "5", "4", "3", "2", 
    "8", "6", "9", "8", "9", "2", "7", "4", "5", "9", "4", "6", 
    "8", "3", "4", "6", "2", "4", "6", "8", "2", "4", "6", "7", 
    "8", "6", "1", "1", "1", "1", "1", "1", "5", "8"), TST = c("9", 
    "8", "7", "8", "9", "4", "7", "6", "2", "7", "5", "4", "5", 
    "6", "4", "2", "3", "2", "4", "5", "4", "5", "7", "4", "4", 
    "3", "7", "4", "3", "7", "7", "3", "4", "4", "3", "8", "4", 
    "6", "3"), PHY = c("8", "9", "9", "9", "5", "9", "9", "7", 
    "9", "5", "9", "9", "7", "8", "6", "6", "4", "4", "5", "6", 
    "3", "4", "3", "3", "5", "6", "3", "5", "9", "9", "8", "9", 
    "6", "9", "6", "9", "6", "1", "1")), row.names = c(NA, 39L
), class = "data.frame")

roskam <- matrix(as.numeric(as.matrix(roskam)), 39, 9)

rosdata <- matrix(0, 351, 4)


k <- 1
for (j in 1:9) {
  for (i in 1:39) {
    rosdata[k, ] <- c(i, j, roskam[i, j], 1)
    k <- k + 1
  }
}
