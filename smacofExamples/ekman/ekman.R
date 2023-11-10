source(data(ekman, package = "smacof")
n <- attr(ekman, "Size")
ekmanRC <- as.vector((1 - ekman) ^ 3)
ekmanR <- as.matrix((1 - ekman) ^ 3)


