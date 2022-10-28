#' @title Single Output Data Generation
#'
#' @description This function is used to simulate the data in a single output
#' scenario.
#'
#' @param N Sample size.
#' @param nX Number of inputs. Possible values: \code{1}, \code{2}, \code{3},
#' \code{4}, \code{5},\code{6}, \code{9}, \code{12} and \code{15}.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rexp
#'
#' @return \code{data.frame} with the simulated data.
#'
#' @export
CobbDouglas <- function(N, nX, seed = ) {
  if (!(nX %in% c(1, 2, 3, 4, 5, 6, 9, 12, 15))) {
    stop(paste(nX, "is not allowed"))
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  for (x in 1:nX) {
    data[, x] <- runif(n = N, min = 0, max = 1)
  }

  u <- rexp(n = N, rate = 1 / 3)

  if (nX == 1) {
    y <- data[, "x1"]**0.5
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 2) {
    y <- (data[, "x1"]**0.4) * (data[, "x2"]**0.1)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 3) {
    y <- (data[, "x1"]**0.3) * (data[, "x2"]**0.1) * (data[, "x3"]**0.1)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 4) {
    y <- (data[, "x1"]**0.3) * (data[, "x2"]**0.1) * (data[, "x3"]**0.08) *
      (data[, "x4"]**0.02)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 5) {
    y <- (data[, "x1"]**0.3) * (data[, "x2"]**0.1) * (data[, "x3"]**0.08) *
      (data[, "x4"]**0.01) * (data[, "x5"]**0.01)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 6) {
    y <- (data[, "x1"]**0.3) * (data[, "x2"]**0.1) * (data[, "x3"]**0.08) *
      (data[, "x4"]**0.01) * (data[, "x5"]**0.006) * (data[, "x6"]**0.004)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y
  } else if (nX == 9) {
    y <- (data[, "x1"]**0.3) * (data[, "x2"]**0.1) * (data[, "x3"]**0.08) *
      (data[, "x4"]**0.005) * (data[, "x5"]**0.004) * (data[, "x6"]**0.001) *
      (data[, "x7"]**0.005) * (data[, "x8"]**0.004) * (data[, "x9"]**0.001)
    data["y"] <- y * exp(-u)
    data["yD"] <- y
  } else if (nX == 12) {
    y <- (data[, "x1"]**0.2) * (data[, "x2"]**0.075) * (data[, "x3"]**0.025) *
      (data[, "x4"]**0.05) * (data[, "x5"]**0.05) * (data[, "x6"]**0.08) *
      (data[, "x7"]**0.005) * (data[, "x8"]**0.004) * (data[, "x9"]**0.001) *
      (data[, "x10"]**0.005) * (data[, "x11"]**0.004) * (data[, "x12"]**0.001)
    data["y"] <- y * exp(-u)
    data["yD"] <- y
  } else {
    y <- (data[, "x1"]**0.15) * (data[, "x2"]**0.025) * (data[, "x3"]**0.025) *
      (data[, "x4"]**0.05) * (data[, "x5"]**0.025) * (data[, "x6"]**0.025) *
      (data[, "x7"]**0.05) * (data[, "x8"]**0.05) * (data[, "x9"]**0.08) *
      (data[, "x10"]**0.005) * (data[, "x11"]**0.004) * (data[, "x12"]**0.001) *
      (data[, "x13"]**0.005) * (data[, "x14"]**0.004) * (data[, "x15"]**0.001)
    data["y"] <- y * exp(-u)
    data["yD"] <- y
  }

  return(data)
}
