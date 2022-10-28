#' @title Gradient Tree Boosting
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as LS-boosting using adapted Multivariate Adaptive Regression Splines (MARS) as base learners.
#'
#' @name EATBoost
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param num.iterations Maximum number of iterations the algorithm will perform
#' @param learning.rate Learning rate that control the overfitting of the algorithm. Value must be in (0,1]
#' @param num.leaves Maximum number of reflected pairs created by the forward algorithm of MARS.
#'
#' @importFrom stats predict
#' @importFrom eat EAT
#'
#' @return A \code{EATBoost} object.
#'
#' @export
EATBoost <- function(data, x, y, num.iterations, num.leaves, learning.rate) {
  if (!is.null(num.iterations) && num.iterations < 1) {
    stop("num.iterations = ", num.iterations, "not valid. Number of iterations must be greater than 1")
  }
  if (!is.null(learning.rate) && (learning.rate <= 0 || learning.rate > 1)) {
    stop("learning.rate = ", learning.rate, "not valid. Learning rate must in (0,1]")
  }
  if (!is.null(num.leaves) && num.leaves < 2) {
    stop("num.leaves = ", num.leaves, "not valid. Maximum number of leaves at each iteration must be greater than or equal to 2")
  }

  # ===========#
  # VARIABLES #
  # ===========#

  # Prepocess
  data <- preProcess(data, x, y)

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)
  nY <- length(y)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - nY)
  y <- (nX + 1):ncol(data)

  # pseudo-residuals
  residuals <- matrix(0, ncol = nY, nrow = nrow(data))

  # list of models created in each iterations
  EAT.models <- list()

  # ===========#
  # FIT MODEL #
  # ===========#

  browser()

  # initial prediction
  f0 <- matrix(rep(max(data[, y]), N), ncol = nY, nrow = nrow(data))
  prediction <- f0

  # prediction at each iteratrion
  for (it in 1:num.iterations) {

    # Calculate pseudo-residuals
    residuals <- data[, y] - prediction

    # Fit forward MARS to pseudo-residuals
    data.q <- as.data.frame(cbind(data[, x], residuals))
    colnames(data.q) <- colnames(data)
    model.q <- EAT(
      data = data.q,
      x = x,
      y = y,
      max.leaves = num.leaves
    )

    EAT.models[[it]] <- model.q

    # Update prediction
    prediction.q <- predict(model.q, data.q, x)
    prediction <- prediction + learning.rate * prediction.q
  }

  # EATBoost object
  EATBoost <- EATBoost_object(
    data, x, y, num.iterations, num.leaves,
    learning.rate, EAT.models, f0, prediction
  )

  return(EATBoost)
}


EATBoost_object <- function(data, x, y, num.iterations, num.leaves,
                            learning.rate, EAT.models, f0, prediction) {
  EATBoost_object <- list(
    "data" = list(
      df = data,
      x = x,
      y = y,
      input_names = names(data)[x],
      output_names = names(data)[y],
      row_names = rownames(data)
    ),
    "control" = list(
      num.iterations = num.iterations,
      num.leaves = num.leaves,
      learning.rate = learning.rate
    ),
    "EAT.models" = EAT.models,
    "f0" = f0,
    "prediction" = prediction
  )

  class(EATBoost_object) <- "EATBoost"

  return(EATBoost_object)
}
