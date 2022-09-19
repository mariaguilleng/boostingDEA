#' @title LS-Boosting with adapted Multivariate Adaptive Frontier Splines (MARS)
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as LS-boosting using adapted Multivariate Adaptive Regression Splines (MARS) as base learners.
#'
#' @name MARSBoost
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param num.iterations Maximum number of iterations the algorithm will perform
#' @param learning.rate Learning rate that control the overfitting of the algorithm. Value must be in (0,1]
#' @param num.terms Maximum number of reflected pairs created by the forward algorithm of MARS.
#'
#' @return A \code{MARSBoost} object.
#'
#' @export
MARSBoost <- function(data, x, y, num.iterations, learning.rate, num.terms) {

  
  if (!is.null(num.iterations) && num.iterations < 1)
    stop("num.iterations = ", num.iterations, "not valid. Number of iterations must be greater than 1")
  if (!is.null(learning.rate) && (learning.rate <= 0 || learning.rate > 1))
    stop("learning.rate = ", learning.rate, "not valid. Learning rate must in (0,1]")
  if (!is.null(num.terms) && num.terms < 2) 
    stop("num.terms = ", num.terms, "not valid. Maximum number of reflected pairs at each iteration must be greater than or equal to 2")
  
  
  #===========#
  # VARIABLES #
  #===========#

  # Data in data[x, y] format.
  data <- preProcess(data = data, x = x, y = y, na.rm = na.rm)

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # pseudo-residuals
  residuals <- matrix(0,ncol = length(y), nrow = nrow(data))
  residuals.smooth <- matrix(0,ncol = length(y), nrow = nrow(data))

  # list of models created in each iterations
  MARS.models <- list()

  #===========#
  # FIT MODEL #
  #===========#

  # initial prediction
  f0 <- matrix(rep(max(data[, y]),N),ncol = length(y), nrow = nrow(data))
  prediction <- f0
  prediction.smooth <- f0

  # prediction at each iteratrion
  for (it in 1:num.iterations) {

    # Calculate pseudo-residuals
    residuals <- data[, y] - prediction
    residuals.smooth <- data[, y] - prediction.smooth

    # Fit forward MARS to pseudo-residuals
    data_q <- as.data.frame(cbind(data[,x],residuals))
    colnames(data_q) <- colnames(data)
    model_q <- MARSAdapted(
      data = data_q,
      x = x,
      y = y,
      nterms = num.terms)

    MARS.models[[it]] <- model_q

    # Update prediction
    predictions_data_q <- predict(model_q, data_q, x)
    prediction <- prediction + learning.rate*predictions_data_q

    predictions_data_q_smooth <- predict(model_q, data_q, x, class = 2)
    prediction.smooth <- prediction.smooth + learning.rate*predictions_data_q_smooth

    if (((ncol(model_q[["MARS.Forward"]][["B"]]) - 1)/2) == 0) {
      break
    }

  }

  # MARSBoost object
  MARSBoost <- MARSBoost_object(data,x,y,num.iterations, learning.rate, num.terms,
                               MARS.models, f0, prediction, prediction.smooth)

  return(MARSBoost)

}


#' @title Create a MARSBoost object
#'
#' @description This function saves information about the LS-Boosted Multivariate Adapative Frontier Splines model.
#'
#' @name MARSBoost
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param num.iterations Maximum number of iterations the algorithm will perform
#' @param learning.rate Learning rate that control overfitting of the algorithm. Value must be in (0,1]
#' @param num.terms Maximum number of reflected pairs created by the forward algorithm of MARS.
#' @param f0 Initial predictions of the model (they correspond to maximum value of each output variable)
#' @param MARS.models List of the adapted forward MARS models created in each iterations
#' @param prediction Final predictions of the original data without applying the smoothing procedure
#' @param prediction.smooth Final predictions of the original data after applying the smoothing procedure
#'
#' @return A \code{MARSBoost} object.
#'
#' @export
MARSBoost_object <- function(data,x,y,num.iterations, learning.rate, num.terms,
                             MARS.models, f0, prediction, prediction.smooth) {

  MARSBoost_object <- list("data" = list(df = data,
                                    x = x,
                                    y = y,
                                    input_names = names(data)[x],
                                    output_names = names(data)[y],
                                    row_names = rownames(data)),
                      "control" = list(num.terms = num.terms,
                                       num.iterations = num.iterations,
                                       learning.rate = learning.rate),
                      "f0" = f0,
                      "MARS.models" = MARS.models,
                      "prediction" = prediction,
                      "prediction.smooth" = prediction.smooth)

  class(MARSBoost_object) <- "MARSBoost"

  return(MARSBoost_object)

}


