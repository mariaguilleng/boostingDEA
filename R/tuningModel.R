#' @title Tuning an MARSBoost model
#'
#' @description This funcion computes the root mean squared error (RMSE) for a
#' set of MARSBoost models built with a grid of given hyperparameters.
#'
#' @param training Training \code{data.frame} or \code{matrix} containing the
#'  variables for model construction.
#' @param test Test \code{data.frame} or \code{matrix} containing the variables
#' for model assessment.
#' @param x Column input indexes in \code{training}.
#' @param y Column output indexes in \code{training}.
#' @param num.iterations Maximum number of iterations the algorithm will perform
#' @param learning.rate Learning rate that control overfitting of the algorithm.
#'  Value must be in (0,1]
#' @param num.terms Maximum number of reflected pairs created by the forward
#' algorithm of MARS.
#' @param verbose Controls the verbosity.
#'
#' @importFrom dplyr arrange %>%
#' @importFrom MLmetrics MSE RMSE
#'
#' @return A \code{data.frame} with the sets of hyperparameters and the root
#' mean squared error (RMSE) associated for each model.
#'
#' @export
bestMARSBoost <- function(training, test, x, y, num.iterations, learning.rate,
                          num.terms, verbose = TRUE) {

  training <- preProcess(training, x, y)
  test <- preProcess(test, x, y)

  if (!identical(sort(names(training)), sort(names(test)))) {
    stop("Different variable names in training and test set")
  }

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(training) - length(y))
  y <- (length(x) + 1):ncol(training)

  # Grid of hyperparameters
  hp <- expand.grid(num.iterations = num.iterations,
                    learning.rate = learning.rate,
                    num.terms = num.terms,
                    RMSE = NA,
                    MSE = NA)

  num.combinations <- nrow(hp)
  for (i in 1:num.combinations) {

    if (verbose) {
      cat("Trying combination ", i, "/", num.combinations)
    }

    MARSBoost.model <- MARSBoost(data = training, x = x, y = y,
                                 num.iterations = hp[i, "num.iterations"],
                                 learning.rate = hp[i, "learning.rate"],
                                 num.terms = hp[i, "num.terms"])

    # RMSE
    pred <- predict(MARSBoost.model, test, x)
    hp[i, "RMSE"] <- RMSE(pred$y_pred, test[,y])
    hp[i, "MSE"] <- MSE(pred$y_pred, test[,y])

    if (verbose) {
      cat(" -- RMSE : ", round(hp[i, "RMSE"],4), "\n")
    }

  }

  hp <- hp %>% arrange(RMSE)

  return(hp)
}


#' @title Tuning an EATBoost model
#'
#' @description This funcion computes the root mean squared error (RMSE) for a
#' set of EATBoost models built with a grid of given hyperparameters.
#'
#' @param training Training \code{data.frame} or \code{matrix} containing the
#'  variables for model construction.
#' @param test Test \code{data.frame} or \code{matrix} containing the variables
#' for model assessment.
#' @param x Column input indexes in \code{training}.
#' @param y Column output indexes in \code{training}.
#' @param num.iterations Maximum number of iterations the algorithm will perform
#' @param learning.rate Learning rate that control overfitting of the algorithm.
#'  Value must be in (0,1]
#' @param num.leaves Maximum number of terminal leaves in each tree at each
#' iteration
#' @param verbose Controls the verbosity.
#'
#' @importFrom dplyr arrange %>%
#' @importFrom MLmetrics MSE RMSE
#'
#' @return A \code{data.frame} with the sets of hyperparameters and the root
#' mean squared error (RMSE) and mean squere error (MSE) associated for each
#' model.
#'
#' @export
bestEATBoost <- function(training, test, x, y, num.iterations, learning.rate,
                         num.leaves, verbose = TRUE) {

  training <- preProcess(training, x, y)
  test <- preProcess(test, x, y)

  if (!identical(sort(names(training)), sort(names(test)))) {
    stop("Different variable names in training and test set")
  }

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(training) - length(y))
  y <- (length(x) + 1):ncol(training)

  # Grid of hyperparameters
  hp <- expand.grid(num.iterations = num.iterations,
                    learning.rate = learning.rate,
                    num.leaves = num.leaves,
                    RMSE = NA,
                    MSE = NA)

  num.combinations <- nrow(hp)
  for (i in 1:num.combinations) {

    if (verbose) {
      cat("Trying combination ", i, "/", num.combinations)
    }

    EATBoost.model <- EATBoost(data = training, x = x, y = y,
                               num.iterations = hp[i, "num.iterations"],
                               learning.rate = hp[i, "learning.rate"],
                               num.leaves = hp[i, "num.leaves"])

    # RMSE
    pred <- predict(EATBoost.model, test, x)
    hp[i, "RMSE"] <- RMSE(pred$y_pred, test[,y])
    hp[i, "MSE"] <- MSE(pred$y_pred, test[,y])

    if (verbose) {
      cat(" -- RMSE : ", round(hp[i, "RMSE"],4), "\n")
    }

  }

  hp <- hp %>% arrange(RMSE)

  return(hp)
}
