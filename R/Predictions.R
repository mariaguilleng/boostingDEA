#' @title Model Prediction for Adapted Multivariate Adaptive Frontier Splines.
#'
#' @description This function predicts the expected output by a \code{MARS}
#' object.
#'
#' @param object A \code{MARSAdapted} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param class Model for prediction. \code{1} MARS Boost without smoothing
#' procedure.\code{2} MARS Boost with smoothing procedure..
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{data.frame} with the predicted values.
predict.MARSAdapted <- function(object, newdata, x, class = 1, ...) {
  train_names <- object[["data"]][["input_names"]]
  test_names <- names(newdata)[x]

  if (!identical(sort(train_names), sort(test_names))) {
    stop("Different variable names in training and test sets.")
  }

  # Select variables and reorder as in training data
  newdata <- newdata[, x, drop = FALSE][train_names]
  N <- nrow(newdata)
  y <- object[["data"]][["y"]]

  # B1
  B <- matrix(rep(1, N), nrow = N)

  # Prediction for normal (without smoothing)
  if (class == 1) {
    model <- object[["MARS.Forward"]]
    # Always in pairs
    knots <- model[["knots"]]

    if (!is.null(knots)) {
      for (i in 1:nrow(knots)) {
        xi <- knots[i, 1]
        t <- knots[i, 2]

        # Basis functions for the new data
        hinge1 <- ifelse(newdata[, xi] > t, newdata[, xi] - t, 0)
        hinge2 <- ifelse(newdata[, xi] < t, t - newdata[, xi], 0)

        # Update B
        B <- cbind(B, hinge1, hinge2)
      }
    }

    # Prediction for smooth
  } else if (class == 2) {
    model <- object[["MARS.Forward.Smooth"]]

    knots <- model[["knots"]]

    if (is.null(knots)) {
      # Only possible in Smooth MARS (from backward)
      B <- B
    } else {
      for (st in c("paired", "not paired")) {
        for (i in 1:length(model[["cubic_knots"]])) {

          # This variables does not split the data
          if (is.null(model[["cubic_knots"]][[i]])) next

          for (j in 1:length(model[["cubic_knots"]][[i]])) {
            if (model[["cubic_knots"]][[i]][[j]][["status"]] != st) next

            # Cubic Knot
            Ct <- model[["cubic_knots"]][[i]][[j]][["t"]]

            # Knot
            t <- Ct[2]

            # Sides of that knot
            side <- c("L", "R")

            # Basis functions for the new data
            B <- CreateCubicBF(newdata, i, Ct, B, side)
          }
        }
      }
    }
  } else {
    stop("Algorithm class not valid.")
  }

  if (is.null(knots)) {
    predictions <- as.data.frame(rep(model[["alpha"]], N))
  } else {
    predictions <- as.data.frame(B %*% model[["alpha"]])
  }
  names(predictions) <- paste(object[["data"]][["output_names"]], "_pred",
                              sep = "")

  return(predictions)
}


#' @title Model Prediction for Efficiency Analysis Trees.
#'
#' @description This function predicts the expected output by an \code{EAT} object.
#'
#' @param object An \code{EAT} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param ... further arguments passed to or from other methods.
#'
#' @importFrom dplyr %>%
#'
#' @return \code{data.frame} with the predicted values.
predict.EAT <- function(object, newdata, x, ...) {

  if (!inherits(object, "EAT")){
    stop(paste(deparse(substitute(object)), "must be an EAT object"))
  }

  train_names <- object[["data"]][["input_names"]]
  test_names <- names(newdata)[x]

  if (!identical(sort(train_names), sort(test_names))) {
    stop("Different variable names in training set and test set.")
  }

  # Select variables and reorder as in training data
  newdata <- newdata[, x, drop = FALSE][train_names]

  y <- object[["data"]][["y"]]

  tree <- object[["tree"]]

  predictions <- c()

  for(register in 1:nrow(newdata)){
    ti <- 1

    while (tree[[ti]][["SL"]] != -1) {
      if (newdata[register, ][tree[[ti]][["xi"]]] < as.data.frame(tree[[ti]][["s"]])) {
        ti <- posIdNode(tree, tree[[ti]][["SL"]])
      } else {
        ti <- posIdNode(tree, tree[[ti]][["SR"]])
      }
    }
    predictions <- append(predictions, unlist(tree[[ti]][["y"]]))
  }

  predictions <- matrix(predictions, ncol = length(y), byrow = T) %>%
    as.data.frame()

  names(predictions) <- paste(object[["data"]][["output_names"]],"_pred", sep = "")

  return(predictions)
}


#' @title Model Prediction for Boosted Multivariate Adaptive Frontier Splines
#'
#' @description This function predicts the expected output by a \code{MARSBoost}
#' object.
#'
#' @param object A \code{MARSBoost} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param class Model for prediction. \code{1} MARS Boost without smoothing.
#' \code{2} MARS Boost with smoothing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{data.frame} with the predicted values.
#'
#' @export
predict.MARSBoost <- function(object, newdata, x, class = 1, ...) {
  train_names <- object[["data"]][["input_names"]]
  test_names <- names(newdata)[x]

  if (!identical(sort(train_names), sort(test_names))) {
    stop("Different variable names in training and test sets.")
  }

  # Select variables and reorder as in training data
  newdata <- newdata[, x, drop = FALSE][train_names]
  N <- nrow(newdata)
  y <- object[["data"]][["y"]]

  num.iterations <- object[["control"]][["num.iterations"]]
  learning.rate <- object[["control"]][["learning.rate"]]
  MARS_models <- object[["MARS.models"]]
  f0 <- object[["f0"]]

  # Get predictions
  predictions <- matrix(rep(f0[1, 1], N), ncol = length(y), nrow = N)

  if (class %in% c(1, 2)) {
    for (it in 1:length(MARS_models)) {
      model_q <- MARS_models[[it]]
      pred_q <- predict(model_q, newdata, x, class)
      predictions <- predictions + learning.rate * pred_q
    }
  } else {
    stop("Class not valid. \n1 - Prediction without smoothing procedure\n2 -
         Prediction with smoothing procedure\n")
  }

  return(predictions)
}


#' @title Model prediction for EATBoost algorithm
#'
#' @description This function predicts the expected output by a \code{EATBoost}
#' object.
#'
#' @param object A \code{EATBoost} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{data.frame} with the predicted values.
#'
#' @export
predict.EATBoost <- function(object, newdata, x, ...) {

  train_names <- object[["data"]][["input_names"]]
  test_names <- names(newdata)[x]

  if (!identical(sort(train_names), sort(test_names))) {
    stop("Different variable names in training and test sets.")
  }

  # Select variables and reorder as in training data
  newdata <- newdata[, x, drop = FALSE][train_names]
  N <- nrow(newdata)
  y <- object[["data"]][["y"]]

  num.iterations <- object[["control"]][["num.iterations"]]
  learning.rate <- object[["control"]][["learning.rate"]]
  EAT.models <- object[["EAT.models"]]
  f0 <- object[["f0"]]

  # Get predictions
  init_pred <- as.numeric(as.vector(f0[1,]))
  prediction <- matrix(rep(init_pred,N), ncol = length(y), nrow = N,
                       byrow = TRUE)

  for (it in 1:length(EAT.models)) {
    model.q <- EAT.models[[it]]
    pred.q <- predict(model.q, newdata, x)
    prediction <- prediction + learning.rate*pred.q
  }

  return(prediction)
}


#' @title Model Prediction for DEA
#'
#' @description This function predicts the expected output by a \code{DEA}
#' object.
#'
#' @param object A \code{DEA} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param y Outputs index.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{data.frame} with the predicted values. Valid measures are:
#' \code{rad.out}.
#'
#' @export
predict.DEA <- function(object, newdata, x, y, ...) {
  scores <- BBC_out(newdata, x, y, object[["data"]][["df"]],
                    object[["data"]][["x"]], object[["data"]][["y"]])
  pred <- as.data.frame(scores * object[["data"]][["df"]][1:nrow(newdata),object[["data"]][["y"]]])
  names(pred) <- paste(object[["data"]][["output_names"]],"_pred", sep = "")
  return(pred)
}


#' @title Model Prediction for FDH
#'
#' @description This function predicts the expected output by a \code{FDH}
#' object.
#'
#' @param object A \code{FDH} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param y Outputs index.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{data.frame} with the predicted values. Valid measures are:
#' \code{rad.out}.
#'
#' @export
predict.FDH <- function(object, newdata, x, y, ...) {
  scores <- BBC_out(newdata, x, y, object[["data"]][["df"]],
                    object[["data"]][["x"]], object[["data"]][["y"]],
                    FDH = TRUE)
  pred <- as.data.frame(scores * object[["data"]][["df"]][1:nrow(newdata),object[["data"]][["y"]]])
  names(pred) <- paste(object[["data"]][["output_names"]],"_pred", sep = "")
  return(pred)
}


