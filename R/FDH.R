#' @title Free Disposal Hull model
#'
#' @description This function estimates a production frontier satisfying Free
#' Disposal HUll axioms using the radial output measure.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in
#' the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return A \code{FDH} object.
#'
#' @export
FDH <- function(data, x, y) {
  result <- predictFDH_BBC_out(data, x, y)

  # FDH object
  FDH <- FDH_object(data, x, y, result$pred, result$score)

  return(FDH)
}


#' @title Create a FDH object
#'
#' @description This function saves information about the FDH model.
#'
#' @name FDH
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in
#' the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param pred Output predictions using the BBC radial output measure
#' @param score Efficiency score using the BBC radial output measure
#'
#' @return A \code{FDH} object.
#'
#' @export
FDH_object <- function(data, x, y, pred, score) {
  FDH_object <- list(
    "data" = list(
      df = data,
      x = x,
      y = y,
      input_names = names(data)[x],
      output_names = names(data)[y],
      row_names = rownames(data)
    ),
    "pred" = pred,
    "score" = score
  )

  class(FDH_object) <- "FDH"

  return(FDH_object)
}



#' @title Model prediction for FDH
#'
#' @description This function predicts the expected output through a FDH model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{data.frame} with the the predicted values through a FDH model
#' and the efficiency score
predictFDH_BBC_out <- function(data, x, y, dataOriginal = data,
                               xOriginal = x, yOriginal = y) {
  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = j, ncol = 1)

  # get scores
  for (d in 1:j) {
    objVal <- matrix(ncol = j + 1, nrow = 1)
    objVal[1] <- 1
    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + 1)
    lp.control(lps, sense = "max")
    set.objfn(lps, objVal)
    # constrain 2.1 and 2.2
    for (xi in 1:nX)
    {
      add.constraint(lps, xt = c(0, xOriginal_k[, xi]), "<=", rhs = x_k[d, xi])
    }
    for (yi in 1:nY)
    {
      add.constraint(lps, xt = c(-y_k[d, yi], yOriginal_k[, yi]), ">=", rhs = 0)
    }
    # Constrain 2.3 - phi = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, j)), type = "=", rhs = 1)
    # Constrain 2.4 - Binary
    set.type(lps, 2:j + 1, "binary")
    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  # get prediction
  pred_FDH <- scores * data[, y]

  result <- data.frame(pred = pred_FDH, score = scores)
  return(result)
}
