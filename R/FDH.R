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

  scores <- BBC_out(data, x, y, FDH = TRUE)
  pred <- scores * data[, y]

  # FDH object
  FDH <- FDH_object(data, x, y, pred, scores)

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
