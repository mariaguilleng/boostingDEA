#' @title Data Envelope Analysis model
#'
#' @description This function estimates a production frontier satisfying Data
#' Envelope Analysis axioms using the radial output measure.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in
#' the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return A \code{DEA} object.
#'
#' @export
DEA <- function(data, x, y) {

  scores <- BBC_out(data, x, y)
  pred_DEA <- scores * as.data.frame(data[,y], row.names = rownames(data))
  names(pred_DEA) <- paste(colnames(data)[y], "_pred", sep = "")


  # DEA object
  DEA <- DEA_object(data, x, y, pred_DEA, scores)

  return(DEA)
}


#' @title Create a DEA object
#'
#' @description This function saves information about the DEA model.
#'
#' @name DEA
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in
#' the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param pred Output predictions using the BBC radial output measure
#' @param score Efficiency score using the BBC radial output measure
#'
#' @return A \code{DEA} object.
#'
#' @export
DEA_object <- function(data, x, y, pred, score) {
  DEA_object <- list(
    "data" = list(
      df = data,
      x = x,
      y = y,
      input_names = names(data)[x],
      output_names = names(data)[y],
      dmu_names = rownames(data)
    ),
    "prediction" = pred
  )

  class(DEA_object) <- "DEA"

  return(DEA_object)
}
