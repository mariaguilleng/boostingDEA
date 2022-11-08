#' @title Calculate efficiency scores
#'
#' @description Calculates the efficiency score corresponding to the given model
#' using the given measure
#'
#' @param model Model object for which efficiency score is computed. Valid classes
#' are: \code{DEA}, \code{FDH}, \code{EATBoost} and \code{MARSBoost}.
#' @param measure Efficiency measure used. Valid measures are: \code{rad.out},
#' \code{rad.in}
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param heuristic Only used if \code{model} is \code{EATBoost}. This indicates
#' whether the heuristic or the exact approach is used.
#' @param direction.vector Only used when \code{measure} is \code{DDF}.Direction vector.
#' Valid values are: \code{dmu} (x_0, y_0), \code{unit} (unit vector),
#' \code{mean} (mean values of each variable) and a user specific vector of
#' the same length as the number of input and output variables
#' @param weights Only used when \code{measure} is \code{WAM}. Weights.
#' Valid values are: \code{MIP} (Measure of Inefficiency
#' Proportions), \code{RAM} (Range Adjusted Measure), \code{BAM} (Bounded
#' Adjusted Measure), \code{normalized} (normalized weighted additive model)
#' and a user specific vector of the same length as the number of input and
#' output variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
#'
#' @export
efficiency <- function(model, measure, data, x, y, heuristic = FALSE,
                       direction.vector, weights) {

  # parameter control
  valid_measures <- c("rad.out", "rad.in", "Russell.out", "Russell.in", "DDF",
                      "WAM")
  valid_models <- c("DEA", "FDH", "EATBoost", "MARSBoost")
  if (! measure %in% valid_measures) {
    stop("Measure not valid. Valid values for measure are: ", valid_measures)
  }
  model_class <- class(model)
  if (! model_class %in% valid_models) {
    stop("Model not valid. Valid classes for model are: ", valid_models)
  }

  # data used to create the model
  dataOriginal <- model[["data"]][["df"]]
  xOriginal <- model[["data"]][["x"]]
  yOriginal <- model[["data"]][["y"]]

  # calculate score
  switch(model_class,
         # DEA model
         DEA = {
           if (measure == "rad.out") {
             score <- BBC_out(data, x, y, dataOriginal, xOriginal, yOriginal)
           } else if (measure == "rad.in") {
             score <- BBC_in(data, x, y, dataOriginal, xOriginal, yOriginal)
           } else if (measure == "Russell.out") {
             score <- Russell_out(data, x, y, dataOriginal, xOriginal, yOriginal)
           } else if (measure == "Russell.in") {
             score <- Russell_in(data, x, y, dataOriginal, xOriginal, yOriginal)
           } else if (measure == "DDF") {
             score <- DDF(data, x, y, dataOriginal, xOriginal, yOriginal,
                          direction.vector = direction.vector)
           } else if (measure == "WAM") {
             score <- WAM(data, x, y, dataOriginal, xOriginal, yOriginal,
                          weights = weights)
           }
         },
         # FDH model
         FDH = {
           if (measure == "rad.out") {
             score <- BBC_out(data, x, y, dataOriginal, xOriginal, yOriginal,
                              FDH = TRUE)
           } else if (measure == "rad.in") {
             score <- BBC_in(data, x, y, dataOriginal, xOriginal, yOriginal,
                              FDH = TRUE)
           } else if (measure == "Russell.out") {
             score <- Russell_out(data, x, y, dataOriginal, xOriginal, yOriginal,
                                  FDH = TRUE)
           } else if (measure == "Russell.in") {
             score <- Russell_in(data, x, y, dataOriginal, xOriginal, yOriginal,
                                 FDH = TRUE)
           } else if (measure == "DDF") {
             score <- DDF(data, x, y, dataOriginal, xOriginal, yOriginal,
                          FDH = TRUE, direction.vector = direction.vector)
           } else if (measure == "WAM") {
             score <- WAM(data, x, y, dataOriginal, xOriginal, yOriginal,
                          FDH = TRUE, weights = weights)
           }
         },
         # EATBoost model
         EATBoost = {
           # heuristic approach
           if(heuristic == TRUE) {
             pred <- predict(model, dataOriginal, xOriginal)
             baseData <- cbind(dataOriginal[, xOriginal], pred)
           #exact scores
           } else {
             stop("Only heuristic approaach")
             # get a
             # get f(a)
             # baseData <- a + f(a)
           }
           if (measure == "rad.out") {
             score <- BBC_out(data, x, y, baseData, xOriginal, yOriginal)
           } else if (measure == "rad.in") {
             score <- BBC_in(data, x, y, baseData, xOriginal, yOriginal)
           } else if (measure == "Russell.out") {
             score <- Russell_out(data, x, y, baseData, xOriginal, yOriginal)
           } else if (measure == "Russell.in") {
             score <- Russell_in(data, x, y, baseData, xOriginal, yOriginal)
           } else if (measure == "DDF") {
             score <- DDF(data, x, y, dataOriginal, xOriginal, yOriginal,
                          direction.vector = direction.vector)
           } else if (measure == "WAM") {
             score <- WAM(data, x, y, dataOriginal, xOriginal, yOriginal,
                          weights = weights)
           }
         }
  )
  df <- as.data.frame(score, row.names = model[["data"]][["row_names"]])
  colnames(df) <- c(paste(model_class, measure, sep = "."))
  return(df)
}


get.a.trees <- function(EATBoost_model) {
  list_a = list()
  EAT.models <- EATBoost_model[["EAT.models"]]
  for (i in 1:length(EAT.models)) {
    list_a[[i]] <- EAT.models[[i]][["model"]][["a"]]
  }
  return(list_a)
}

get.b.trees <- function(EATBoost_model) {
  list_b = list()
  EAT.models <- EATBoost_model[["EAT.models"]]
  for (i in 1:length(EAT.models)) {
    b <- matrix(ncol = length(EATBoost_model[["data"]][["x"]]),
                nrow = EATBoost_model[["control"]][["num.leaves"]])
    leave <- 1
    for (j in 1:length(EAT.models[[i]][["tree"]])) {
      if (EAT.models[[i]][["tree"]][[j]][["SL"]] == -1) {
        b[leave,] <- EAT.models[[i]][["tree"]][[j]][["b"]]
        leave <- leave + 1
      }
    }
    list_b[[i]] <- b
  }
  return(list_b)
}

#' @export
get.a.EATBoost <- function(EATBoost_model) {

  list_a <- get.a.trees(EATBoost_model)
  list_b <- get.b.trees(EATBoost_model)

  final_a = matrix(nrow = 0, ncol = length(EATBoost_model[["data"]][["x"]]))

  num.iterations <- EATBoost_model[["control"]][["num.iterations"]]
  num.leaves <- EATBoost_model[["control"]][["num.leaves"]]

  for (i in 1:(num.iterations-1)) {

    # seleccionar a actual
    for (j in 1:num.leaves) {
      current_a <- list_a[[i]][j,]
      current_b <- list_b[[i]][j,]

      # comparar con el resto de a
      for (k in (i+1):num.iterations) {
        for (m in 1:num.leaves) {
          new_a <-list_a[[k]][m,]
          new_b <-list_b[[k]][m,]
          cat("--- comparando arbol", i, "nodo", j, "con arbol", k, "nodo", m, "\n")

          # intersection_a <- get.intersection.a(current_a, current_b, new_a, new_b)
          # TODO: calcular la interseccion como max {a} y min {b} y comprobar
          # si es correcta ( a < b en todas las coordenadas)

          # rbind(final_a, intersection_a)
          # si la interseccion es correcta, añadir a la solución
        }
      }
    }
  }
}





