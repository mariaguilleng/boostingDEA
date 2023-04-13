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
efficiency <- function(model, measure = "rad.out", data, x, y, heuristic = TRUE,
                       direction.vector = NULL, weights = NULL) {

  # parameter control
  valid_measures <- c("rad.out", "rad.in", "Russell.out", "Russell.in", "DDF",
                      "WAM", "ERG")
  valid_models <- c("DEA", "FDH", "EATBoost", "MARSBoost")
  if (! measure %in% valid_measures) {
    stop("Measure not valid. Valid values for measure are: ", valid_measures)
  }
  model_class <- class(model)
  if (! model_class %in% valid_models) {
    stop("Model not valid. Valid classes for model are: ", valid_models)
  }

  if (measure != "WAM" && !is.null(weights)) {
    warning("Value of weights parameter will be ignored")
  }
  if (measure != "DDF" && !is.null(direction.vector)) {
    warning("Value of direction.vector parameter will be ignored")
  }

  # MARSBoost
  if (model_class == "MARSBoost"){
    if (measure == "rad.out") {
      score <- model$prediction / data[,y]
      df <- as.data.frame(score, row.names = rownames(data))
      colnames(df) <- c(paste(model_class, measure, sep = "."))
      return(df)
    } else {
      stop("Measure not valid for MARSBoost model. Only valid measure is rad.out")
    }
  }

  # data used to create the model
  fdh <- (model_class == "FDH" ||
          (model_class == "EATBoost" && !measure %in% c("WAM", "ERG")))
  xOriginal <- model[["data"]][["x"]]
  yOriginal <- model[["data"]][["y"]]
  dataOriginal <- model[["data"]][["df"]]
  if (model_class %in% c("DEA", "FDH")) {
    baseData <- dataOriginal
  } else if (model_class == "EATBoost" && heuristic == TRUE) {
    pred <- predict(model, dataOriginal, xOriginal)
    baseData <- cbind(dataOriginal[, xOriginal], pred)
  } else { # model_class == "EATBoost" && heuristic == FALSE
    cat("Calculating EATBoost", measure, "efficiency measure. This migth take a while...\n")
    final_a <- get.a.EATBoost(model)
    colnames(final_a) <- model[["data"]][["input_names"]]
    pred_a <- predict(model, final_a, 1:ncol(final_a))
    baseData <- cbind(final_a, pred_a)
  }

  # calculate score
  if (measure == "rad.out") {
    score <- BBC_out(data, x, y, baseData, xOriginal, yOriginal, FDH = fdh)
  } else if (measure == "rad.in") {
    score <- BBC_in(data, x, y, baseData, xOriginal, yOriginal, FDH = fdh)
  } else if (measure == "Russell.out") {
    score <- Russell_out(data, x, y, baseData, xOriginal, yOriginal, FDH = fdh)
  } else if (measure == "Russell.in") {
    score <- Russell_in(data, x, y, baseData, xOriginal, yOriginal, FDH = fdh)
  } else if (measure == "DDF") {
    score <- DDF(data, x, y, baseData, xOriginal, yOriginal,
                 direction.vector = direction.vector, FDH = fdh)
  } else if (measure == "WAM") {
    score <- WAM(data, x, y, baseData, xOriginal, yOriginal,
                 weights = weights, FDH = fdh)
  } else { # measure == ERG
    score <- ERG(data, x, y, baseData, xOriginal, yOriginal, FDH = fdh)
  }

  # return score
  df <- as.data.frame(score, row.names = rownames(data))
  if (model_class == "EATBoost" && heuristic) {
    colnames(df) <- c(paste(model_class, "heu", measure, sep = "."))
  } else {
    colnames(df) <- c(paste(model_class, measure, sep = "."))
  }

  return(df)
}

#' @title Get the inferior corner of the leave support from all trees
#' of \code{EATBoost}
#'
#' @description Calculates the inferior corner of the support of all leave nodes
#' of every tree created in the \code{EATBoost} model
#'
#' @param EATBoost_model Model from class \code{EATBoost} from which the data
#' are obtained
#'
#' @return \code{list} of \code{matrix}. The length of the list is equal to
#' the \code{num.iterations} of the \code{EATBoost_model}. Each \code{matrix}
#' corresponds to a tree,  where the number of columns is the number of input
#'  variables and the number of rows to the number of leaves
#'
get.a.trees <- function(EATBoost_model) {
  list_a = list()
  EAT.models <- EATBoost_model[["EAT.models"]]
  for (i in 1:length(EAT.models)) {
    list_a[[i]] <- EAT.models[[i]][["model"]][["a"]]
  }
  return(list_a)
}

#' @title Get the superior corner of the leave support from all trees
#' of \code{EATBoost}
#'
#' @description Calculates the superior corner of the support of all leave nodes
#' of every tree created in the \code{EATBoost} model
#'
#' @param EATBoost_model Model from class \code{EATBoost} from which the data
#' are obtained
#'
#' @return \code{list} of \code{matrix}. The length of the list is equal to
#' the \code{num.iterations} of the \code{EATBoost_model}. Each \code{matrix}
#' corresponds to a tree,  where the number of columns is the number of input
#'  variables and the number of rows to the number of leaves
#'
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

#' @title Get intersection between two leaves supports
#'
#' @description Calculates the intersection between two leave nodes from
#' different trees of a \code{EATBoost} model.
#'
#' @param comb_a_actual Inferior corner of first leave support
#' @param comb_b_actual Superior corner of first leave support
#'
#' @return \code{vector} with the intersection. \code{NULL} if intersection
#' is not valid.
#'
get.intersection.a <- function(comb_a_actual, comb_b_actual) {
  intersection_a <- apply(comb_a_actual, 2, max)
  intersection_b <- apply(comb_b_actual, 2, min)
  if (sum(intersection_a < intersection_b) == length(intersection_a)) {
    return(intersection_a)
  }
}

#' @title Get \code{EATBoost} leaves supports
#'
#' @description Calculates the inferior corner of the leaves supports of a
#' \code{EATBoost} model.
#'
#' @param EATBoost_model Model from class \code{EATBoost} from which the data
#' are obtained
#'
#' @import utils
#'
#' @return \code{data.frame} with the leave supports
#'
get.a.EATBoost <- function(EATBoost_model) {

  nX <- length(EATBoost_model[["data"]][["x"]])
  num.iterations <- EATBoost_model[["control"]][["num.iterations"]]
  num.leaves <- EATBoost_model[["control"]][["num.leaves"]]

  # list_a : matrix of size num.iterations x num.leaves x nX
  list_a <- get.a.trees(EATBoost_model)
  # list_b : matrix of size num.iterations x num.leaves x nX
  list_b <- get.b.trees(EATBoost_model)

  final_a = matrix(nrow = 0, ncol = nX)

  # indices de la ultima combinacion posible
  comb_index_max <- rep(num.leaves, num.iterations)
  # indices de la combinacion actual
  comb_index_actual <- rep(1, num.iterations)

  # progress bar
  max.iterartions <- num.leaves^num.iterations
  ii <- 0
  pb <- txtProgressBar(min = 0, max = max.iterartions, style = 3)

  # mientras que los indices de la actual no mayores que los de la ultima
  repeat {

    comb_a_actual <- matrix(nrow = 0, ncol = nX)
    comb_b_actual <- matrix(nrow = 0, ncol = nX)

    # hacer intersecciones
    for (i in 1:num.iterations) {
      comb_a_actual <- rbind(comb_a_actual, list_a[[i]][comb_index_actual[i],])
      comb_b_actual <- rbind(comb_b_actual, list_b[[i]][comb_index_actual[i],])
    }
    a <- get.intersection.a(comb_a_actual, comb_b_actual)
    final_a <- rbind(final_a, a)

    # progress bar
    ii = ii + 1
    setTxtProgressBar(pb, ii)

    # condicion de parada
    if (sum(comb_index_actual) == sum(comb_index_max)) {
      break
    }

    # actualizar indices
    for (i in 1:num.iterations) {
      comb_index_actual[i] = comb_index_actual[i]+1
      if (comb_index_actual[i] <= comb_index_max[i]) {
        break
      } else {
        comb_index_actual[i] = 1
      }
    }

  }
  # progress bar
  close(pb)

  final_a <- as.matrix(final_a[!duplicated(final_a),])
  rownames(final_a) <- NULL
  return(as.data.frame(final_a))
}





