#' @title Linear programming model for radial output measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
BBC_out <- function(data, x, y, dataOriginal = data,
                    xOriginal = x, yOriginal = y, FDH = FALSE) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for (d in 1:numDMU) {
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
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 2:(j + 1), "binary")
    }
    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}



#' @title Linear programming model for radial input measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
BBC_in <- function(data, x, y, dataOriginal = data,
                   xOriginal = x, yOriginal = y, FDH = FALSE) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for (d in 1:numDMU) {
    objVal <- matrix(ncol = j + 1, nrow = 1)
    objVal[1] <- 1
    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + 1)
    lp.control(lps, sense = "min")
    set.objfn(lps, objVal)
    # constrain 2.1 and 2.2
    for (xi in 1:nX)
    {
      add.constraint(lps, xt = c(- x_k[d, xi], xOriginal_k[, xi]), "<=", rhs = 0)
    }
    for (yi in 1:nY)
    {
      add.constraint(lps, xt = c(0, yOriginal_k[, yi]), ">=", rhs = y_k[d, yi])
    }
    # Constrain 2.3 - phi = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, j)), type = "=", rhs = 1)
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 2:(j + 1), "binary")
    }
    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}


#' @title Linear programming model for Russell output measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
Russell_out <- function(data, x, y, dataOriginal = data,
                        xOriginal = x, yOriginal = y, FDH = FALSE) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for (d in 1:numDMU) {

    objVal <- matrix(ncol = j + nY, nrow = 1)
    objVal[1:nY] <- 1 / nY

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + nY)
    lp.control(lps, sense = "max")
    set.objfn(lps, objVal)
    # constrain 2.1 and 2.2
    for (xi in 1:nX)
    {
      add.constraint(lps, xt = c(rep(0,nY), xOriginal_k[, xi]), "<=", rhs = x_k[d, xi])
    }
    for (yi in 1:nY)
    {
      vec <- c()
      vec[yi] <- - y_k[d, yi]
      vec[(1:nY)[- yi]] <- 0
      vec[(nY + 1):(nY + j)] <- yOriginal_k[, yi]
      add.constraint(lps, xt = vec, ">=", rhs = 0)
    }
    # Constrain 2.3 - sum(lambdas) = 1
    add.constraint(lprec = lps, xt = c(rep(0, nY), rep(1, j)), type = "=", rhs = 1)
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 1:j + nY, "binary")
    }
    # Constrain 2.5 - phi_m >= 1
    set.bounds(lps, columns = 1:nY, lower = rep(1.0,nY), upper = rep(Inf,nY))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}


#' @title Linear programming model for Russell input measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
Russell_in <- function(data, x, y, dataOriginal = data,
                       xOriginal = x, yOriginal = y, FDH = FALSE) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for (d in 1:numDMU) {

    objVal <- matrix(ncol = j + nX, nrow = 1)
    objVal[1:nX] <- 1 / nX

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + nX)
    lp.control(lps, sense = "min")
    set.objfn(lps, objVal)
    # constrain 2.1 and 2.2
    for (xi in 1:nX)
    {
      vec <- c()
      vec[xi] <- - x_k[d, xi]
      vec[(1:nX)[- xi]] <- 0
      vec[(nX + 1):(nX + j)] <- xOriginal_k[, xi]

      add.constraint(lps, xt = vec, "<=",  rhs = 0)
    }
    for (yi in 1:nY)
    {
      add.constraint(lps, xt = c(rep(0, nX), yOriginal_k[, yi]), ">=", rhs = y_k[d, yi])
    }
    # Constrain 2.3 - sum(lambdas) = 1
    add.constraint(lprec = lps, xt = c(rep(0, nX), rep(1, j)), type = "=", rhs = 1)
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 1:j + nX, "binary")
    }
    # Constrain 2.5 -  0 <= phi_m <= 1
    set.bounds(lps, columns = 1:nX, lower = rep(0,nX), upper = rep(1.0,nX))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}


#' @title Linear programming model for Directional Distance Function measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#' @param direction.vector Direction vector. Valid values are: \code{dmu} (x_0, y_0),
#' \code{unit} (unit vector), \code{mean} (mean values of each variable) and
#' a user specific vector of the same length as the number of input and output
#' variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#' @importFrom methods is
#'
#' @return \code{matrix} with the the predicted score
DDF <- function(data, x, y, dataOriginal = data,
                xOriginal = x, yOriginal = y, FDH = FALSE, direction.vector) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  dmu <- FALSE # flag for direction vector as (x_0, y_0)
  g <- direction.vector
  if (is.null(direction.vector)) {
    stop("Direction vector g must be specified")
  }
  if (is(direction.vector,"numeric")) {
    if (length(direction.vector) != length(x) + length(y)) {
      stop("Length of direction vector must be equal to the number of input and
           output variables")
    }
  } else if (direction.vector == "unit") {
    g <- rep(1, length(x) + length(y))
  } else if ( direction.vector == "mean") {
    g <- sapply(data[,c(x,y)], mean)
  } else if (direction.vector == "dmu") {
    dmu <- TRUE
  } else {
    stop("Invalid value of direction vector")
  }


  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for(d in 1:numDMU){

    if(dmu) {
      g <- c()
      for(var in c(x,y)) {
        g[var] <- data[d, var]
      }
    }

    objVal <- matrix(ncol = 1 + j, nrow = 1) # beta + lambdas
    objVal[1] <- 1 # beta

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + 1)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX)
    { # beta[g-], a, <=, x
      add.constraint(lps, xt = c(g[xi], xOriginal_k[, xi]), "<=",  rhs = x_k[d, xi])
    }
    for(yi in 1:nY)
    { # - y, d(a), >=, beta[g+]
      add.constraint(lps, xt = c(- g[nX + yi], yOriginal_k[, yi]), ">=", rhs = y_k[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, j)), type = "=", rhs = 1)

    # Constrain 2.4
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 2:(j + 1), "binary")
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}


#' @title Linear programming model for Weighted Additive Model
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#' @param weights Weights. Valid values are: \code{MIP} (Measure of Inefficiency
#' Proportions), \code{RAM} (Range Adjusted Measure), \code{BAM} (Bounded
#' Adjusted Measure), \code{normalized} (normalized weighted additive model)
#' and a user specific vector of the same length as the number of input and
#' output variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#' @importFrom stats sd
#' @importFrom methods is
#'
#' @return \code{matrix} with the the predicted score
WAM <- function(data, x, y, dataOriginal = data,
                xOriginal = x, yOriginal = y, FDH = FALSE, weights) {

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  if (is.null(weights)) {
    stop("Weigths must be specified")
  }
  if (is(weights, "numeric")) {
    if (length(weights) != length(x) + length(y)) {
      stop("Length of direction vector must be equal to the number of input and
           output variables")
    }
  } else if (weights == "RAM") {
    inputRanges <- apply(xOriginal_k, 2, max) - apply(xOriginal_k, 2, min)
    outputRanges <- apply(yOriginal_k, 2, max) - apply(yOriginal_k, 2, min)
    ranges <- 1 / ((nX + nY) * c(inputRanges, outputRanges))
  } else if (weights == "BAM") {
    min_x <- apply(xOriginal_k, 2, min)
    max_y <- apply(yOriginal_k, 2, max)
  } else if (weights == "normalized") {
    input_sd <- apply(xOriginal_k, 2, sd)
    output_sd <- apply(yOriginal_k, 2, sd)
    sd_vector <- 1 / c(input_sd, output_sd)
  } else if (weights  != "MIP") {
    stop("Invalid value of weights")
  }

  # get scores
  for(d in 1:numDMU){

    objVal <- matrix(ncol = nX + nY + j, nrow = 1)

    if (is(weights,"numeric")) {
      objVal[1:(nX + nY)] <- weights
    } else if (weights == "MIP") {
      objVal[1:(nX + nY)] <- c(1 / x_k[d, ], 1 / y_k[d, ])
    } else if (weights == "RAM"){
      objVal[1:(nX + nY)] <- ranges
    } else if (weights == "BAM") {
      objVal[1:(nX + nY)] <- c(1 / ((nX+nY) * (x_k[d, ] - min_x)),
                               1 / ((nX+nY) * (max_y - y_k[d, ])))
      objVal <- replace(objVal, objVal==Inf, 0)
    } else if (weights == "normalized") {
      objVal[1:(nX + nY)] <- sd_vector
    }

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nX + nY + j)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX)
    {
      vec <- c()
      vec[xi] <- 1
      vec[(1:nX)[- xi]] <- 0
      vec[(nX + 1):(nX + nY)] <- 0
      vec[(nX + nY + 1):(nY + nX + j)] <- xOriginal_k[, xi]

      add.constraint(lps, xt = vec, "<=",  rhs = x_k[d, xi])
    }

    for(yi in 1:nY)
    {
      vec <- c()
      vec[1:nX] <- 0
      vec[nX + yi] <- - 1
      vec[((nX + 1):(nX + nY))[- yi]] <- 0
      vec[(nX + nY + 1):(nY + nX + j)] <- yOriginal_k[, yi]

      add.constraint(lps, xt = vec, ">=", rhs = y_k[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(rep(0, nY + nX), rep(1, j)),
                   type = "=", rhs = 1)

    # Constrain 2.4
    if (FDH) {
      set.type(lps, columns = 1:j + (nX + nY), type = c("binary"))
    }

    status <- solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}


#' @title Enhanced Russell Graph measure
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the new variables
#' in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#' @param dataOriginal \code{data.frame} or \code{matrix} containing the
#' original
#' variables used to create the model.
#' @param xOriginal Vector. Column input indexes in original data.
#' @param yOriginal Vector. Column output indexes in original data.
#' @param FDH Binary decision variables
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type
#' set.bounds get.objective
#'
#' @return \code{matrix} with the the predicted score
ERG <- function(data, x, y, dataOriginal = data,
                xOriginal = x, yOriginal = y, FDH = FALSE) {

  options(scipen=999)

  if (length(x) != length(xOriginal) || length(y) != length(yOriginal)) {
    stop("Size of inputs or outputs does not match original data sample")
  }

  # variables
  j <- nrow(dataOriginal)
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  xOriginal_k <- as.matrix(dataOriginal[, x])
  yOriginal_k <- as.matrix(dataOriginal[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = numDMU, ncol = 1)

  # get scores
  for (d in 1:numDMU) {

    objVal <- matrix(ncol = j + 1 + nX + nY, nrow = 1)
    objVal[1] <- 1
    objVal[(1:nX) + 1] <- c(- 1 / (nX * x_k[d, ]))

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY + 2, ncol = j + 1 + nX + nY)
    lp.control(lps, sense = "min")
    set.objfn(lps, objVal)
    # contraint 2.1
    # TODO: modificar -- esta mal
    vec <- c()
    vec[1] <- 1
    vec[(1:nX)+1] <- 0
    vec[(nX + 2):(nX + nY + 1)] <- c(1 / (nY * y_k[d, ]))
    vec[(nX + nY + 2):(nY + nX + j + 1)] <- 0

    add.constraint(lps, xt = vec, "=", rhs = 1)
    # constrain 2.2
    for (xi in 1:nX)
    {
      vec <- c()
      vec[1] <- - x_k[d, xi]
      vec[(1:(nX+nY)) + 1] <- 0
      vec[xi+1] <- 1
      vec[(nX + nY + 2):(nX + nY + 1 + j)] <- xOriginal_k[, xi]

      add.constraint(lps, xt = vec, "=",  rhs = 0)
    }
    # constrain 2.3
    for(yi in 1:nY)
    {
      vec <- c()
      vec[1] <- - y_k[d, yi]
      vec[(1:(nX+nY)) + 1] <- 0
      vec[1 + nX + yi] <- - 1
      vec[(nX + nY + 2):(nY + nX + j + 1)] <- yOriginal_k[, yi]

      add.constraint(lps, xt = vec, "=", rhs = 0)
    }
    # Constrain 2.4 - sum(lambdas) = beta
    add.constraint(lprec = lps, xt = c(-1, rep(0, nX + nY), rep(1, j)),
                   type = "=", rhs = 0)
    if (FDH) {
      # Constrain 2.4 - Binary
      set.type(lps, 1:j + (nX+nY+1), "binary")
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}
