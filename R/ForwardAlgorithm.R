#' @title Mean Squared Error
#'
#' @description This function computes the mean squared error between two
#' numeric vectors.
#'
#' @param y Vector of actual data.
#' @param yPred Vector of predicted values.
#'
#' @return Mean Squared Error.
mse <- function(y, yPred) {
  error <- sum((y - yPred)^2) / length(y)
  return(round(error, 4))
}


#' @title Add a new pair of Basis Functions
#'
#' @description This function adds the best pair of basis functions to the model
#'
#' @param data data \code{data.frame} or \code{matrix} containing the variables
#' in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param ForwardModel \code{list} containing the set of basis functions and the
#' B matrix.
#' @param knots_list \code{list} containing the set of selected knots.
#' @param Kp Maximum degree of interaction allowed.
#' @param minspan \code{integer}. Minimum number of observations between knots.
#' When \code{minspan = 0}, it is calculated as in Friedman's MARS paper section
#' 3.8 with alpha = 0.05.
#' @param Le \code{integer} Minimum number of observations before the first and
#' after the final knot.
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly
#' @param err_min Minimun error in the split.
#'
#' @importFrom dplyr %>%
#'
#' @return A \code{list} containing the matrix of basis functions (\code{B}), a
#' \code{list} of basis functions (\code{BF}), a \code{list} of selected knots
#' (\code{knots_list}) and the minimun error (\code{err_min}).
AddBF <- function(data, x, y, ForwardModel, knots_list, Kp, minspan, Le,
                  linpreds, err_min) {
  N <- nrow(data)
  nX <- length(x)

  # Set of basis functions
  BF_list <- ForwardModel[["BF"]]
  nBF <- length(BF_list)

  # Matrix of basis functions
  B <- ForwardModel[["B"]]

  signal <- 0 # improvement

  for (p in 1:nBF) {

    # Basis function for the division
    bf <- BF_list[[p]]

    # =================== #
    # Condition I: degree #
    # =================== #

    # If the number of variables exceeds the maximum allowed degree,
    # the division cannot be carried out by this basis function.
    if (length(bf[["xi"]]) >= Kp && !all(bf[["xi"]] == -1)) next

    for (xi in 1:nX) {

      # ===================== #
      # Condition II: product #
      # ===================== #

      # The same xi cannot appear twice in a product.
      if (any(bf[["xi"]] == xi)) next

      # Observations in the space.
      index <- bf[["Bp"]] != 0

      # ===================== #
      # Condition III: L & Le #
      # ===================== #

      # Set to FALSE first and last "Le" observations
      index[c(1:Le, (length(index) - Le + 1):length(index))] <- FALSE

      # Minspan
      if (minspan == 0) {
        L <- ceiling(-log2(-(1 / (nX * sum(index))) * log(1 - 0.05)) / 2.5)
      } else {
        L <- minspan
      }

      # Set to FALSE based on minspan
      if (!is.null(knots_list[[xi]])) {
        # Used indexes
        used_indexes <- sapply(knots_list[[xi]], "[[", "index")
        minspan_indexes <- c()

        for (idx in 1:length(used_indexes)) {
          # L observations of distance between knots
          ind <- c((used_indexes[idx] - L):(used_indexes[idx] + L))
          minspan_indexes <- append(minspan_indexes, ind)
        }

        # index must be greater than 0 and lower than N
        minspan_indexes <- unique(sort(minspan_indexes[minspan_indexes > 0][minspan_indexes <= N]))
        index[minspan_indexes] <- FALSE
      }

      # Knots
      knots <- data[index, xi] %>%
        unique() %>%
        sort()

      if (linpreds == TRUE) {
        knots <- c(0, knots)
      }

      if (length(knots) <= 1) next

      for (i in 1:length(knots)) {
        # Update B
        New_B <- CreateBF(data, xi, knots[i], B, p)

        # Estimate coefficients
        coefs <- EstimCoeffsForward(New_B, data[, y])

        # Predictions
        y_hat <- New_B %*% coefs

        # mse
        err <- mse(data[, y], y_hat)

        if (err < err_min) {

          # Model has improved
          signal <- 1
          err_min <- err
          Best_B <- New_B

          # index
          observation <- which(data[, xi] == knots[i])

          # New pair of basis functions
          bf1 <- bf2 <- bf

          # id
          bf1[["id"]] <- nBF + 1
          bf2[["id"]] <- nBF + 2

          # status
          bf1[["status"]] <- bf2[["status"]] <- "paired"

          # side
          bf1[["side"]] <- "R"
          bf2[["side"]] <- "L"

          # Bp
          bf1[["Bp"]] <- New_B[, ncol(New_B) - 1]
          bf2[["Bp"]] <- New_B[, ncol(New_B)]

          # xi
          if (all(bf[["xi"]] == -1)) {
            bf1[["xi"]] <- bf2[["xi"]] <- c(xi)
          } else {
            bf1[["xi"]] <- bf2[["xi"]] <- c(bf[["xi"]], xi)
          }

          # t
          bf1[["t"]] <- bf2[["t"]] <- knots[i]

          # R
          bf1[["R"]] <- bf2[["R"]] <- err_min

          # alpha
          bf1[["alpha"]] <- bf2[["alpha"]] <- coefs
        }
      }
    }
  }

  if (signal) {

    # Append new basis functions
    BF_list <- append(BF_list, list(bf1))
    BF_list <- append(BF_list, list(bf2))

    # bf2 & bf1 have the same xi and t
    sze_knt_xi <- length(bf1[["xi"]]) # in case of degree > 1 --> to select last xi
    knots_list[[bf1[["xi"]][sze_knt_xi]]] <- append(
      knots_list[[bf1[["xi"]][sze_knt_xi]]],
      list(list(t = bf1[["t"]], index = observation))
    )

    return(list(Best_B, BF_list, knots_list, err_min))
  } else {
    return(list(B, BF_list, knots_list, err_min))
  }
}

#' @title Generate a new pair of Basis Functions
#'
#' @description This function generates two new basis functions from a variable
#' and a knot.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in
#' the model.
#' @param xi \code{integer}. Variable index of the new basis function(s).
#' @param knt Knot for creating the new basis function(s).
#' @param B \code{matrix} of basis functions on which the new pair of functions
#' is added.
#' @param p \code{integer}. Parent basis function index.
#'
#' @return Matrix of basis function (\code{B}) updated with the new basis
#' functions.
CreateBF <- function(data, xi, knt, B, p) {

  # Create (xi-t)+ and (t-xi)+
  hinge1 <- ifelse(data[, xi] > knt, data[, xi] - knt, 0)
  hinge2 <- ifelse(data[, xi] < knt, knt - data[, xi], 0)

  # two new basis functions
  bf1 <- B[, p] * hinge1
  bf2 <- B[, p] * hinge2

  # update B
  B <- cbind(B, bf1, bf2)

  return(B)
}

#' @title Estimate Coefficients in Multivariate Adaptive Frontier Splines during
#' Forward Procedure.
#'
#' @description This function solves a Quadratic Programming Problem to obtain a
#' set of coefficients.
#'
#' @param B \code{matrix} of basis functions.
#' @param y Output \code{vector} in data.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return \code{vector} with the coefficients estimated.
EstimCoeffsForward <- function(B, y) {
  n <- nrow(B)
  p <- ncol(B)

  # vars: c(alpha_0, alpha_1, ..., alpha_P, e_1, ... , e_n)

  objVal <- c(rep(0, p), rep(1, n))

  # Structure for lpSolve
  # LP: lps <- make.lp(nrow = n + (p - 1) / 2 + p - 1, ncol = n + p)
  # LP: lp.control(lps, sense = 'min')
  # LP: set.objfn(lps, objVal)

  # = #
  # A #
  # = #
  # Equality constraints: y_hat - e = y
  Amat1 <- cbind(B, diag(rep(-1, n), n))

  # Concavity
  Amat2 <- matrix(0, nrow = (p - 1) / 2, ncol = p + n)

  pairs <- 1:p
  pairs <- pairs[lapply(pairs, "%%", 2) == 0]

  for (i in pairs) {
    Amat2[i / 2, c(i, i + 1)] <- -1
  }

  # Increasing monotony
  a0 <- rep(0, p - 1)
  alps <- diag(x = c(1, -1), p - 1)
  mat0 <- matrix(0, p - 1, n)

  Amat3 <- cbind(a0, alps, mat0)

  Amat <- rbind(Amat1, Amat2, Amat3)

  # = #
  # b #
  # = #
  bvec <- c(y, rep(0, (p - 1) / 2), rep(0, p - 1))

  # Direction of inequality
  # Rglpk
  dir <- c(rep("==", n), rep(">=", nrow(Amat) - n))

  # Bounds
  # Rglpk
  bounds <- list(lower = list(ind = 1L:p, val = rep(-Inf, p)))

  # LP
  # for(i in 1:n){
  # add.constraint(lps, xt = Amat[i, ], "=",  rhs = bvec[i])
  # }

  # for(i in (n + 1):nrow(Amat)){
  # add.constraint(lps, xt = Amat[i, ], ">=",  rhs = bvec[i])
  # }

  # set.bounds(lps, lower = c(rep(- Inf, p), upper = rep(0, n)))

  # Solve
  # Rglpk
  ans <- Rglpk_solve_LP(objVal, Amat, dir, bvec, bounds)

  alpha <- ans$solution[1:p]

  # LP: solve(lps)
  # LP: alpha <- get.variables(lps)[1:p]

  return(alpha)
}
