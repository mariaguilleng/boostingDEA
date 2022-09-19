#' @title Adapted Multivariate Adaptive Frontier Splines
#'
#' @description Create an adapted version of Multivariate Adaptive Regression Splines (MARS) model to estimate a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity.
#'
#' @name MARSAdapted
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param nterms Maximum number of reflected pairs created by the forward algorithm of MARS.
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}. If it is set to \code{-1}, \code{GCV = RSS / n}.
#' @param err_red Minimum reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param minspan Minimum number of observations between knots. When \code{minspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param endspan Minimum number of observations before the first and after the final knot. When \code{endspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @return An \code{AdaptedMARS} object.
#'
#' @export
MARSAdapted <- function(data, x, y, nterms, Kp = 1, d = 2, err_red = 0.01,
                  minspan = 0, endspan = 0, linpreds = FALSE,
                  na.rm = TRUE) {

  # Data in data[x, y] format.
  data <- preProcess(data = data, x = x, y = y, na.rm = na.rm)

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # ================= #
  # FORWARD ALGORITHM #
  # ================= #

  # basis function
    #     id: index
    # status: intercept / paired / not paired
    #   side: E (entire) / R (right) / L (left)
    #     Bp: basis function
    #     xi: variable for splitting
    #      t: knot for splitting
    #      R: mean squared error between true data and predicted data (B %*% alpha)
    #  alpha: regression coefficients
  bf <- list(
    "id" = 1,
    "status" = "intercept",
    "side" = "E",
    "Bp" = rep(1, N),
    "xi" = c(-1),
    "t" = c(-1),
    "R" = mse(data[, y], matrix(rep(1, N)) %*% max(data[, y])),
    "alpha" = max(data[, y])
  )

  # Set of knots. It saves indexes of data used as knots.
  knots_list <- vector("list", nX)

  # Set of basis functions and B Matrix
  ForwardModel <- list(BF = list(bf),
                       B = matrix(rep(1, N)))

  # Endspan
  if (endspan == 0){
    Le <- ceiling(3 - log2(0.05 / nX))
  } else {
    Le <- endspan
  }

  # Error of first basis function
  err <- bf[["R"]]

  while(length(ForwardModel[["BF"]]) + 2 < nterms) {

    # Divide the node
    B_BF_knots_err <- AddBF(data, x, y, ForwardModel, knots_list,
                            Kp, minspan, Le, linpreds, err_min = err)

     ForwardModel[["B"]] <- B_BF_knots_err[[1]]
    ForwardModel[["BF"]] <- B_BF_knots_err[[2]]
              knots_list <- B_BF_knots_err[[3]]
                 New_err <- B_BF_knots_err[[4]]

    # New minimun error
    if (New_err < err * (1 - err_red)) {
      err <- New_err

    } else {
      break
    }
  }

  # ==
  # Forward MARS
  # ==

  xi <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["xi"]]))
  t <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["t"]]))
  alpha <- ForwardModel[["BF"]][[length(ForwardModel[["BF"]])]][["alpha"]]

  knotsForward <- unique(cbind(xi, t))

  MARS.Forward = list(
       BF = ForwardModel[["BF"]],
        B = ForwardModel[["B"]],
    knots = knotsForward,
    alpha = alpha)

  # ==
  # Smooth forward
  # ==
  MARS.Forward.Smooth <- MARSAdaptedSmooth(data,nX,knotsForward,data[,y])


  # MARSAdapted object
  MARSAdapted <- MARSAdapted_object(data, x, y, rownames(data), nterms, Kp, d, err_red, minspan, 
               endspan, na.rm, MARS.Forward, MARS.Forward.Smooth)

  return(MARSAdapted)
}


#' @title Create an MARSAdapted object
#'
#' @description This function saves information about the adapetd Multivariate Adapative Frontier Splines model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param rownames \code{string}. Data rownames.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}. If set to \code{-1}, \code{GCV = RSS / n}.
#' @param err_red Minimun reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param minspan Minimum number of observations between knots. When \code{minspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param endspan Minimum number of observations before the first and after the final knot. When \code{endspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' @param MARS.Forward The Multivariate Adaptive Frontier Splines model after applying the forward algorithm without the smoothing procedures
#' @param MARS.Forward.Smooth The Multivariate Adaptive Frontier Splines model after applying the forward algorithm after appling the smoothing procedure
#'
#' @return A \code{MARSAdapted} object.
MARSAdapted_object <- function(data, x, y, rownames, nterms, Kp, d, err_red, minspan, endspan, na.rm, MARS.Forward, MARS.Forward.Smooth) {

  MARSAdapted_object <- list("data" = list(df = data,
                                    x = x,
                                    y = y,
                                    input_names = names(data)[x],
                                    output_names = names(data)[y],
                                    row_names = rownames),
                       "control" = list(nterms = nterms,
                                        Kp = Kp,
                                        d = d,
                                        err_red = err_red,
                                        minspan = minspan,
                                        endspan = endspan,
                                        na.rm = na.rm),
                       "MARS.Forward" = MARS.Forward,
                       "MARS.Forward.Smooth" = MARS.Forward.Smooth)

  class(MARSAdapted_object) <- "MARSAdapted"

  return(MARSAdapted_object)

}

