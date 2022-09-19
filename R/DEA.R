#' @title Model prediction for DEA
#'
#' @description This function predicts the expected output through a DEA model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Vector. Column input indexes in data.
#' @param y Vector. Column output indexes in data.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @export
#'
#' @return \code{data.frame} with the original data and the predicted values through a DEA model.
predictDEA <- function(data,x,y) {
  
  # variables
  j <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_k <- as.matrix(data[, y])
  nX <- length(x)
  nY <- length(y)
  scores <- matrix(nrow = j, ncol = 1)
  
  # get scores
  for(d in 1:j){
    objVal <- matrix(ncol = j + 1, nrow = 1)
    objVal[1] <- 1
    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + 1)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)
    # constrain 2.1 and 2.2
    for(xi in 1:nX)
    {
      add.constraint(lps, xt = c(0, x_k[, xi]), "<=",  rhs = x_k[d, xi])
    }
    for(yi in 1:nY)
    {
      add.constraint(lps, xt = c(- y_k[d, yi], y_k[, yi]), ">=", rhs = 0)
    }
    # Constrain 2.3 - phi = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, j)), type = "=", rhs = 1)
    solve(lps)
    scores[d, ] <- get.objective(lps)
  }
  
  # get prediction
  pred_DEA <- scores * data[, y]
  
  result <- data.frame(pred = pred_DEA, score = scores)
  return(result)
}
