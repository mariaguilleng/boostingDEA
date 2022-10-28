#' @title Data Pre-processing for Multivariate Adaptive Frontier Splines.
#'
#' @description This function arranges the data in the required format and
#' displays error messages.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables
#' in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param na.rm \code{logical} If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @importFrom stats na.omit
#'
#' @return It returns a \code{data.frame} in the required format.
preProcess <- function(data, x, y, na.rm = TRUE) {

  # x and y well / bad introduced

  cols <- 1:length(data)
  if (!(all(x %in% cols))) {
    stop("x index(es) are not in data.")

    if (!(all(y %in% cols))) {
      stop("y index(es) are not in data.")
    }
  }

  # Dataframe
  # List with variables
  # Matrix

  if (is.list(data) && !is.data.frame(data)) {

    # Data names?

    ifelse(is.null(names(data)),
      nms <- 1:length(data), # if not 1:x
      nms <- names(data)
    )

    data <- data.frame(matrix(unlist(data), ncol = length(nms), byrow = FALSE))
    names(data) <- nms
  } else if (is.matrix(data) || is.data.frame(data)) {
    data <- data.frame(data)
  }

  # Classes
  varClass <- unlist(sapply(data, class))

  # Output classes
  outClass <- varClass[y] %in% c("numeric", "double", "integer")

  # Error
  if (!all(outClass)) {
    stop(paste(names(data)[y][!outClass][1],
               "is not a numeric or integer vector"))
  }

  # Input classes
  # Ordered --> numeric
  for (i in x) {
    if (is.ordered(data[, i])) {
      data[, i] <- as.numeric(data[, i])
    }
  }

  # Define classes again
  varClass <- unlist(sapply(data, class))

  inpClass <- varClass[x] %in% c("numeric", "double", "integer")

  # Error
  if (!all(inpClass)) {
    stop(paste(names(data)[x][!inpClass][1],
               "is not a numeric, integer or ordered vector"))
  }

  data <- data[, c(x, y)]

  # NA values

  if (any(is.na(data))) {
    if (na.rm == TRUE) {
      data <- na.omit(data)
      warning("Rows with NA values have been omitted .\n")
    } else {
      stop("Please, detele or impute NA registers or set na.rm = TRUE to
           omit them. \n")
    }
  }

  return(data)
}
