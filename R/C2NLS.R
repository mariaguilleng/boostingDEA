library(quadprog)
library(dplyr)

# =================================== #
# Concave nonparametric least squares #
# =================================== #
# Equivalente a DEA gracias a R4
CNLS_DEA <- function(data,x,y) {

  n <- dim(data)[1]
  nX <- length(x)

  # num_vars = n + n*nX + n
  # vars: c(alpha_0, ..., alpha_N, beta_1_1,...,beta_1_n,...,beta_n_1,...,beta_n_n, e_1, ... , e_n)
  colnames <- c(paste("alpha_", 1:n, sep = ""),
                paste(paste("beta_", sort(rep(1:n,nX)), sep = ""), rep(1:nX,n), sep = "_"),
                paste("ep_", 1:n, sep = ""))

  # ================= #
  # D: Quadratic part #
  # ================= #

  # Identity matrix (1 for epsilons)
  Dmat <- diag(rep(1e-12, 2*n+n*nX))

  # Near 0 for alpha and beta
  Dmat[(n+n*nX+1):(2*n+n*nX), (n+n*nX+1):(2*n+n*nX)] <- diag(rep(1, n), n)
  colnames(Dmat) <- colnames

  # ============== #
  # d: Linear part #
  # ============== #
  dvec <- rep(0, 2*n + n*nX)

  # ================================== #
  # A: matrix defining the constraints #
  # ================================== #

  # R1 --> alpha_i + beta_i * x_i + ep_i = y_i
  # n restricciones de tipo R1
  Amat1 <- cbind(diag(rep(1, n), n), matrix(rep(diag(0,n),nX),n), diag(rep(1, n), n))
  colnames(Amat1) <- colnames
  cont <- n+1
  for (i in 1:n) {
    for (j in 1:nX) {
      Amat1[i,cont] <- data[i,j]
      cont <- cont + 1
    }
  }

  # R2 --> alpha_i - alpha_h + beta_i * x_i - beta_h * x_i + ep_i - ep_h = 0
  # n * n restricciones tipo R2
  Amat2 <- c()
  for (i in 1:n) {
    Amat2_N <- cbind(diag(rep(1, n), n), matrix(rep(diag(0,n),nX),n), diag(rep(0, n), n) )
    Amat2_N[,i] <- Amat2_N[,i] - 1
    cont <- n+1
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k,cont] <- data[i,j]
        cont <- cont + 1
      }
    }
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k,n+j+(i-1)*nX] <- Amat2_N[k,n+j+(i-1)*nX] - data[i,j]
      }
    }
    Amat2 <- rbind(Amat2, Amat2_N)
  }

  # R3 --> b_i >= 0
  # n restricciones de tipo R3
  Amat3 <- cbind(matrix(data = 0, ncol = n, nrow = n*nX), diag(rep(1, n*nX), n*nX), matrix(data = 0, ncol = n, nrow = n*nX))

  # R4 --> - e_i >= 0
  # n restricciones de tipo R4
  Amat4 <- cbind(diag(rep(0, n), n), matrix(rep(diag(0,n),nX),n), diag(rep(-1, n), n))

  # n + n*n + n*nX + n restricciones
  Amat <- rbind(Amat1, Amat2, Amat3, Amat4)
  colnames(Amat) <- colnames

  # ================================ #
  # b: right size of the constraints #
  # ================================ #
  bvec <- c(data[,y], rep(0, n*n), rep(0, n*nX), rep(0, n))

  # ============================ #
  # Solve the quadratic problems #
  # ============================ #
  s <- solve.QP(Dmat = Dmat, dvec = dvec,
                Amat = t(Amat), bvec = bvec,
                meq = n)

  # Solutions
  alpha <- s$solution[1:n]
  beta <- matrix(nrow = n, ncol = nX)
  cont <- n + 1
  for (i in 1:n) {
    for (j in 1:nX) {
      beta[i,j] = s$solution[cont]
      cont <- cont + 1
    }
  }
  ep <- tail(s$solution, n)
  df <- data.frame(alpha = alpha, beta = beta, epsilon = ep)
  df <- df %>% mutate_if(is.numeric, round, digits=10)
  f <- alpha + rowSums(beta*data[,x])
  df <- cbind(f,df)
  return(df)
}


# ================================ #
# Corrected ordinary least squares #
# ================================ #
# Regresión a la media
CNLS <- function(data,x,y) {

  n <- dim(data)[1]
  nX <- length(x)

  # num_vars = n + n*nX + n
  # vars: c(alpha_0, ..., alpha_N, beta_1, ..., beta_n, e_1, ... , e_n)
  colnames <- c(paste("alpha_", 1:n, sep = ""),
                paste(paste("beta_", sort(rep(1:n,nX)), sep = ""), rep(1:nX,n), sep = "_"),
                paste("ep_", 1:n, sep = ""))

  # ================= #
  # D: Quadratic part #
  # ================= #

  # Identity matrix (1 for epsilons)
  Dmat <- diag(rep(1e-12, 2*n+n*nX))

  # Near 0 for alpha and beta
  Dmat[(n+n*nX+1):(2*n+n*nX), (n+n*nX+1):(2*n+n*nX)] <- diag(rep(1, n), n)
  colnames(Dmat) <- colnames

  # ============== #
  # d: Linear part #
  # ============== #
  dvec <- rep(0, 2*n + n*nX)

  # ================================== #
  # A: matrix defining the constraints #
  # ================================== #

  # R1 --> alpha_i + beta_i * x_i + ep_i = y_i
  # n restricciones de tipo R1
  Amat1 <- cbind(diag(rep(1, n), n), matrix(rep(diag(0,n),nX),n), diag(rep(1, n), n))
  colnames(Amat1) <- colnames
  cont <- n+1
  for (i in 1:n) {
    for (j in 1:nX) {
      Amat1[i,cont] <- data[i,j]
      cont <- cont + 1
    }
  }

  # R2 --> alpha_i - alpha_h + beta_i * x_i - beta_h * x_i + ep_i - ep_h = 0
  # n * n restricciones tipo R2
  Amat2 <- c()
  for (i in 1:n) {
    Amat2_N <- cbind(diag(rep(1, n), n), matrix(rep(diag(0,n),nX),n), diag(rep(0, n), n) )
    Amat2_N[,i] <- Amat2_N[,i] - 1
    cont <- n+1
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k,cont] <- data[i,j]
        cont <- cont + 1
      }
    }
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k,n+j+(i-1)*nX] <- Amat2_N[k,n+j+(i-1)*nX] - data[i,j]
      }
    }
    Amat2 <- rbind(Amat2, Amat2_N)
  }

  # R3 --> b_i >= 0
  # n restricciones de tipo R3
  Amat3 <- cbind(matrix(data = 0, ncol = n, nrow = n*nX), diag(rep(1, n*nX), n*nX), matrix(data = 0, ncol = n, nrow = n*nX))

  # n + n*n + n*nX
  Amat <- rbind(Amat1, Amat2, Amat3)
  colnames(Amat) <- colnames

  # ================================ #
  # b: right size of the constraints #
  # ================================ #
  bvec <- c(data[,y], rep(0, n*n), rep(0, n*nX))

  # ============================ #
  # Solve the quadratic problems #
  # ============================ #
  s <- solve.QP(Dmat = Dmat, dvec = dvec,
                Amat = t(Amat), bvec = bvec,
                meq = n)

  # Solutions
  alpha <- s$solution[1:n]
  beta <- matrix(nrow = n, ncol = nX)
  cont <- n + 1
  for (i in 1:n) {
    for (j in 1:nX) {
      beta[i,j] = s$solution[cont]
      cont <- cont + 1
    }
  }
  ep <- tail(s$solution, n)
  df <- data.frame(alpha = alpha, beta = beta, epsilon = ep)
  df <- df %>% mutate_if(is.numeric, round, digits=10)
  f <- alpha + rowSums(beta*data[,x])
  df <- cbind(f,df)
  return(df)
}


# ============================================= #
# Corrected Concave nonparametric least squares #
# ============================================= #
# StoNEZD
C2NLS <- function(data,x,y) {
  n <- nrow(data)
  nX <- length(x)
  df_CNLS <- CNLS(data,x,y)
  ep <- df_CNLS$epsilon - max(df_CNLS$epsilon)
  alpha <- df_CNLS$alpha + max(df_CNLS$epsilon)
  beta <- matrix(nrow = n, ncol = nX)
  for (i in 1:n) {
    for (j in 1:nX) {
      beta[i,j] = df_CNLS[i,j+2] # f,alpha,betas,nep
    }
  }
  f <- alpha + rowSums(beta*data[,x])
  C2NLS <- C2NLS_object(data,x,y,alpha, beta, ep, f)
  return(C2NLS)
}

C2NLS_object <- function(data, x, y, alpha, beta, epsilon, f) {

  C2NLS_object <- list("data" = list(data = data,
                                    x = x,
                                    y = y,
                                    input_names = names(data)[x],
                                    output_names = names(data)[y]),
                      "alpha" = alpha,
                      "beta" = beta,
                      "epsilon" = epsilon,
                      "f" = f)

  class(C2NLS_object) <- "C2NLS"

  return(C2NLS_object)
}

predict.C2NLS <- function(C2NLS_object, newdata, x) {

  # Check that both data samples (train and predict) have the same number of inputs
  train_x <- C2NLS_object[["data"]][["x"]]
  if (length(x) != length(train_x)) {
    stop("Different variable names in training and test sets.")
  }
  n_pred <- nrow(newdata)

  # Get alpha and betas from training
  n_train <- nrow(C2NLS_object[["data"]][["data"]])
  alpha <- C2NLS_object[["alpha"]]
  beta <- C2NLS_object[["beta"]]
  epsilon <- C2NLS_object[["epsilon"]]

  # Predict: f(x) = min {alpha_i + beta_i * x}
  f <- matrix(nrow = n_pred, ncol = 1)
  for (k in 1:n_pred) {
    min = Inf
    for (i in 1:n_train) {
      if (length(train_x) == 1) { # mono-input
        comb_actual = alpha[k] + beta[k]*newdata[k,x]
      } else { # multi-input
        comb_actual = alpha[k] + rowSums(beta[k,]*newdata[k,x])
      }
      if (comb_actual < min) {
        min <- comb_actual
      }
    }
    f[k] = min
  }

  return(data.frame(x = newdata[,x], f = f))
}


