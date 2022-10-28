library(dplyr)
library(quadprog)
library(lpSolveAPI)
library(tictoc)
library(ggplot2)
library(compare)
library(Rglpk)
library(eat)

# ==
# Generation of the data
# ==
nX <- 1
x <- 1:nX
y <- nX+1
N <- 100
seed <- 1234
data <- CobbDouglas(N,nX)

nterms <- 5
num_iterarions <- 10
learning_rate <- 0.5

ggplot(data) +
  geom_point(aes(x = x1, y = y, colour = 'y')) +
  geom_point(aes(x = x1, y = yD, colour = 'yD'))

# ==
# MAFSBoost
# ==
model_boost <- MARSBoost(data, x, y,
                         num_iterarions, learning_rate, nterms)
# without smoothin
pred_boost <- predict(model_boost, data, x, class = 1)
compare(pred_boost,model_boost$prediction)

# with smoothing
pred_boost_smooth <- predict(model_boost, data, x, class = 2)
compare(pred_boost_smooth,model_boost$prediction.smooth)

mse(data$yD, pred_boost)
mse(data$yD, pred_boost_smooth)

# ==
# DEA
# ==
DEA_model <- DEA(data,x,y)
pred_DEA <- predict(DEA_model, data, x, y)
compare(pred_DEA,DEA_model$pred)
mse(data$yD, pred_DEA)

# ==
# FDH
# ==
FDH_model <- FDH(data,x,y)
pred_FDH <- predict(FDH_model, data, x, y)
compare(pred_FDH,FDH_model$pred)
mse(data$yD, pred_FDH)

# ==
# EATBoost
# ==
# eat_model <- EAT(data, x, y)
# eatboost_model <- EATBoost(data, x, y, num_iterarions, nterms, learning_rate)
# pred_eat <- predict(eat_model, data, x)
#
# eatboost_model$prediction
# eatboost_model$f0
#
# data.q <- read.csv("data_q.csv", colClasses=c("NULL","double","double"))
#
# eat_model <- EAT(data.q, 1, 2, max.leaves = 5)
# predict(eat_model, data.q,1)


# =========
# EAT
# =========
eat_model <- EAT(data, x, y)
pred_eat <- predict(eat_model, data, x)
mse(data$yD, pred_eat)

# ==
# PLOT
# ==
# plot
dataplt <- data.frame(x = data$x1,
                      y = data$y,
                      yD = data$yD,
                      pred_boost = pred_boost$y_pred,
                      pred_boost_smooth = pred_boost_smooth$y_pred,
                      pred_dea = pred_DEA,
                      pred_fdh = pred_FDH,
                      pred_eat = pred_eat$y_pred)

ggplot(dataplt) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = yD, colour = 'yD')) +
  geom_line(aes(x = x, y = pred_boost, colour = 'MARSBoost')) +
  geom_line(aes(x = x, y = pred_boost_smooth, colour = 'Smooth')) +
  geom_line(aes(x = x, y = pred_dea, colour = 'DEA')) +
  geom_line(aes(x = x, y = pred_fdh, colour = 'FDH')) +
  geom_line(aes(x = x, y = pred_eat, colour = 'EAT'))
