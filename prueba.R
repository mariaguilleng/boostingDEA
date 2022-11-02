library(ggplot2)
library(compare)
library(eat)
library(MLmetrics)

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
                         num_iterarions, learning_rate, 5)
# without smoothin
pred_boost <- predict(model_boost, data, x, class = 1)
sum(pred_boost == model_boost$prediction) == N

# with smoothing
pred_boost_smooth <- predict(model_boost, data, x, class = 2)
sum(pred_boost_smooth == model_boost$prediction.smooth) == N

MSE(pred_boost$y_pred, data$yD)
MSE(pred_boost_smooth$y_pred, data$yD)

# tuning
selected <- sample(1:N, N * 0.8) # Training indexes
training <- data[selected, ] # Training set
test <- data[- selected, ] # Test set
result <- bestMARSBoost(training, test, x, y,
                        num.iterations = c(8,9,10,11,12),
                        learning.rate = c(0.4, 0.5, 0.6),
                        num.terms = c(6,7,8,9))
result[1,]
model_boost_tuned <- MARSBoost(data, x, y,
                               num.iterations = result[1, "num.iterations"],
                               learning.rate = result[1, "learning.rate"],
                               num.terms = result[1, "num.terms"])
# without smoothin
pred_boost <- predict(model_boost_tuned, data, x, class = 1)
sum(pred_boost == model_boost_tuned$prediction) == N

# with smoothing
pred_boost_smooth <- predict(model_boost_tuned, data, x, class = 2)
sum(pred_boost_smooth == model_boost_tuned$prediction.smooth) == N

MSE(pred_boost$y_pred, data$yD)
MSE(pred_boost_smooth$y_pred, data$yD)


# ==
# DEA
# ==
DEA_model <- DEA(data,x,y)
pred_DEA <- predict(DEA_model, data, x, y)
compare(pred_DEA,DEA_model$pred)
MSE(pred_DEA, data$yD)



# ==
# FDH
# ==
FDH_model <- FDH(data,x,y)
pred_FDH <- predict(FDH_model, data, x, y)
compare(pred_FDH,FDH_model$pred)
MSE(pred_FDH, data$yD)


# ==
# EATBoost
# ==
eatboost_model <- EATBoost(data, x, y, num_iterarions, nterms, learning_rate)
pred_eatboost <- predict(eatboost_model, data, x)
sum(eatboost_model$prediction == pred_eatboost) == N
sum(eatboost_model$prediction == eatboost_model$f0) == N


# =========
# EAT
# =========
eat_model <- EAT(data, x, y)
pred_eat <- predict(eat_model, data, x)
MSE(pred_eat$y_pred, data$yD)


# ==
# PLOT
# ==
dataplt <- data.frame(x = data$x1,
                      y = data$y,
                      yD = data$yD,
                      pred_boost = pred_boost$y_pred,
                      pred_boost_smooth = pred_boost_smooth$y_pred,
                      pred_dea = pred_DEA,
                      pred_fdh = pred_FDH,
                      pred_eat = pred_eat$y_pred,
                      pred_eatboost = pred_eatboost$y_pred)

ggplot(dataplt) +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = yD, colour = 'yD')) +
  geom_line(aes(x = x, y = pred_boost, colour = 'MARSBoost')) +
  geom_line(aes(x = x, y = pred_boost_smooth, colour = 'Smooth')) +
  geom_line(aes(x = x, y = pred_dea, colour = 'DEA')) +
  geom_line(aes(x = x, y = pred_fdh, colour = 'FDH')) +
  geom_line(aes(x = x, y = pred_eat, colour = 'EAT')) +
  geom_line(aes(x = x, y = pred_eatboost, colour = 'EATBoost'))


data("PISAindex")
PISAindex
